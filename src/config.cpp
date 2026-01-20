#include "config.hpp"

namespace {

/**
 * Generic template function for parsing per-subsystem settings.
 * Handles two formats:
 * 1. A single table that applies to all subsystems
 * 2. An array of tables with 'subsystem' field for per-subsystem specification
 *
 * @tparam SettingsType The type of settings to parse (e.g., ControlParameterizationSettings)
 * @tparam ParseFunc The type of the parsing function
 * @param toml The TOML table containing the configuration
 * @param key The key name in the TOML table (e.g., "control_parameterization")
 * @param num_subsystems The total number of subsystems
 * @param default_settings The default settings to use for all subsystems initially
 * @param parse_func Function to parse a single table into SettingsType
 * @param logger Logger for error reporting
 * @return Vector of settings, one per subsystem
 */
template<typename SettingsType, typename ParseFunc>
std::vector<SettingsType> parsePerSubsystemSettings(const toml::table& toml, const std::string& key, size_t num_subsystems, const SettingsType& default_settings, ParseFunc parse_func, const MPILogger& logger) {

  // Case 1: Single table applies to all subsystems
  if (toml[key].is_table()) {
    auto* settings_table = toml[key].as_table();
    // Warn if 'subsystem' is specified in a single table - it will be ignored
    if (settings_table->contains("subsystem")) {
      logger.log("# Warning: '" + key + "' is a single table, so 'subsystem' field is ignored. Use array of tables to specify per-subsystem settings.\n");
    }
    SettingsType parsed_settings = parse_func(*settings_table);
    return std::vector<SettingsType>(num_subsystems, parsed_settings);
  }
  // Case 2: Array of tables for per-subsystem specification
  else if (toml[key].is_array()) {
    std::vector<SettingsType> settings;
    auto* settings_array = toml[key].as_array();
    for (const auto& elem : *settings_array) {
      if (!elem.is_table()) {
        logger.exitWithError(key + " array elements must be tables");
      }
      auto* elem_table = elem.as_table();

      // Get the subsystem index, or pair of indices for coupling parameters
      const std::string subsystem_key = "subsystem";
      if (!elem_table->contains(subsystem_key)) {
        logger.exitWithError(key + " array element must have 'subsystem' field");
      }
      // For coupling parameters, subsystem is a pair of indices, otherwise its a single index
      size_t index;
      if (elem_table->get(subsystem_key)->is_array()) {
        // Coupling parameter case: subsystem is an array of two indices
        auto subsys_array = validators::vectorField<size_t>(*elem_table, subsystem_key).hasLength(2).value();
        size_t i = subsys_array[0];;
        size_t j = subsys_array[1];
        if (i >= num_subsystems || j >= num_subsystems) {
          throw validators::ValidationError(key, "subsystem index out of range for key '" + key + "'");
        }
        // Compute unique index for pair (i,j) with i<j: Convert to linear index: (0,1), (0,2), ..., (0,n-1), (1,2), ..., (1,n-1), ..., (n-2,n-1)
        if (i > j) std::swap(i, j);
        index = i * (num_subsystems-1) - i * (i - 1) / 2 + (j - i - 1);
        size_t num_pairs = (num_subsystems * (num_subsystems - 1)) / 2;
        settings.resize(num_pairs, default_settings);
      } else if (elem_table->get(subsystem_key)->is_value()) {
        // Single subsystem index case
        index = validators::field<size_t>(*elem_table, subsystem_key).greaterThanEqual(0).lessThan(num_subsystems).value();
        settings.resize(num_subsystems, default_settings);
      } else {
        logger.exitWithError("subsystem field must be an integer index or an array of two indices");
      }
      // Parse the settings for this entry
      settings[index] = parse_func(*elem_table);
    }
    return settings;
  } else {
    logger.exitWithError(key + " must be a table (applies to all) or an array of tables (per-subsystem specification)");
  }
  return std::vector<SettingsType>();
}

} // namespace


Config::Config(const MPILogger& logger, const toml::table& toml) : logger(logger) {
  try {
    // Get section tables - only [system] is required
    const auto* system_table = toml["system"].as_table();
    if (system_table == nullptr) {
      logger.exitWithError("[system] table is required");
    }

    // Other tables are optional - use empty table as fallback
    const auto control_table = toml["control"].is_table() ? *toml["control"].as_table() : toml::table{};
    const auto optimization_table = toml["optimization"].is_table() ? *toml["optimization"].as_table() : toml::table{};
    const auto output_table = toml["output"].is_table() ? *toml["output"].as_table() : toml::table{};
    const auto solver_table = toml["solver"].is_table() ? *toml["solver"].as_table() : toml::table{};

    // Parse system options from [system] table

    nlevels = validators::vectorField<size_t>(*system_table, "nlevels").minLength(1).positive().value();
    size_t num_osc = nlevels.size();

    nessential = validators::scalarOrVectorOr<size_t>(*system_table, "nessential", num_osc, nlevels);

    ntime = validators::field<size_t>(*system_table, "ntime").positive().value();

    dt = validators::field<double>(*system_table, "dt").positive().value();

    transfreq = validators::scalarOrVector<double>(*system_table, "transfreq", num_osc);

    selfkerr = validators::scalarOrVectorOr<double>(*system_table, "selfkerr", num_osc, std::vector<double>(num_osc, ConfigDefaults::SELFKERR));

    // Parse crosskerr and Jkl coupling: either one value (all-to-all coupling) or array of tables with 'subsystem = [i,j]' field for i-j coupling)
    size_t num_pairs = (num_osc - 1) * num_osc / 2;
    crosskerr.assign(num_pairs, ConfigDefaults::CROSSKERR);
    Jkl.assign(num_pairs, ConfigDefaults::JKL);
    // Overwrite for crosskerr
    if (system_table->contains("crosskerr")) {
      if ((*system_table)["crosskerr"].is_value()) {
        double single_val = validators::field<double>(*system_table, "crosskerr").value();
        crosskerr.assign(num_pairs, single_val);
      } else {
      auto parseFunc = [](const toml::table& t) { return validators::field<double>(t, "value").value(); };
      crosskerr = parsePerSubsystemSettings<double>(*system_table, "crosskerr", num_osc, ConfigDefaults::CROSSKERR, parseFunc, logger);
      }
    }
    // Overwrite for Jkl
    if (system_table->contains("Jkl")) {
      if ((*system_table)["Jkl"].is_value()) {
        double single_val = validators::field<double>(*system_table, "Jkl").value();
        Jkl.assign(num_pairs, single_val);
      } else {
      auto parseFunc = [](const toml::table& t) { return validators::field<double>(t, "value").value(); };
      Jkl = parsePerSubsystemSettings<double>(*system_table, "Jkl", num_osc, ConfigDefaults::JKL, parseFunc, logger);
      }
    }

    rotfreq = validators::scalarOrVectorOr<double>(*system_table, "rotfreq", num_osc, std::vector<double>(num_osc, ConfigDefaults::ROTFREQ));

    hamiltonian_file_Hsys = validators::getOptional<std::string>((*system_table)["hamiltonian_file_Hsys"]);
    hamiltonian_file_Hc = validators::getOptional<std::string>((*system_table)["hamiltonian_file_Hc"]);

    // Parse decoherence setting
    decoherence_type = ConfigDefaults::DECOHERENCE_TYPE;
    decay_time = std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME);
    dephase_time = std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME);
    if (system_table->contains("decoherence")) {
      auto* decoherence_table = (*system_table)["decoherence"].as_table();
      if (!decoherence_table) {
        logger.exitWithError("decoherence must be a table");
      }
      auto type_str = validators::field<std::string>(*decoherence_table, "type").valueOr("none");
      decoherence_type = parseEnum(type_str, DECOHERENCE_TYPE_MAP, ConfigDefaults::DECOHERENCE_TYPE);
      decay_time = validators::scalarOrVectorOr<double>(*decoherence_table, "decay_time", num_osc, std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
      dephase_time = validators::scalarOrVectorOr<double>(*decoherence_table, "dephase_time", num_osc, std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
    }

    // Parse initial condition table
    auto init_cond_table = validators::getRequiredTable(*system_table, "initial_condition");
    auto type_opt = parseEnum(validators::field<std::string>(init_cond_table, "type").value(), INITCOND_TYPE_MAP);
    if (!type_opt.has_value()) {
      logger.exitWithError("initial condition type not found.");
    }
    initial_condition.type = type_opt.value();
    initial_condition.levels = validators::getOptionalVector<size_t>(init_cond_table["levels"]);
    initial_condition.filename = validators::getOptional<std::string>(init_cond_table["filename"]);
    initial_condition.subsystem= validators::getOptionalVector<size_t>(init_cond_table["subsystem"]);

    // Parse control options from [control] table
    control_zero_boundary_condition = control_table["zero_boundary_condition"].value_or(ConfigDefaults::CONTROL_ZERO_BOUNDARY_CONDITION);

    // Parse control parameterization, either table (applies to all) or array (per-oscillator)
    ControlParameterizationSettings default_param;
    control_parameterizations.assign(num_osc, default_param);
    if (control_table.contains("parameterization")) {
      auto parseParamFunc = [this](const toml::table& t) { return parseControlParameterizationSpecs(t); };
      control_parameterizations = parsePerSubsystemSettings<ControlParameterizationSettings>(control_table, "parameterization", num_osc, default_param, parseParamFunc, logger);
    }

    // Parse control initialization, either as a table (applies to all) or per-oscillator table
    ControlInitializationSettings default_init;
    control_initializations.assign(num_osc, default_init);
    if (control_table.contains("initialization")) {
      auto parseInitFunc = [this](const toml::table& t) { return parseControlInitializationSpecs(t); };
      control_initializations = parsePerSubsystemSettings<ControlInitializationSettings>(control_table, "initialization", num_osc, default_init, parseInitFunc, logger);
    }

    // Parse optional control bounds: either single value or per-oscillator array
    control_amplitude_bounds = validators::scalarOrVectorOr<double>(control_table, "amplitude_bound", num_osc, std::vector<double>(num_osc, ConfigDefaults::CONTROL_AMPLITUDE_BOUND));

    // Parse carrier frequencies: either one vector (applies to all oscillators) or per-oscillator array of tables
    std::vector<double> default_carrier_freq = {ConfigDefaults::CARRIER_FREQ};
    carrier_frequencies.assign(num_osc, default_carrier_freq);
    if (control_table.contains("carrier_frequency")) {
      // Check if carrier_frequency is a direct array of values (shorthand for applying to all)
      auto* carrier_freq_array = control_table["carrier_frequency"].as_array();
      if (carrier_freq_array && !carrier_freq_array->empty() && !carrier_freq_array->front().is_table()) {
        // Direct array format: carrier_frequency = [1.0, 2.0]
        auto values = validators::vectorField<double>(control_table, "carrier_frequency").value();
        carrier_frequencies.assign(num_osc, values);
      } else {
        // Table or array of tables format
        auto parseFunc = [](const toml::table& t) { return validators::vectorField<double>(t, "value").value(); };
        carrier_frequencies = parsePerSubsystemSettings<std::vector<double>>(control_table, "carrier_frequency", num_osc, default_carrier_freq, parseFunc, logger);
      }
    }

    // Parse optimization options from [optimization] table
    optim_target = parseOptimTarget(optimization_table, num_osc);

    optim_objective = parseEnum(optimization_table["objective"].value<std::string>(), OBJECTIVE_TYPE_MAP, ConfigDefaults::OPTIM_OBJECTIVE);

    // Parse optional weights
    optim_weights = validators::vectorField<double>(optimization_table, "weights").valueOr({ConfigDefaults::OPTIM_WEIGHT});

    // Parse optional optimization tolerances
    if (!optimization_table.contains("tolerance")) {
      optim_tol_grad_abs = ConfigDefaults::OPTIM_TOL_GRAD_ABS;
      optim_tol_grad_rel = ConfigDefaults::OPTIM_TOL_GRAD_REL;
      optim_tol_finalcost = ConfigDefaults::OPTIM_TOL_FINALCOST;
      optim_tol_infidelity = ConfigDefaults::OPTIM_TOL_INFIDELITY;
    } else {
      // Parse tolerance table
      auto* tol_table = optimization_table["tolerance"].as_table();
      if (!tol_table) {
        logger.exitWithError("tolerance must be a table");
      }
      optim_tol_grad_abs = validators::field<double>(*tol_table, "grad_abs").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_ABS);
      optim_tol_grad_rel = validators::field<double>(*tol_table, "grad_rel").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_REL);
      optim_tol_finalcost = validators::field<double>(*tol_table, "final_cost").positive().valueOr(ConfigDefaults::OPTIM_TOL_FINALCOST);
      optim_tol_infidelity = validators::field<double>(*tol_table, "infidelity").positive().valueOr(ConfigDefaults::OPTIM_TOL_INFIDELITY);
    }

    optim_maxiter = validators::field<size_t>(optimization_table, "maxiter").greaterThanEqual(0).valueOr(ConfigDefaults::OPTIM_MAXITER);

    // Parse tikhonov inline table
    if (!optimization_table.contains("tikhonov")) {
      optim_tikhonov_coeff = ConfigDefaults::OPTIM_TIKHONOV_COEFF;
      optim_tikhonov_use_x0 = ConfigDefaults::OPTIM_TIKHONOV_USE_X0;
    } else {
      auto regul_table = optimization_table["tikhonov"].as_table();
      if (!regul_table) {
        logger.exitWithError("tikhonov must be a table");
      }
      optim_tikhonov_coeff= validators::field<double>(*regul_table, "coeff").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_TIKHONOV_COEFF);
      optim_tikhonov_use_x0 = validators::field<bool>(*regul_table, "use_x0").valueOr(ConfigDefaults::OPTIM_TIKHONOV_USE_X0);
    }

    // Parse penalty table
    if (!optimization_table.contains("penalty")) {
      // No penalty table is specified, use defaults
      optim_penalty_leakage = ConfigDefaults::OPTIM_PENALTY_LEAKAGE;
      optim_penalty_weightedcost = ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST;
      optim_penalty_weightedcost_width = ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH;
      optim_penalty_dpdm = ConfigDefaults::OPTIM_PENALTY_DPDM;
      optim_penalty_energy = ConfigDefaults::OPTIM_PENALTY_ENERGY;
      optim_penalty_variation = ConfigDefaults::OPTIM_PENALTY_VARIATION;
    } else {
      // Parse penalty table
      auto penalty_table = optimization_table["penalty"].as_table();
      if (!penalty_table) {
        logger.exitWithError("penalty must be a table");
      }
      optim_penalty_leakage = validators::field<double>(*penalty_table, "leakage").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
      optim_penalty_weightedcost = validators::field<double>(*penalty_table, "weightedcost").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
      optim_penalty_weightedcost_width = validators::field<double>(*penalty_table, "weightedcost_width").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
      optim_penalty_dpdm = validators::field<double>(*penalty_table, "dpdm").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_DPDM);
      optim_penalty_energy = validators::field<double>(*penalty_table, "energy").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_ENERGY);
      optim_penalty_variation = validators::field<double>(*penalty_table, "variation").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_VARIATION);
    }

    // Parse output options from [output] table
    output_directory = output_table["directory"].value_or(ConfigDefaults::OUTPUT_DIRECTORY);

    // Parse observables as an array of strings (defaults to empty array)
    output_observables.clear();
    if (auto output_observables_array = output_table["observables"].as_array()) {
      for (auto&& elem : *output_observables_array) {
        if (auto str = elem.value<std::string>()) {
          auto enum_val = parseEnum(*str, OUTPUT_TYPE_MAP);
          if (!enum_val.has_value()) {
            logger.exitWithError("Unknown output type: " + *str);
          }
          output_observables.push_back(enum_val.value());
        } else {
          logger.exitWithError("output type array must contain strings");
        }
      }
    }

    output_timestep_stride = validators::field<size_t>(output_table, "timestep_stride").greaterThanEqual(0).valueOr(ConfigDefaults::OUTPUT_TIMESTEP_STRIDE);

    output_optimization_stride = validators::field<size_t>(output_table, "optimization_stride").greaterThanEqual(0).valueOr(ConfigDefaults::OUTPUT_OPTIMIZATION_STRIDE);

    // Parse solver options from [solver] table
    runtype = parseEnum(solver_table["runtype"].value<std::string>(), RUN_TYPE_MAP, ConfigDefaults::RUNTYPE);

    usematfree = solver_table["usematfree"].value_or(ConfigDefaults::USEMATFREE);

    // Parse linearsolver as an inline table
    if (!solver_table.contains("linearsolver")) {
      // No linearsolver table specified, use defaults
      linearsolver_type = ConfigDefaults::LINEARSOLVER_TYPE;
      linearsolver_maxiter = ConfigDefaults::LINEARSOLVER_MAXITER;
    } else {
      auto* linearsolver_table_inner = solver_table["linearsolver"].as_table();
      if (!linearsolver_table_inner) {
        logger.exitWithError("linearsolver must be a table");
      }
      linearsolver_type = parseEnum(validators::field<std::string>(*linearsolver_table_inner, "type").value(), LINEAR_SOLVER_TYPE_MAP, ConfigDefaults::LINEARSOLVER_TYPE);
      linearsolver_maxiter = validators::field<size_t>(*linearsolver_table_inner, "maxiter").positive().valueOr(ConfigDefaults::LINEARSOLVER_MAXITER);
    }

    timestepper_type = parseEnum(solver_table["timestepper"].value<std::string>(), TIME_STEPPER_TYPE_MAP, ConfigDefaults::TIMESTEPPER_TYPE);

    int rand_seed_ = solver_table["rand_seed"].value_or(ConfigDefaults::RAND_SEED);
    setRandSeed(rand_seed_);

  } catch (const validators::ValidationError& e) {
    logger.exitWithError(std::string(e.what()));
  }

  // Finalize interdependent settings, then validate
  finalize();
  validate();
}

Config::Config(const MPILogger& logger, const ParsedConfigData& settings) : logger(logger) {

  if (!settings.nlevels.has_value()) {
    logger.exitWithError("nlevels cannot be empty");
  }
  nlevels = settings.nlevels.value();
  size_t num_osc = nlevels.size();

  nessential = settings.nessential.value_or(nlevels);
  copyLast(nessential, num_osc);

  if (!settings.ntime.has_value()) {
    logger.exitWithError("ntime cannot be empty");
  }
  ntime = settings.ntime.value();
  if (ntime <= 0) {
    logger.exitWithError("ntime must be positive, got " + std::to_string(ntime));
  }

  if (!settings.dt.has_value()) {
    logger.exitWithError("dt cannot be empty");
  }
  dt = settings.dt.value();
  if (dt <= 0) {
    logger.exitWithError("dt must be positive, got " + std::to_string(dt));
  }

  if (!settings.transfreq.has_value()) {
    logger.exitWithError("transfreq cannot be empty");
  }

  transfreq = settings.transfreq.value();
  copyLast(transfreq, num_osc);

  selfkerr = settings.selfkerr.value_or(std::vector<double>(num_osc, ConfigDefaults::SELFKERR));
  copyLast(selfkerr, num_osc);

  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  crosskerr = settings.crosskerr.value_or(std::vector<double>(num_pairs_osc, ConfigDefaults::CROSSKERR));
  copyLast(crosskerr, num_pairs_osc);

  Jkl = settings.Jkl.value_or(std::vector<double>(num_pairs_osc, ConfigDefaults::JKL));
  copyLast(Jkl, num_pairs_osc);

  rotfreq = settings.rotfreq.value_or(std::vector<double>(num_osc, ConfigDefaults::ROTFREQ));
  copyLast(rotfreq, num_osc);

  hamiltonian_file_Hsys = settings.hamiltonian_file_Hsys;
  hamiltonian_file_Hc = settings.hamiltonian_file_Hc;

  decoherence_type = settings.decoherence_type.value_or(ConfigDefaults::DECOHERENCE_TYPE);

  decay_time = settings.decay_time.value_or(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
  copyLast(decay_time, num_osc);

  dephase_time = settings.dephase_time.value_or(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
  copyLast(dephase_time, num_osc);

  if (!settings.initialcondition.has_value()) {
    logger.exitWithError("initialcondition cannot be empty");
  }
  initial_condition.type = settings.initialcondition.value().type;
  initial_condition.filename = settings.initialcondition.value().filename;
  initial_condition.levels = settings.initialcondition.value().levels;
  initial_condition.subsystem = settings.initialcondition.value().subsystem;

  // Control and optimization parameters
  control_zero_boundary_condition = settings.control_zero_boundary_condition.value_or(ConfigDefaults::CONTROL_ZERO_BOUNDARY_CONDITION);

  control_parameterizations = parseControlParameterizationsCfg(settings.indexed_control_parameterizations);

  // Control initialization
  if (settings.indexed_control_init.has_value()) {
    auto init_map = settings.indexed_control_init.value();
    // First check for global file initialization and populate to all oscillators if present
    if (init_map.find(0) != init_map.end() && init_map[0].filename.has_value()) {
      control_initializations.resize(num_osc);
      std::string control_initialization_file = init_map[0].filename.value();
      for (size_t i = 0; i < num_osc; i++) {
        control_initializations[i] = ControlInitializationSettings{ControlInitializationType::FILE, std::nullopt, std::nullopt, control_initialization_file};
      }
    } else {
      control_initializations = parseControlInitializationsCfg(settings.indexed_control_init);
    }
  } else {
    // Initialize with defaults when no control initialization is provided
    control_initializations.resize(num_osc);
    for (size_t i = 0; i < num_osc; i++) {
      control_initializations[i] = ControlInitializationSettings{
        ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt};
    }
  }

  // Parse control_amplitude_bounds from CFG format (returns vector of vectors, but we need vector)
  auto control_amplitude_bounds_cfg = parseOscillatorSettingsCfg<double>(settings.indexed_control_amplitude_bounds, control_parameterizations.size(), {ConfigDefaults::CONTROL_AMPLITUDE_BOUND});
  control_amplitude_bounds.resize(control_amplitude_bounds_cfg.size());
  for (size_t i = 0; i < control_amplitude_bounds_cfg.size(); ++i) {
    control_amplitude_bounds[i] = control_amplitude_bounds_cfg[i].empty() ? ConfigDefaults::CONTROL_AMPLITUDE_BOUND : control_amplitude_bounds_cfg[i][0];
  }

  carrier_frequencies.resize(num_osc);
  carrier_frequencies = parseOscillatorSettingsCfg<double>(settings.indexed_carrier_frequencies, num_osc, {ConfigDefaults::CARRIER_FREQ});

  if (settings.optim_target.has_value()) {
    optim_target = settings.optim_target.value();
    optim_target.gate_rot_freq = settings.gate_rot_freq.value_or(std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));
  } else {
    // No optim_target specified, use default (no target)
    OptimTargetSettings default_target;
    optim_target = default_target;
  }

  optim_objective = settings.optim_objective.value_or(ConfigDefaults::OPTIM_OBJECTIVE);

  optim_weights = settings.optim_weights.value_or(std::vector<double>{ConfigDefaults::OPTIM_WEIGHT});

  optim_tol_grad_abs = settings.optim_tol_grad_abs.value_or(ConfigDefaults::OPTIM_TOL_GRAD_ABS);
  optim_tol_grad_rel = settings.optim_tol_grad_rel.value_or(ConfigDefaults::OPTIM_TOL_GRAD_REL);
  optim_tol_finalcost = settings.optim_tol_finalcost.value_or(ConfigDefaults::OPTIM_TOL_FINALCOST);
  optim_tol_infidelity = settings.optim_tol_infidelity.value_or(ConfigDefaults::OPTIM_TOL_INFIDELITY);
  optim_maxiter = settings.optim_maxiter.value_or(ConfigDefaults::OPTIM_MAXITER);

  optim_tikhonov_coeff = settings.optim_regul.value_or(ConfigDefaults::OPTIM_TIKHONOV_COEFF);
  optim_tikhonov_use_x0 = settings.optim_regul_tik0.value_or(ConfigDefaults::OPTIM_TIKHONOV_USE_X0);
  if (settings.optim_regul_interpolate.has_value()) {
    // Handle deprecated optim_regul_interpolate logic
    optim_tikhonov_use_x0 = settings.optim_regul_interpolate.value();
    logger.log("# Warning: 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
  }

  optim_penalty_leakage = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
  optim_penalty_weightedcost = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
  optim_penalty_weightedcost_width = settings.optim_penalty_param.value_or(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
  optim_penalty_dpdm = settings.optim_penalty_dpdm.value_or(ConfigDefaults::OPTIM_PENALTY_DPDM);
  optim_penalty_energy = settings.optim_penalty_energy.value_or(ConfigDefaults::OPTIM_PENALTY_ENERGY);
  optim_penalty_variation = settings.optim_penalty_variation.value_or(ConfigDefaults::OPTIM_PENALTY_VARIATION);

  // Output parameters
  output_directory = settings.datadir.value_or(ConfigDefaults::OUTPUT_DIRECTORY);

  // Convert old per-oscillator output to global output_observables (apply to all oscillators)
  auto indexed_output_vec = parseOscillatorSettingsCfg<OutputType>(settings.indexed_output, num_osc);
  output_observables.clear();
  // Collect unique output types from all oscillators
  std::set<OutputType> unique_types;
  for (const auto& osc_output : indexed_output_vec) {
    for (const auto& type : osc_output) {
      unique_types.insert(type);
    }
  }
  // Convert set to vector
  output_observables.assign(unique_types.begin(), unique_types.end());

  output_timestep_stride = settings.output_timestep_stride.value_or(ConfigDefaults::OUTPUT_TIMESTEP_STRIDE);
  output_optimization_stride = settings.output_optimization_stride.value_or(ConfigDefaults::OUTPUT_OPTIMIZATION_STRIDE);
  runtype = settings.runtype.value_or(ConfigDefaults::RUNTYPE);
  usematfree = settings.usematfree.value_or(ConfigDefaults::USEMATFREE);
  linearsolver_type = settings.linearsolver_type.value_or(ConfigDefaults::LINEARSOLVER_TYPE);
  linearsolver_maxiter = settings.linearsolver_maxiter.value_or(ConfigDefaults::LINEARSOLVER_MAXITER);
  timestepper_type = settings.timestepper_type.value_or(ConfigDefaults::TIMESTEPPER_TYPE);
  setRandSeed(settings.rand_seed.value_or(ConfigDefaults::RAND_SEED));

  // Finalize interdependent settings, then validate
  finalize();
  validate();
}

Config::~Config() {}

Config Config::fromFile(const std::string& filename, const MPILogger& logger) {
  if (hasSuffix(filename, ".toml")) {
    return Config::fromToml(filename, logger);
  } else {
    // TODO cfg: delete this when .cfg format is removed.
    logger.log(
        "# Warning: Config file does not have .toml extension. "
        "The deprecated .cfg format will be removed in future versions.\n");
    return Config::fromCfg(filename, logger);
  }
}

Config Config::fromToml(const std::string& filename, const MPILogger& logger) {
  toml::table toml = toml::parse_file(filename);
  return Config(logger, toml);
}

Config Config::fromTomlString(const std::string& toml_content, const MPILogger& logger) {
  toml::table toml = toml::parse(toml_content);
  return Config(logger, toml);
}

Config Config::fromCfg(const std::string& filename, const MPILogger& logger) {
  CfgParser parser(logger);
  ParsedConfigData settings = parser.parseFile(filename);
  return Config(logger, settings);
}

Config Config::fromCfgString(const std::string& cfg_content, const MPILogger& logger) {
  CfgParser parser(logger);
  ParsedConfigData settings = parser.parseString(cfg_content);
  return Config(logger, settings);
}

namespace {

std::string formatDouble(double value) {
  std::ostringstream oss;
  oss << std::setprecision(std::numeric_limits<double>::digits10);
  oss << value;
  // Format e.g. 0 as 0.0
  std::string str = oss.str();
  if (str.find('.') == std::string::npos && str.find('e') == std::string::npos) {
    str += ".0";
  }
  return str;
}

std::string toStringCoupling(const std::vector<double>& couplings, size_t num_osc) {
  if (couplings.empty()) return "[]";

  // If all couplings are the same, print single value
  bool all_equal = std::adjacent_find(couplings.begin(), couplings.end(), std::not_equal_to<double>{}) == couplings.end();
  if (all_equal) {
    return formatDouble(couplings[0]);
  }

  // Collect non-zero couplings with their pair indices
  std::vector<std::pair<std::pair<size_t, size_t>, double>> nonzero_couplings;
  size_t pair_idx = 0;
  for (size_t i = 0; i < num_osc - 1; i++) {
    for (size_t j = i + 1; j < num_osc; j++) {
      if (pair_idx < couplings.size() && couplings[pair_idx] != 0.0) {
        nonzero_couplings.push_back({{i, j}, couplings[pair_idx]});
      }
      pair_idx++;
    }
  }

  // Build TOML table format
  std::string result = "[\n";
  for (size_t i = 0; i < nonzero_couplings.size(); ++i) {
    auto [pair, value] = nonzero_couplings[i];
    auto [first, second] = pair;
    result += " { subsystem = [" + std::to_string(first) + "," + std::to_string(second) + "], value = " + formatDouble(value) + "}";
    if (i < nonzero_couplings.size() - 1) {
      result += ", ";
    }
    result += "\n";
  }
  result += "]";
  return result;
}

template <typename T>
std::string printVector(const std::vector<T>& vec) {
  if (vec.empty()) return "[]";

  std::string result = "[";
  if constexpr (std::is_floating_point_v<T>) {
    result += formatDouble(vec[0]);
    for (size_t i = 1; i < vec.size(); ++i) {
      result += ", " + formatDouble(vec[i]);
    }
  } else {
    result += std::to_string(vec[0]);
    for (size_t i = 1; i < vec.size(); ++i) {
      result += ", " + std::to_string(vec[i]);
    }
  }
  result += "]";
  return result;
}


std::string toString(const InitialConditionSettings& initial_condition) {
  auto type_str = "type = \"" + enumToString(initial_condition.type, INITCOND_TYPE_MAP) + "\"";
  switch (initial_condition.type) {
    case InitialConditionType::FROMFILE:
      return "{" + type_str + ", filename = \"" + initial_condition.filename.value() + "\"}";
    case InitialConditionType::PRODUCT_STATE: {
      std::string out = "{" + type_str + ", levels = ";
      out += printVector(initial_condition.levels.value());
      out += "}";
      return out;
    }
    case InitialConditionType::ENSEMBLE: {
      std::string out = "{" + type_str + ", subsystem = ";
      out += printVector(initial_condition.subsystem.value());
      out += "}";
      return out;
    }
    case InitialConditionType::DIAGONAL: {
      std::string out = "{" + type_str + ", subsystem = ";
      out += printVector(initial_condition.subsystem.value());
      out += "}";
      return out;
    }
    case InitialConditionType::BASIS: {
      std::string out = "{" + type_str + ", subsystem = ";
      out += printVector(initial_condition.subsystem.value());
      out += "}";
      return out;
    }
    case InitialConditionType::THREESTATES:
    case InitialConditionType::NPLUSONE:
    case InitialConditionType::PERFORMANCE:
      return "{" + type_str + "}";
  }
  return "unknown";
}

std::string toString(const OptimTargetSettings& optim_target) {
  auto type_str = "type = \"" + enumToString(optim_target.type, TARGET_TYPE_MAP) + "\"";
  switch (optim_target.type) {
    case TargetType::GATE: {
      std::string out = "{" + type_str;
      if (optim_target.gate_type.has_value()) {
        out += ", gate_type = \"" + enumToString(optim_target.gate_type.value(), GATE_TYPE_MAP) + "\"";
      }
      if (optim_target.filename.has_value() && !optim_target.filename.value().empty()) {
        out += ", filename = \"" + optim_target.filename.value() + "\"";
      }
      if (optim_target.gate_rot_freq.has_value()) {
        out += ", gate_rot_freq = " + printVector(optim_target.gate_rot_freq.value());
      }
      out += "}";
      return out;
    }
    case TargetType::STATE: {
      std::string out = "{" + type_str;
      if (optim_target.levels.has_value()) {
        out += ", levels = " + printVector(optim_target.levels.value());
      }
      if (optim_target.filename.has_value() && !optim_target.filename.value().empty()) {
        out += ", filename = \"" + optim_target.filename.value() + "\"";
      }
      out += "}";
      return out;
    }
    case TargetType::NONE:
      return "{" + type_str + "}";
  }
  return "unknown";
}

// Template helper for toString functions that output either a single item or an array with per-item overrides
template <typename T, typename PrintFunc, typename CompareFunc>
std::string toStringWithOptionalPerSubsystem(const std::vector<T>& items, PrintFunc printItems, CompareFunc areEqual) {

  if (items.empty()) return "[]";

  // Check if all items are the same
  bool all_equal = std::adjacent_find(items.begin(), items.end(), [&areEqual](const auto& a, const auto& b) { return !areEqual(a, b); }) == items.end();

  if (all_equal) {
    std::string out = printItems(items.front());
    // If items are not wrapped in either {...} or [...], add {} here.
    if (out.front() != '{' && out.front() != '[') {
      out = "{" + out + "}";
    }
    return out;
  } else {
    // Output as array with per-subsystem overrides
    std::string out = "[\n";
    for (size_t i = 0; i < items.size(); ++i) {
      out += "  { subsystem = " + std::to_string(i) + ", " + printItems(items[i])+ "}";
      if (i < items.size() - 1) {
        out += ",";
      }
      out += "\n";
    }
    out += "]";
    return out;
  }
}

std::string toString(const std::vector<ControlInitializationSettings>& control_initializations) {
  // Helper function to print all items of a single ControlInitializationSettings
  auto printItems = [](const ControlInitializationSettings& init) {
    std::string out = "";
    out += "type = \"" + enumToString(init.type, CONTROL_INITIALIZATION_TYPE_MAP) + "\"";
    out += init.filename.has_value() ? ", filename = \"" + init.filename.value() + "\"" : "";
    out += init.amplitude.has_value() ? ", amplitude = " + formatDouble(init.amplitude.value()) : "";
    out += init.phase.has_value() ? ", phase = " + formatDouble(init.phase.value()) : "";
    return out;
  };

  // Helper function to compare two ControlInitializationSettings items
  auto areEqual = [](const ControlInitializationSettings& a, const ControlInitializationSettings& b) {
    return a.type == b.type && a.amplitude == b.amplitude && a.phase == b.phase;
  };

  return toStringWithOptionalPerSubsystem(control_initializations, printItems, areEqual);
}

std::string toString(const std::vector<ControlParameterizationSettings>& control_parameterizations) {
  // Helper function to print all items of a single ControlParameterizationSetting
  auto printItems = [](const ControlParameterizationSettings& param) {
    std::string out = "";
    out += "type = \"" + enumToString(param.type, CONTROL_TYPE_MAP) + "\"";
    out += param.nspline.has_value() ? ", num = " + std::to_string(param.nspline.value()) : "";
    out += param.tstart.has_value() ? ", tstart = " + formatDouble(param.tstart.value()) : "";
    out += param.tstop.has_value() ? ", tstop = " + formatDouble(param.tstop.value()) : "";
    out += param.scaling.has_value() ? ", scaling = " + formatDouble(param.scaling.value()) : "";
    return out;
  };

  // Helper function to compare two ControlParameterizationSettings items
  auto areEqual = [](const ControlParameterizationSettings& a, const ControlParameterizationSettings& b) {
    return a.type == b.type && a.nspline == b.nspline &&
           a.tstart == b.tstart && a.tstop == b.tstop && a.scaling == b.scaling;
  };

  return toStringWithOptionalPerSubsystem(control_parameterizations, printItems, areEqual);
}

std::string toString(const std::vector<std::vector<double>>& carrier_frequencies) {
  // Helper function to print all items of a single vector<double>
  auto printItems = [](const std::vector<double>& freqs) {
    return "value = " + printVector(freqs);
  };

  // Helper function to compare two vector<double> items
  auto areEqual = [](const std::vector<double>& a, const std::vector<double>& b) {
    return a == b;
  };

  return toStringWithOptionalPerSubsystem(carrier_frequencies, printItems, areEqual);
}

// Prints a single double value if all vector elements are equal, otherwise prints the vector
std::string toString(const std::vector<double>& vec) {
  if (vec.empty()) return "[]";

  bool all_equal = std::adjacent_find(vec.begin(), vec.end(), std::not_equal_to<double>{}) == vec.end();

  if (all_equal) {
    return formatDouble(vec[0]);
  }
  return printVector(vec);
}

} // namespace

// Print config as toml
// Decided to do this manually instead of with the tomlplusplus library so we could control the ordering and comments.
void Config::printConfig(std::stringstream& log) const {
  log << "[system]\n";

  // System parameters
  log << "nlevels = " << printVector(nlevels) << "\n";
  log << "nessential = " << printVector(nessential) << "\n";
  log << "ntime = " << ntime << "\n";
  log << "dt = " << dt << "\n";
  log << "transfreq = " << printVector(transfreq) << "\n";
  log << "selfkerr = " << printVector(selfkerr) << "\n";
  log << "crosskerr = " << toStringCoupling(crosskerr, nlevels.size()) << "\n";
  log << "Jkl = " << toStringCoupling(Jkl, nlevels.size()) << "\n";
  log << "rotfreq = " << printVector(rotfreq) << "\n";
  log << "decoherence = {\n";
  log << "  type = \"" << enumToString(decoherence_type, DECOHERENCE_TYPE_MAP) << "\",\n";
  log << "  decay_time = " << printVector(decay_time) << ",\n";
  log << "  dephase_time = " << printVector(dephase_time) << "\n";
  log << "}\n";
  log << "initial_condition = " << toString(initial_condition) << "\n";
  if (hamiltonian_file_Hsys.has_value()) {
    log << "hamiltonian_file_Hsys = \"" << hamiltonian_file_Hsys.value() << "\"\n";
  }
  if (hamiltonian_file_Hc.has_value()) {
    log << "hamiltonian_file_Hc = \"" << hamiltonian_file_Hc.value() << "\"\n";
  }

  log << "\n";
  log << "[control]\n";

  log << "parameterization = " << toString(control_parameterizations) << "\n";
  log << "carrier_frequency = " << toString(carrier_frequencies) << "\n";
  log << "initialization = " << toString(control_initializations) << "\n";
  log << "amplitude_bound = " << toString(control_amplitude_bounds) << "\n";
  log << "zero_boundary_condition = " << (control_zero_boundary_condition ? "true" : "false") << "\n";

  log << "\n";
  log << "[optimization]\n";

  log << "target = " << toString(optim_target) << "\n";
  log << "objective = \"" << enumToString(optim_objective, OBJECTIVE_TYPE_MAP) << "\"\n";
  bool uniform_weights = std::adjacent_find(optim_weights.begin(), optim_weights.end(), std::not_equal_to<double>{}) == optim_weights.end();
  if (!uniform_weights) {
    log << "weights = " << printVector(optim_weights) << "\n";
  }
  log << "tolerance = { grad_abs = " << optim_tol_grad_abs
      << ", grad_rel = " << optim_tol_grad_rel
      << ", final_cost = " << optim_tol_finalcost
      << ", infidelity = " << optim_tol_infidelity << " }\n";
  log << "maxiter = " << optim_maxiter << "\n";
  log << "tikhonov = { coeff = " << optim_tikhonov_coeff
      << ", use_x0 = " << (optim_tikhonov_use_x0 ? "true" : "false") << " }\n";
  log << "penalty = { leakage = " << optim_penalty_leakage
      << ", energy = " << optim_penalty_energy
      << ", dpdm = " << optim_penalty_dpdm
      << ", variation = " << optim_penalty_variation
      << ", weightedcost = " << optim_penalty_weightedcost
      << ", weightedcost_width = " << optim_penalty_weightedcost_width << " }\n";

  log << "\n";
  log << "[output]\n";

  log << "directory = \"" << output_directory << "\"\n";
  log << "observables = [";
  for (size_t j = 0; j < output_observables.size(); ++j) {
    log << "\"" << enumToString(output_observables[j], OUTPUT_TYPE_MAP) << "\"";
    if (j < output_observables.size() - 1) log << ", ";
  }
  log << "]\n";
  log << "timestep_stride = " << output_timestep_stride << "\n";
  log << "optimization_stride = " << output_optimization_stride << "\n";

  log << "\n";
  log << "[solver]\n";

  log << "runtype = \"" << enumToString(runtype, RUN_TYPE_MAP) << "\"\n";
  log << "usematfree = " << (usematfree ? "true" : "false") << "\n";
  log << "linearsolver = { type = \"" << enumToString(linearsolver_type, LINEAR_SOLVER_TYPE_MAP) << "\", maxiter = " << linearsolver_maxiter << " }\n";
  log << "timestepper = \"" << enumToString(timestepper_type, TIME_STEPPER_TYPE_MAP) << "\"\n";
  log << "rand_seed = " << rand_seed << "\n";
}

void Config::finalize() {
  // Hamiltonian file + matrix-free compatibility check
  if ((hamiltonian_file_Hsys.has_value() || hamiltonian_file_Hc.has_value()) && usematfree) {
    logger.log(
        "# Warning: Matrix-free solver cannot be used when Hamiltonian is read from file. Switching to sparse-matrix "
        "version.\n");
    usematfree = false;
  }

  if (usematfree && nlevels.size() > 5) {
    logger.log(
        "Warning: Matrix free solver is only implemented for systems with 2, 3, 4, or 5 oscillators."
        "Switching to sparse-matrix solver now.\n");
    usematfree = false;
  }

  // DIAGONAL and BASIS initial conditions in the Schroedinger case are the same. Overwrite it to DIAGONAL
  if (decoherence_type == DecoherenceType::NONE && initial_condition.type == InitialConditionType::BASIS) {
    initial_condition.type = InitialConditionType::DIAGONAL;
  }

  // For BASIS, ENSEMBLE, and DIAGONAL, or default to all oscillators IDs
  if (initial_condition.type == InitialConditionType::BASIS || initial_condition.type == InitialConditionType::ENSEMBLE || initial_condition.type == InitialConditionType::DIAGONAL) {
    if (!initial_condition.subsystem.has_value()) {
      initial_condition.subsystem = std::vector<size_t>(nlevels.size());
      for (size_t i = 0; i < nlevels.size(); i++) {
        initial_condition.subsystem->at(i) = i;
      }
    }
  }

  // Compute number of initial conditions
  n_initial_conditions = computeNumInitialConditions(initial_condition, nlevels, nessential, decoherence_type);

  // overwrite decay or dephase times with zeros, if the decoherence type is only one of them, or none.
  if (decoherence_type == DecoherenceType::DECAY) {
    std::fill(dephase_time.begin(), dephase_time.end(), 0);
  } else if (decoherence_type == DecoherenceType::DEPHASE) {
    std::fill(decay_time.begin(), decay_time.end(), 0);
  } else if (decoherence_type == DecoherenceType::NONE) {
    std::fill(decay_time.begin(), decay_time.end(), 0);
    std::fill(dephase_time.begin(), dephase_time.end(), 0);
  }

  // Scale optimization weights such that they sum up to one
  // If a single value was provided, replicate it for all initial conditions
  if (optim_weights.size() == 1) {
    // TODO remove this when removing cfg format
    copyLast(optim_weights, n_initial_conditions);
  } else if (optim_weights.size() != n_initial_conditions) {
    logger.exitWithError("optim_weights vector has length " + std::to_string(optim_weights.size()) + " but must have length " + std::to_string(n_initial_conditions) + " (number of initial conditions)");
  }
  // Scale the weights so that they sum up to 1
  double scaleweights = 0.0;
  for (size_t i = 0; i < optim_weights.size(); i++) scaleweights += optim_weights[i];
  for (size_t i = 0; i < optim_weights.size(); i++) optim_weights[i] = optim_weights[i] / scaleweights;

  // Set weightedcost width to zero if weightedcost penalty is zero
  if (optim_penalty_weightedcost == 0.0) {
    optim_penalty_weightedcost_width = 0.0;
  }

  // Set control variation penalty to zero if not using 2nd order Bspline parameterization
  for (size_t i = 0; i < control_parameterizations.size(); i++) {
    if (control_parameterizations[i].type != ControlType::BSPLINE0) {
      optim_penalty_variation = 0.0;
      break;
    }
  }

  // Unset control initialization phase paremeter, unless BSPLINEAMP parameterization is used
  for (size_t i = 0; i < control_initializations.size(); i++) {
    if (control_parameterizations[i].type != ControlType::BSPLINEAMP) {
      control_initializations[i].phase = std::nullopt;
    }
  }
}

void Config::validate() const {

  // Validate essential levels don't exceed total levels
  if (nessential.size() != nlevels.size()) {
    logger.exitWithError("nessential size must match nlevels size");
  }
  for (size_t i = 0; i < nlevels.size(); i++) {
    if (nessential[i] > nlevels[i]) {
      logger.exitWithError("nessential[" + std::to_string(i) + "] = " + std::to_string(nessential[i]) + " cannot exceed nlevels[" + std::to_string(i) + "] = " + std::to_string(nlevels[i]));
    }
  }

  /* Sanity check for Schrodinger solver initial conditions */
  if (decoherence_type == DecoherenceType::NONE) {
    if (initial_condition.type == InitialConditionType::ENSEMBLE ||
        initial_condition.type == InitialConditionType::THREESTATES ||
        initial_condition.type == InitialConditionType::NPLUSONE) {
      logger.exitWithError(
          "\n\n ERROR for initial condition setting: \n When running Schroedingers solver,"
          " the initial condition needs to be either 'state' or 'file' or 'diagonal' or "
          "'basis'."
          " Note that 'diagonal' and 'basis' in the Schroedinger case are the same (all unit vectors).\n\n");
    }
  }

  // Validate control bounds are positive
  for (size_t i = 0; i < control_amplitude_bounds.size(); i++) {
    if (control_amplitude_bounds[i] <= 0.0) {
      logger.exitWithError("control_amplitude_bounds[" + std::to_string(i) + "] must be positive");
    }
  }

  // Validate initial condition settings
  if (initial_condition.type == InitialConditionType::FROMFILE) {
    if (!initial_condition.filename.has_value()) {
      logger.exitWithError("initialcondition of type FROMFILE must have a filename");
    }
  }
  if (initial_condition.type == InitialConditionType::PRODUCT_STATE) {
    if (!initial_condition.levels.has_value()) {
      logger.exitWithError("initialcondition of type PRODUCT_STATE must have 'levels'");
    }
    if (initial_condition.levels->size() != nlevels.size()) {
      logger.exitWithError("initialcondition of type PRODUCT_STATE must have exactly " + std::to_string(nlevels.size()) + " parameters, got " + std::to_string(initial_condition.levels->size()));
    }
    for (size_t k = 0; k < initial_condition.levels->size(); k++) {
      if (initial_condition.levels->at(k) >= nlevels[k]) {
        logger.exitWithError("ERROR in config setting. The requested product state initialization " + std::to_string(initial_condition.levels->at(k)) + " exceeds the number of allowed levels for that oscillator (" + std::to_string(nlevels[k]) + ").\n");
      }
    }
  }
  if (initial_condition.type == InitialConditionType::BASIS ||
      initial_condition.type == InitialConditionType::DIAGONAL ||
      initial_condition.type == InitialConditionType::ENSEMBLE) {
    if (!initial_condition.subsystem.has_value()) {
      logger.exitWithError("initialcondition of type BASIS, DIAGONAL, or ENSEMBLE must have 'subsystem'");
    }
    if (initial_condition.subsystem->back() >= nlevels.size()) {
      logger.exitWithError("Last element in initialcondition params exceeds number of oscillators");
    }
    for (size_t i = 1; i < initial_condition.subsystem->size() - 1; i++) {
      if (initial_condition.subsystem->at(i) + 1 != initial_condition.subsystem->at(i + 1)) {
        logger.exitWithError("List of oscillators for ensemble initialization should be consecutive!\n");
      }
    }
  }
}

size_t Config::computeNumInitialConditions(InitialConditionSettings init_cond_settings, std::vector<size_t> nlevels, std::vector<size_t> nessential, DecoherenceType decoherence_type) const {
  size_t n_initial_conditions = 0;
  switch (init_cond_settings.type) {
    case InitialConditionType::FROMFILE:
    case InitialConditionType::PRODUCT_STATE:
    case InitialConditionType::PERFORMANCE:
    case InitialConditionType::ENSEMBLE:
      n_initial_conditions = 1;
      break;
    case InitialConditionType::THREESTATES:
      n_initial_conditions = 3;
      break;
    case InitialConditionType::NPLUSONE:
      // compute system dimension N
      n_initial_conditions = 1;
      for (size_t i = 0; i < nlevels.size(); i++) {
        n_initial_conditions *= nlevels[i];
      }
      n_initial_conditions += 1;
      break;
    case InitialConditionType::DIAGONAL:
      /* Compute ninit = dim(subsystem defined by list of oscil IDs) */
      if (!init_cond_settings.subsystem.has_value()) {
        logger.exitWithError("expected diagonal initial condition to have list of subsystems ");
      }
      n_initial_conditions = 1;
      for (size_t oscilID : init_cond_settings.subsystem.value()) {
        if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
      }
      break;
    case InitialConditionType::BASIS:
      /* Compute ninit = dim(subsystem defined by list of oscil IDs) */
      if (!init_cond_settings.subsystem.has_value()) {
        logger.exitWithError("expected diagonal initial condition to have list of subsystems");
      }
      n_initial_conditions = 1;
      for (size_t oscilID : init_cond_settings.subsystem.value()) {
        if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
      }
      // if Schroedinger solver: ninit = N, do nothing.
      // else Lindblad solver: ninit = N^2
      if (decoherence_type != DecoherenceType::NONE) {
        n_initial_conditions = (size_t)pow(n_initial_conditions, 2.0);
      }
      break;
  }
  logger.log("Number of initial conditions: " + std::to_string(n_initial_conditions) + "\n");
  return n_initial_conditions;
}

void Config::setRandSeed(int rand_seed_) {
  rand_seed = rand_seed_;
  if (rand_seed < 0) {
    std::random_device rd;
    rand_seed = rd(); // random non-reproducable seed
  }
}

ControlParameterizationSettings Config::parseControlParameterizationSpecs(const toml::table& param_table) const {
  std::string type_str = validators::field<std::string>(param_table, "type").value();
  auto type_enum = parseEnum(type_str, CONTROL_TYPE_MAP);
  if (!type_enum.has_value()) {
    logger.exitWithError("Unknown control parameterization type: " + type_str);
  }

  ControlParameterizationSettings param;
  param.type = type_enum.value();

  switch (param.type) {
    case ControlType::BSPLINE:
    case ControlType::BSPLINE0:
      param.nspline = validators::field<size_t>(param_table, "num").value();
      param.tstart = validators::getOptional<double>(param_table["tstart"]);
      param.tstop = validators::getOptional<double>(param_table["tstop"]);
      break;

    case ControlType::BSPLINEAMP:
      param.nspline = validators::field<size_t>(param_table, "num").value();
      param.scaling = validators::field<double>(param_table, "scaling").value();
      param.tstart = validators::getOptional<double>(param_table["tstart"]);
      param.tstop = validators::getOptional<double>(param_table["tstop"]);
      break;

    case ControlType::NONE:
      break;
  }

  return param;
}


ControlInitializationSettings Config::parseControlInitializationSpecs(const toml::table& init_table) const {
  std::string type = validators::field<std::string>(init_table, "type").value();
  auto type_enum = parseEnum(type, CONTROL_INITIALIZATION_TYPE_MAP);
  if (!type_enum.has_value()) {
    logger.exitWithError("Unknown control initialization type: " + type);
  }

  ControlInitializationSettings init;
  init.type = type_enum.value();

  if (init.type == ControlInitializationType::FILE) {
    init.filename = validators::field<std::string>(init_table, "filename").value();
    if (!init.filename.has_value()) {
      logger.exitWithError("control_initialization of type 'file' must have a 'filename' parameter");
    }
  } else {
    init.amplitude = validators::field<double>(init_table, "amplitude").valueOr(ConfigDefaults::CONTROL_INIT_AMPLITUDE);
    init.phase = validators::field<double>(init_table, "phase").greaterThanEqual(0.0).valueOr(ConfigDefaults::CONTROL_INIT_PHASE);
  }

  return init;
}


OptimTargetSettings Config::parseOptimTarget(const toml::table& toml, size_t num_osc) const {
  OptimTargetSettings optim_target;

  if (toml.contains("target")) {
    if (!toml["target"].as_table()) {
      logger.exitWithError("target must be a table");
    }
    const auto* target_table = toml["target"].as_table();

    // Get the target type, or default to NONE
    auto type_str = validators::field<std::string>(*target_table, "type").valueOr("none");
    auto type_opt = parseEnum(type_str, TARGET_TYPE_MAP);
    if (!type_opt.has_value()) {
      logger.exitWithError("Unknown optim_target type: " + type_str);
    }
    optim_target.type = type_opt.value();

    // Parse other settings based on type
    if (optim_target.type == TargetType::GATE) {
      // For Gate target: Either gate_type or filename needs to be provided
      auto gate_type_str = validators::field<std::string>(*target_table, "gate_type").valueOr("none");
      optim_target.gate_type = parseEnum(gate_type_str, GATE_TYPE_MAP);
      optim_target.filename = validators::getOptional<std::string>((*target_table)["filename"]);
      // Make sure either gate_type or filename is provided
      if (!optim_target.gate_type.has_value() && !optim_target.filename.has_value()) {
        logger.exitWithError("For optim_target of type 'gate', either gate_type or filename must be specified");
      }
      // Prioritize gate from file.
      if (optim_target.filename.has_value()) {
        optim_target.gate_type = GateType::FILE;
      }

      // For gate, check for optional gate rotation frequencies
      optim_target.gate_rot_freq = validators::scalarOrVectorOr<double>(*target_table, "gate_rot_freq", num_osc, std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));

    } else if (optim_target.type == TargetType::STATE) {
      // State target: Either levels for product state or filename needs to be provided
      optim_target.filename = validators::getOptional<std::string>((*target_table)["filename"]);
      optim_target.levels = validators::getOptionalVector<size_t>((*target_table)["levels"]);
      // Validate levels, if provided
      if (optim_target.levels.has_value()) {
        if (optim_target.levels->size() != nlevels.size()) {
          logger.exitWithError("optim_target levels size does not match number of oscillators");
        }
        for (size_t i = 0; i < nlevels.size(); i++) {
          if (optim_target.levels->at(i) >= nlevels[i]) {
            logger.exitWithError("ERROR in config setting. The requested product state target |" + std::to_string(optim_target.levels->at(i)) +"> exceeds the number of modeled levels for that oscillator (" + std::to_string(nlevels[i]) + ").\n");
          }
        }
      }
      // Make sure either levels or filename is provided
      if (!optim_target.levels.has_value() && !optim_target.filename.has_value()) {
        logger.exitWithError("For optim_target of type 'state', either levels or filename must be specified");
      }
      // Prioritize state from file.
      if (optim_target.filename.has_value()) {
        optim_target.levels = std::nullopt;
      }
    }
  }

  return optim_target;
}


// CFG parsing helpers
// TODO cfg: delete these when .cfg format is removed.

template <typename T>
std::vector<std::vector<T>> Config::parseOscillatorSettingsCfg(
    const std::optional<std::map<int, std::vector<T>>>& indexed, size_t num_entries,
    const std::vector<T>& default_values) const {
  // Start with all defaults
  std::vector<std::vector<T>> result(num_entries, default_values);

  // Overwrite with specified values
  if (indexed.has_value()) {
    for (const auto& [idx, vals] : *indexed) {
      if (idx >= 0 && static_cast<size_t>(idx) < num_entries) {
        result[idx] = vals;
      }
    }
  }
  return result;
}

std::vector<ControlParameterizationSettings> Config::parseControlParameterizationsCfg(const std::optional<std::map<int, ControlParameterizationData>>& parameterizations_map) const {
  // Use default-initialized struct (defaults provided in struct definition)
  ControlParameterizationSettings default_parameterization;

  // Populate default if paramterization is not specified
  if (!parameterizations_map.has_value()) {
    return std::vector<ControlParameterizationSettings>(nlevels.size(), default_parameterization);
  }

  // Otherwise, parse specified parameterizations for each oscillator
  auto parsed_parameterizations = std::vector<ControlParameterizationSettings>(nlevels.size(), default_parameterization);
  for (size_t i = 0; i < parsed_parameterizations.size(); i++) {
    if (parameterizations_map.value().find(static_cast<int>(i)) != parameterizations_map.value().end()) {

      // auto parameterization = parseControlParameterizationCfg(parameterizations_map.value().at(i));
      auto oscil_config = parameterizations_map.value().at(static_cast<int>(i));
      const auto& params = oscil_config.parameters;

      // Create and store the parameterization
      ControlParameterizationSettings parameterization;
      parameterization.type = oscil_config.control_type;
      if (oscil_config.control_type == ControlType::BSPLINE || oscil_config.control_type == ControlType::BSPLINE0) {
        assert(params.size() >= 1); // nspline is required, should be validated in CfgParser
        parameterization.nspline = static_cast<size_t>(params[0]);
        parameterization.tstart = params.size() > 1 ? std::optional<double>(params[1]) : std::nullopt;
        parameterization.tstop = params.size() > 2 ? std::optional<double>(params[2]) : std::nullopt;
      } else if (oscil_config.control_type == ControlType::BSPLINEAMP) {
        assert(params.size() >= 2); // nspline and scaling are required, should be validated in CfgParser
        parameterization.nspline = static_cast<size_t>(params[0]);
        parameterization.scaling = static_cast<double>(params[1]);
        parameterization.tstart = params.size() > 2 ? std::optional<double>(params[2]) : std::nullopt;
        parameterization.tstop = params.size() > 3 ? std::optional<double>(params[3]) : std::nullopt;
      }
      parsed_parameterizations[i] = parameterization;
    }
  }
  return parsed_parameterizations;
}


std::vector<ControlInitializationSettings> Config::parseControlInitializationsCfg(const std::optional<std::map<int, ControlInitializationSettings>>& init_configs) const {

  ControlInitializationSettings default_init;

  std::vector<ControlInitializationSettings> control_initializations(nlevels.size(), default_init);

  if (init_configs.has_value()) {
    for (size_t i = 0; i < nlevels.size(); i++) {
      if (init_configs->find(static_cast<int>(i)) != init_configs->end()) {
        control_initializations[i] = init_configs->at(static_cast<int>(i));
      }
    }
  }

  return control_initializations;
}
