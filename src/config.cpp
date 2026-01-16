#include "config.hpp"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

#include "config_defaults.hpp"
#include "config_validators.hpp"
#include "util.hpp"

// Common TOML table key constants
namespace {
const std::string OSC_ID_KEY = "oscID";

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
      // Check if subsystem field is an array or a simple index
      size_t index;
      if (elem_table->get(subsystem_key)->is_array()) {
        // Coupling parameter case: subsystem is an array of two indices
        auto subsys_array = validators::vectorField<size_t>(*elem_table, subsystem_key).hasLength(2).value();
        size_t i = subsys_array[0];;
        size_t j = subsys_array[1];
        // if (i >= num_entries || j >= num_entries) {
        //   throw validators::ValidationError(key, "subsystem index out of range for key '" + key + "'");
        // }
        // Compute unique index for pair (i,j) with i<j: Convert to linear index: (0,1), (0,2), ..., (0,n-1), (1,2), ..., (1,n-1), ..., (n-2,n-1)
        if (i > j) std::swap(i, j);
        index = i * num_subsystems - i - i * (i - 1) / 2 + (j - i - 1);
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
    // General options

    nlevels = validators::vectorField<size_t>(toml, "nlevels").minLength(1).positive().value();
    size_t num_osc = nlevels.size();

    nessential = validators::vectorField<size_t>(toml, "nessential").minLength(1).positive().valueOr(nlevels);
    copyLast(nessential, num_osc);

    ntime = validators::field<size_t>(toml, "ntime").positive().value();

    dt = validators::field<double>(toml, "dt").positive().value();

    transfreq = validators::vectorField<double>(toml, "transfreq").minLength(1).value();
    copyLast(transfreq, num_osc);

    selfkerr = validators::vectorField<double>(toml, "selfkerr").minLength(1).valueOr(std::vector<double>(num_osc, ConfigDefaults::SELFKERR));
    copyLast(selfkerr, num_osc);

    // Parse crosskerr and Jkl coupling: either one value (all-to-all coupling) or array of tables with 'subsystem = [i,j]' field for i-j coupling)
    size_t num_pairs = (num_osc - 1) * num_osc / 2;
    crosskerr.assign(num_pairs, ConfigDefaults::CROSSKERR);
    Jkl.assign(num_pairs, ConfigDefaults::JKL);
    // Overwrite for crosskerr
    if (toml.contains("crosskerr")) {
      // Check if crosskerr is a value (all-to-all coupling)
      if (toml["crosskerr"].is_value()) {
        double single_val = validators::field<double>(toml, "crosskerr").value();
        crosskerr.assign(num_pairs, single_val);
      } else {
      // Parse as array of tables format
      auto parseCouplingFunc = [this](const toml::table& t) { return parseCouplingParameterSpecs(t, "crosskerr"); };
      crosskerr = parsePerSubsystemSettings<double>(toml, "crosskerr", num_osc, ConfigDefaults::CROSSKERR, parseCouplingFunc, logger);
      }
    }
    // Overwrite for Jkl
    if (toml.contains("Jkl")) {
      // Check if Jkl is a value (all-to-all coupling)
      if (toml["Jkl"].is_value()) {
        double single_val = validators::field<double>(toml, "Jkl").value();
        Jkl.assign(num_pairs, single_val);
      } else {
      // Parse as array of tables format
      auto parseCouplingFunc = [this](const toml::table& t) { return parseCouplingParameterSpecs(t, "Jkl"); };
      Jkl = parsePerSubsystemSettings<double>(toml, "Jkl", num_osc, ConfigDefaults::JKL, parseCouplingFunc, logger);
      }
    }

    rotfreq = validators::vectorField<double>(toml, "rotfreq").minLength(1).valueOr(std::vector<double>(num_osc, ConfigDefaults::ROTFREQ));
    copyLast(rotfreq, num_osc);

    hamiltonian_file_Hsys = validators::getOptional<std::string>(toml["hamiltonian_file_Hsys"]);
    hamiltonian_file_Hc = validators::getOptional<std::string>(toml["hamiltonian_file_Hc"]);

    if (!toml.contains("decoherence")) {
      // No decoherence table specified, use defaults
      decoherence_type = ConfigDefaults::DECOHERENCE_TYPE;
      decay_time = std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME);
      dephase_time = std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME);
    } else {
      // Parse decoherence table
      auto* decoherence_table = toml["decoherence"].as_table();
      if (!decoherence_table) {
        logger.exitWithError("decoherence must be a table");
      }
      auto type_str = validators::field<std::string>(*decoherence_table, "type").valueOr("none");
      decoherence_type = parseEnum(type_str, DECOHERENCE_TYPE_MAP, ConfigDefaults::DECOHERENCE_TYPE);
      decay_time = validators::vectorField<double>(*decoherence_table, "decay_time").valueOr(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
      dephase_time = validators::vectorField<double>(*decoherence_table, "dephase_time").valueOr(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
    }

    // Parse initial condition table
    auto init_cond_table = validators::getRequiredTable(toml, "initial_condition");
    auto type_opt = parseEnum(validators::field<std::string>(init_cond_table, "type").value(), INITCOND_TYPE_MAP);
    if (!type_opt.has_value()) {
      logger.exitWithError("initial condition type not found.");
    }
    initial_condition.type = type_opt.value();
    initial_condition.levels = validators::getOptionalVector<size_t>(init_cond_table["levels"]);
    initial_condition.filename = validators::getOptional<std::string>(init_cond_table["filename"]);
    initial_condition.osc_IDs = validators::getOptionalVector<size_t>(init_cond_table["oscIDs"]);

    // Optimization options
    control_enforceBC = toml["control_enforceBC"].value_or(ConfigDefaults::CONTROL_ENFORCE_BC);

    // Parse control parameterization, either table (applies to all) or array (per-oscillator)
    ControlParameterizationSettings default_param;
    default_param.type = ConfigDefaults::CONTROL_TYPE;
    default_param.nspline = ConfigDefaults::CONTROL_SPLINE_COUNT;
    default_param.tstart = std::nullopt;
    default_param.tstop = std::nullopt;
    control_parameterizations.assign(num_osc, default_param);
    if (toml.contains("control_parameterization")) {
      auto parseParamFunc = [this](const toml::table& t) { return parseControlParameterizationSpecs(t); };
      control_parameterizations = parsePerSubsystemSettings<ControlParameterizationSettings>(toml, "control_parameterization", num_osc, default_param, parseParamFunc, logger);
    }

    // Parse control initialization, either as a table (applies to all) or per-oscillator table
    ControlInitializationSettings default_init;
    default_init.type = ConfigDefaults::CONTROL_INIT_TYPE; 
    default_init.amplitude = ConfigDefaults::CONTROL_INIT_AMPLITUDE; 
    default_init.phase = std::nullopt; 
    default_init.filename = std::nullopt;
    control_initializations.assign(num_osc, default_init);
    if (toml.contains("control_initialization")) {
      auto parseInitFunc = [this](const toml::table& t) { return parseControlInitializationSpecs(t); };
      control_initializations = parsePerSubsystemSettings<ControlInitializationSettings>(toml, "control_initialization", num_osc, default_init, parseInitFunc, logger);
    }

    // Parse optional control bounds: either single value or per-oscillator array
    control_bounds.assign(num_osc, ConfigDefaults::CONTROL_BOUND);
    if (toml.contains("control_bounds")) {
      if (toml["control_bounds"].as_array()) {
        // Get control_bounds from array
        control_bounds = validators::vectorField<double>(toml, "control_bounds").minLength(1).value();
        copyLast(control_bounds, num_osc);
        control_bounds.resize(num_osc);
      } else if (toml["control_bounds"].is_value()){
        // Get single control_bounds value
        auto single_val = validators::field<double>(toml, "control_bounds").value();
        control_bounds = std::vector<double>(num_osc, single_val);
      } else {
        logger.exitWithError("control_bounds must be either a single value (applies to all oscillators), or an array of values (one per oscillator).");
      }
    }

    // Parse carrier frequencies: either one vector (applies to all oscillators) or per-oscillator array of tables
    std::vector<double> default_carrier_freq = {ConfigDefaults::CARRIER_FREQ};
    carrier_frequencies.assign(num_osc, default_carrier_freq);
    if (toml.contains("carrier_frequency")) {
      // Check if carrier_frequency is a direct array of values (shorthand for applying to all)
      auto* carrier_freq_array = toml["carrier_frequency"].as_array();
      if (carrier_freq_array && !carrier_freq_array->empty() && !carrier_freq_array->front().is_table()) {
        // Direct array format: carrier_frequency = [1.0, 2.0]
        auto values = validators::vectorField<double>(toml, "carrier_frequency").value();
        carrier_frequencies.assign(num_osc, values);
      } else {
        // Table or array of tables format
        auto parseCarrierFunc = [this](const toml::table& t) { return parseCarrierFrequencySpecs(t); };
        carrier_frequencies = parsePerSubsystemSettings<std::vector<double>>(toml, "carrier_frequency", num_osc, default_carrier_freq, parseCarrierFunc, logger);
      }
    } 

    // Parse optimization target:
    optim_target = parseOptimTarget(toml, num_osc);

    optim_objective = parseEnum(toml["optim_objective"].value<std::string>(), OBJECTIVE_TYPE_MAP, ConfigDefaults::OPTIM_OBJECTIVE);

    // Parse optional optim_weights: can be a single value or a vector
    if (!toml.contains("optim_weights")) {
      optim_weights = std::vector<double>{ConfigDefaults::OPTIM_WEIGHT};
    } else {
      if (toml["optim_weights"].as_array()) {        
        // Get optim_weights from array
        optim_weights = validators::vectorField<double>(toml, "optim_weights").minLength(1).value();
      } else if (toml["optim_weights"].is_value()) { 
        // Get single value
        auto single_val = validators::field<double>(toml, "optim_weights").value();
        optim_weights = std::vector<double>{single_val};
      } else {
        logger.exitWithError("optim_weights must be either a single value (applies to all initial conditions), or an array of values (one for each initial condition).");
      }
    }

    // Parse optional optimization tolerances
    if (!toml.contains("optim_tolerance")) {
      optim_tol_grad_abs = ConfigDefaults::OPTIM_TOL_GRAD_ABS;
      optim_tol_grad_rel = ConfigDefaults::OPTIM_TOL_GRAD_REL;
      optim_tol_finalcost = ConfigDefaults::OPTIM_TOL_FINALCOST;
      optim_tol_infidelity = ConfigDefaults::OPTIM_TOL_INFIDELITY;
    } else {
      // Parse optim_tolerance table
      auto* tol_table = toml["optim_tolerance"].as_table();
      if (!tol_table) {
        logger.exitWithError("optim_tolerance must be a table");
      }
      optim_tol_grad_abs = validators::field<double>(*tol_table, "grad_abs").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_ABS);
      optim_tol_grad_rel = validators::field<double>(*tol_table, "grad_rel").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_REL);
      optim_tol_finalcost = validators::field<double>(*tol_table, "final_cost").positive().valueOr(ConfigDefaults::OPTIM_TOL_FINALCOST);
      optim_tol_infidelity = validators::field<double>(*tol_table, "infidelity").positive().valueOr(ConfigDefaults::OPTIM_TOL_INFIDELITY);
    }

    optim_maxiter = validators::field<size_t>(toml, "optim_maxiter").greaterThanEqual(0).valueOr(ConfigDefaults::OPTIM_MAXITER);

    // Parse optim_tikhonov inline table
    if (!toml.contains("optim_tikhonov")) {
      optim_tikhonov_coeff = ConfigDefaults::OPTIM_TIKHONOV_COEFF;
      optim_tikhonov_use_x0 = ConfigDefaults::OPTIM_TIKHONOV_USE_X0;
    } else {
      auto regul_table = toml["optim_tikhonov"].as_table();
      if (!regul_table) {
        logger.exitWithError("optim_tikhonov must be a table");
      }
      optim_tikhonov_coeff= validators::field<double>(*regul_table, "coeff").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_TIKHONOV_COEFF);
      optim_tikhonov_use_x0 = validators::field<bool>(*regul_table, "use_x0").valueOr(ConfigDefaults::OPTIM_TIKHONOV_USE_X0);
    }

    // Parse optim_penalty table
    if (!toml.contains("optim_penalty")) {
      // No optim_penalty table is specified, use defaults
      optim_penalty_leakage = ConfigDefaults::OPTIM_PENALTY_LEAKAGE;
      optim_penalty_weightedcost = ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST;
      optim_penalty_weightedcost_width = ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH;
      optim_penalty_dpdm = ConfigDefaults::OPTIM_PENALTY_DPDM;
      optim_penalty_energy = ConfigDefaults::OPTIM_PENALTY_ENERGY;
      optim_penalty_variation = ConfigDefaults::OPTIM_PENALTY_VARIATION;
    } else {
      // Parse optim_penalty table
      auto penalty_table = toml["optim_penalty"].as_table();
      if (!penalty_table) {
        logger.exitWithError("optim_penalty must be a table");
      }
      optim_penalty_leakage = validators::field<double>(*penalty_table, "leakage").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
      optim_penalty_weightedcost = validators::field<double>(*penalty_table, "weightedcost").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
      optim_penalty_weightedcost_width = validators::field<double>(*penalty_table, "weightedcost_width").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
      optim_penalty_dpdm = validators::field<double>(*penalty_table, "dpdm").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_DPDM);
      optim_penalty_energy = validators::field<double>(*penalty_table, "energy").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_ENERGY);
      optim_penalty_variation = validators::field<double>(*penalty_table, "variation").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_VARIATION);
    }

    datadir = toml["datadir"].value_or(ConfigDefaults::DATADIR);

    // Parse output_type as an array of strings (defaults to empty array)
    output_type.clear();
    if (auto output_type_array = toml["output_type"].as_array()) {
      for (auto&& elem : *output_type_array) {
        if (auto str = elem.value<std::string>()) {
          auto enum_val = parseEnum(*str, OUTPUT_TYPE_MAP);
          if (!enum_val.has_value()) {
            logger.exitWithError("Unknown output type: " + *str);
          }
          output_type.push_back(enum_val.value());
        } else {
          logger.exitWithError("output_type array must contain strings");
        }
      }
    }

    output_timestep_stride = validators::field<size_t>(toml, "output_timestep_stride").greaterThanEqual(0).valueOr(ConfigDefaults::OUTPUT_TIMESTEP_STRIDE);

    output_optimization_stride = validators::field<size_t>(toml, "output_optimization_stride").greaterThanEqual(0).valueOr(ConfigDefaults::OUTPUT_OPTIMIZATION_STRIDE);

    runtype = parseEnum(toml["runtype"].value<std::string>(), RUN_TYPE_MAP, ConfigDefaults::RUNTYPE);

    usematfree = toml["usematfree"].value_or(ConfigDefaults::USEMATFREE);

    // Parse linearsolver as an inline table
    if (!toml.contains("linearsolver")) {
      // No linearsolver table specified, use defaults
      linearsolver_type = ConfigDefaults::LINEARSOLVER_TYPE;
      linearsolver_maxiter = ConfigDefaults::LINEARSOLVER_MAXITER;
    } else {
      auto* linearsolver_table = toml["linearsolver"].as_table();
      if (!linearsolver_table) {
        logger.exitWithError("linearsolver must be a table");
      }
      linearsolver_type = parseEnum(validators::field<std::string>(*linearsolver_table, "type").value(), LINEAR_SOLVER_TYPE_MAP, ConfigDefaults::LINEARSOLVER_TYPE);
      linearsolver_maxiter = validators::field<size_t>(*linearsolver_table, "maxiter").positive().valueOr(ConfigDefaults::LINEARSOLVER_MAXITER);
    }

    timestepper_type = parseEnum(toml["timestepper"].value<std::string>(), TIME_STEPPER_TYPE_MAP, ConfigDefaults::TIMESTEPPER_TYPE);

    int rand_seed_ = toml["rand_seed"].value_or(ConfigDefaults::RAND_SEED);
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
  initial_condition.osc_IDs = settings.initialcondition.value().osc_IDs;

  // Control and optimization parameters
  control_enforceBC = settings.control_enforceBC.value_or(ConfigDefaults::CONTROL_ENFORCE_BC);

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

  // Parse control_bounds from CFG format (returns vector of vectors, but we need vector)
  auto control_bounds_cfg = parseOscillatorSettingsCfg<double>(settings.indexed_control_bounds, control_parameterizations.size(), {ConfigDefaults::CONTROL_BOUND});
  control_bounds.resize(control_bounds_cfg.size());
  for (size_t i = 0; i < control_bounds_cfg.size(); ++i) {
    control_bounds[i] = control_bounds_cfg[i].empty() ? ConfigDefaults::CONTROL_BOUND : control_bounds_cfg[i][0];
  }

  carrier_frequencies.resize(num_osc);
  carrier_frequencies = parseOscillatorSettingsCfg<double>(settings.indexed_carrier_frequencies, num_osc, {ConfigDefaults::CARRIER_FREQ});

  if (settings.optim_target.has_value()) {
    optim_target = settings.optim_target.value();
    optim_target.gate_rot_freq = settings.gate_rot_freq.value_or(std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));
  } else {
    // No optim_target specified, use default (no target)
    OptimTargetSettings default_target;
    default_target.type = ConfigDefaults::OPTIM_TARGET;
    default_target.gate_type = std::nullopt;
    default_target.gate_rot_freq = std::nullopt;
    default_target.levels = std::nullopt;
    default_target.filename = std::nullopt;
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
  datadir = settings.datadir.value_or(ConfigDefaults::DATADIR);
  
  // Convert old per-oscillator output to global output_type (apply to all oscillators)
  auto indexed_output_vec = parseOscillatorSettingsCfg<OutputType>(settings.indexed_output, num_osc);
  output_type.clear();
  // Collect unique output types from all oscillators
  std::set<OutputType> unique_types;
  for (const auto& osc_output : indexed_output_vec) {
    for (const auto& type : osc_output) {
      unique_types.insert(type);
    }
  }
  // Convert set to vector
  output_type.assign(unique_types.begin(), unique_types.end());
  
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

std::string toStringCoupling(const std::vector<double>& couplings, size_t num_osc) {
  if (couplings.empty()) return "[]";

  // If all couplings are the same, print single value
  bool all_same = true;
  for (size_t i = 1; i < couplings.size(); ++i) {
    if (couplings[i] != couplings[0]) {
      all_same = false;
      break;
    }
  }
  if (all_same) {
    std::ostringstream oss;
    oss << std::setprecision(std::numeric_limits<double>::max_digits10);
    oss << couplings[0];
    return oss.str();
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
  std::ostringstream oss;
  oss << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (size_t i = 0; i < nonzero_couplings.size(); ++i) {
    auto [pair, value] = nonzero_couplings[i];
    auto [first, second] = pair;
    oss.str("");  // Clear the stringstream for reuse
    oss << value;
    result += " { subsystem = [" + std::to_string(first) + "," + std::to_string(second) + "], value = " + oss.str() + "}";
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

  std::ostringstream oss;
  if constexpr (std::is_floating_point_v<T>) {
    oss << std::setprecision(std::numeric_limits<T>::max_digits10);
  }

  oss << "[" << vec[0];
  for (size_t i = 1; i < vec.size(); ++i) {
    oss << ", " << vec[i];
  }
  oss << "]";
  return oss.str();
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
      std::string out = "{" + type_str + ", oscIDs = ";
      out += printVector(initial_condition.osc_IDs.value());
      out += "}";
      return out;
    }
    case InitialConditionType::DIAGONAL: {
      std::string out = "{" + type_str + ", oscIDs = ";
      out += printVector(initial_condition.osc_IDs.value());
      out += "}";
      return out;
    }
    case InitialConditionType::BASIS: {
      std::string out = "{" + type_str + ", oscIDs = ";
      out += printVector(initial_condition.osc_IDs.value());
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
  bool all_same = true;
  const auto& first = items[0];
  for (size_t i = 1; i < items.size(); ++i) {
    if (!areEqual(items[i], first)) {
      all_same = false;
      break;
    }
  }
  
  if (all_same) {
    std::string out = printItems(first);
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
    out += init.amplitude.has_value() ? ", amplitude = " + std::to_string(init.amplitude.value()) : "";
    out += init.phase.has_value() ? ", phase = " + std::to_string(init.phase.value()) : "";
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
    out += param.tstart.has_value() ? ", tstart = " + std::to_string(param.tstart.value()) : "";
    out += param.tstop.has_value() ? ", tstop = " + std::to_string(param.tstop.value()) : "";
    out += param.scaling.has_value() ? ", scaling = " + std::to_string(param.scaling.value()) : "";
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
    return printVector(freqs);
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

  bool all_equal = true;
  for (size_t i = 1; i < vec.size(); ++i) {
    if (vec[i] != vec[0]) {
      all_equal = false;
      break;
    }
  }

  if (all_equal) {
    return std::to_string(vec[0]);
  } else {
    return printVector(vec);
  }
}

} // namespace

void Config::printConfig(std::stringstream& log) const {
  log << "# Configuration settings\n";
  log << "# =============================================\n\n";

  // Section 1: Root table settings

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
  // Control pulse settings
  log << "control_parameterization = " << toString(control_parameterizations) << "\n";
  log << "carrier_frequency = " << toString(carrier_frequencies) << "\n";
  log << "control_initialization = " << toString(control_initializations) << "\n";
  log << "control_bounds = " << toString(control_bounds) << "\n";
  log << "control_enforceBC = " << (control_enforceBC ? "true" : "false") << "\n";

  // Optimization parameters
  log << "optim_target = " << toString(optim_target) << "\n";
  log << "optim_objective = \"" << enumToString(optim_objective, OBJECTIVE_TYPE_MAP) << "\"\n";
  log << "optim_weights = " << toString(optim_weights) << "\n";
  log << "optim_tolerance = { grad_abs = " << optim_tol_grad_abs
      << ", grad_rel = " << optim_tol_grad_rel
      << ", final_cost = " << optim_tol_finalcost
      << ", infidelity = " << optim_tol_infidelity << " }\n";
  log << "optim_maxiter = " << optim_maxiter << "\n";
  log << "optim_tikhonov = { coeff = " << optim_tikhonov_coeff   
      << ", use_x0 = " << (optim_tikhonov_use_x0 ? "true" : "false") << " }\n";
  log << "optim_penalty = { leakage = " << optim_penalty_leakage
      << ", energy = " << optim_penalty_energy
      << ", dpdm = " << optim_penalty_dpdm
      << ", variation = " << optim_penalty_variation
      << ", weightedcost = " << optim_penalty_weightedcost
      << ", weightedcost_width = " << optim_penalty_weightedcost_width << " }\n";

  // Output  options 
  log << "datadir = \"" << datadir << "\"\n";
  log << "output_type = [";
  for (size_t j = 0; j < output_type.size(); ++j) {
    log << "\"" << enumToString(output_type[j], OUTPUT_TYPE_MAP) << "\"";
    if (j < output_type.size() - 1) log << ", ";
  }
  log << "]\n";
  log << "output_timestep_stride = " << output_timestep_stride << "\n";
  log << "output_optimization_stride = " << output_optimization_stride << "\n";

  // Solver options
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
    if (!initial_condition.osc_IDs.has_value()) {
      initial_condition.osc_IDs = std::vector<size_t>(nlevels.size());
      for (size_t i = 0; i < nlevels.size(); i++) {
        initial_condition.osc_IDs->at(i) = i;
      }
    }
  }

  // Compute number of initial conditions
  n_initial_conditions = computeNumInitialConditions(initial_condition, nlevels, nessential, decoherence_type);

  // overwrite decay or dephase times with zeros, if the decoherence type is only one of them.
  if (decoherence_type == DecoherenceType::DECAY) {
    std::fill(dephase_time.begin(), dephase_time.end(), 0);
  } else if (decoherence_type == DecoherenceType::DEPHASE) {
    std::fill(decay_time.begin(), decay_time.end(), 0);
  }

  // Scale optimization weights such that they sum up to one
  // If a single value was provided, replicate it for all initial conditions
  // If a vector was provided, verify it has the correct length
  if (optim_weights.size() == 1) {
    copyLast(optim_weights, n_initial_conditions);
  } else if (optim_weights.size() != n_initial_conditions) {
    logger.exitWithError("optim_weights vector has length " + std::to_string(optim_weights.size()) + " but must have length " + std::to_string(n_initial_conditions) + " (number of initial conditions)");
  }
  optim_weights.resize(n_initial_conditions);
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
  for (size_t i = 0; i < control_bounds.size(); i++) {
    if (control_bounds[i] <= 0.0) {
      logger.exitWithError("control_bounds[" + std::to_string(i) + "] must be positive");
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
    if (!initial_condition.osc_IDs.has_value()) {
      logger.exitWithError("initialcondition of type BASIS, DIAGONAL, or ENSEMBLE must have 'osc_IDs'");
    }
    if (initial_condition.osc_IDs->back() >= nlevels.size()) {
      logger.exitWithError("Last element in initialcondition params exceeds number of oscillators");
    }
    for (size_t i = 1; i < initial_condition.osc_IDs->size() - 1; i++) {
      if (initial_condition.osc_IDs->at(i) + 1 != initial_condition.osc_IDs->at(i + 1)) {
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
      if (!init_cond_settings.osc_IDs.has_value()) {
        logger.exitWithError("expected diagonal initial condition to have list of oscIDs");
      }
      n_initial_conditions = 1;
      for (size_t oscilID : init_cond_settings.osc_IDs.value()) {
        if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
      }
      break;
    case InitialConditionType::BASIS:
      /* Compute ninit = dim(subsystem defined by list of oscil IDs) */
      if (!init_cond_settings.osc_IDs.has_value()) {
        logger.exitWithError("expected diagonal initial condition to have list of oscIDs");
      }
      n_initial_conditions = 1;
      for (size_t oscilID : init_cond_settings.osc_IDs.value()) {
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

void Config::validateTableKeys(const toml::table& table, const std::set<std::string>& allowed_keys, const std::string& table_name) const {
  for (const auto& [key, _] : table) {
    if (allowed_keys.find(std::string(key.str())) == allowed_keys.end()) {
      logger.exitWithError("Unknown key '" + std::string(key.str()) + "' in " + table_name + ".");
    }
  }
}

double Config::parseCouplingParameterSpecs(const toml::table& table, const std::string& key) const {
  const std::string subsystem_key = "subsystem";
  const std::string value_key = "value";
  const std::set<std::string> allowed_keys = {subsystem_key, value_key};
  
  validateTableKeys(table, allowed_keys, key);
  
  return validators::field<double>(table, value_key).value();
}

std::vector<double> Config::parseCarrierFrequencySpecs(const toml::table& table) const {
  const std::string subsystem_key = "subsystem";
  const std::string values_key = "value";
  const std::set<std::string> allowed_keys = {subsystem_key, values_key};
  
  validateTableKeys(table, allowed_keys, "carrier_frequency");
  
  return validators::vectorField<double>(table, values_key).value();
}

ControlParameterizationSettings Config::parseControlParameterizationSpecs(const toml::table& param_table) const {
  
  const std::string subsystem_key = "subsystem";
  const std::string type_key = "type";
  const std::string num_key = "num";
  const std::string scaling_key = "scaling";
  const std::string tstart_key = "tstart";
  const std::string tstop_key = "tstop";
  const std::set<std::string> allowed_keys = {subsystem_key, type_key, num_key, scaling_key, tstart_key, tstop_key};

  validateTableKeys(param_table, allowed_keys, "control_parameterization");
  
  // Parse type of parameterization
  std::string type_str = validators::field<std::string>(param_table, type_key).value();
  auto type_enum = parseEnum(type_str, CONTROL_TYPE_MAP);
  if (!type_enum.has_value()) {
    logger.exitWithError("Unknown control parameterization type: " + type_str);
  }

  ControlParameterizationSettings param;
  param.type = type_enum.value();
  
  // Parse other parameters based on type
  switch (param.type) {
    case ControlType::BSPLINE:
    case ControlType::BSPLINE0:
      param.nspline = validators::field<size_t>(param_table, num_key).value();
      param.tstart = validators::getOptional<double>(param_table[tstart_key]);
      param.tstop = validators::getOptional<double>(param_table[tstop_key]);
      break;
      
    case ControlType::BSPLINEAMP:
      param.nspline = validators::field<size_t>(param_table, num_key).value();
      param.scaling = validators::field<double>(param_table, scaling_key).value();
      param.tstart = validators::getOptional<double>(param_table[tstart_key]);
      param.tstop = validators::getOptional<double>(param_table[tstop_key]);
      break;
      
    case ControlType::NONE:
      break;
  }
  
  return param;
}


ControlInitializationSettings Config::parseControlInitializationSpecs(const toml::table& init_table) const {

  const std::string subsystem_key = "subsystem";
  const std::string type_key = "type";
  const std::string filename_key = "filename";
  const std::string amplitude_key = "amplitude";
  const std::string phase_key = "phase";
  const std::set<std::string> allowed_keys = {subsystem_key, type_key, filename_key, amplitude_key, phase_key};

  validateTableKeys(init_table, allowed_keys, "control_initialization");
  
  // Parse type of initialization
  std::string type = validators::field<std::string>(init_table, type_key).value();
  auto type_enum = parseEnum(type, CONTROL_INITIALIZATION_TYPE_MAP);
  if (!type_enum.has_value()) {
    logger.exitWithError("Unknown control initialization type: " + type);
  }

  ControlInitializationSettings init;
  init.type = type_enum.value();

  // Parse other parameters based on type
  if (init.type == ControlInitializationType::FILE) {
    init.filename = validators::field<std::string>(init_table, filename_key).value();
    if (!init.filename.has_value()){
      logger.exitWithError("control_initialization of type 'file' must have a 'filename' parameter");
    }
  } else {
    init.amplitude = validators::field<double>(init_table, amplitude_key).valueOr(ConfigDefaults::CONTROL_INIT_AMPLITUDE);
    init.phase = validators::field<double>(init_table, phase_key).greaterThanEqual(0.0).valueOr(ConfigDefaults::CONTROL_INIT_PHASE);
  }

  return init;
}


OptimTargetSettings Config::parseOptimTarget(const toml::table& toml, size_t num_osc) const {

  const std::string type_key = "type";
  const std::string gate_type_key = "gate_type";
  const std::string filename_key = "filename";
  const std::string gate_rot_freq_key = "gate_rot_freq";
  const std::string levels_key = "levels";
  const std::set<std::string> allowed_keys = {type_key, gate_type_key, filename_key, gate_rot_freq_key, levels_key};

  // Initialize the result with defaults (no target)
  OptimTargetSettings optim_target{ConfigDefaults::OPTIM_TARGET, std::nullopt, std::nullopt, std::nullopt, std::nullopt};

  // If optim_target is specified, parse the provided table
  if (toml.contains("optim_target")) {
    if (!toml["optim_target"].as_table()) {
      logger.exitWithError("optim_target must be a table");
    }
    auto* target_table = toml["optim_target"].as_table();
    validateTableKeys(*target_table, allowed_keys, "optim_target");

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
      auto gate_rot_freq_opt = validators::getOptionalVector<double>((*target_table)["gate_rot_freq"]);
      optim_target.gate_rot_freq = gate_rot_freq_opt.value_or(std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));
      copyLast(optim_target.gate_rot_freq.value(), num_osc);

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
  // Set up a default parameterization for one oscillator
  ControlParameterizationSettings default_parameterization;
  default_parameterization.type = ConfigDefaults::CONTROL_TYPE;
  default_parameterization.nspline = ConfigDefaults::CONTROL_SPLINE_COUNT;

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

  ControlInitializationSettings default_init = ControlInitializationSettings{ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt};

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
