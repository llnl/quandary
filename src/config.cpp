#include "config.hpp"

namespace {

// TOML extraction helpers

/**
 * @brief Helper to get readable type names for error messages
 */
template <typename T>
std::string getTypeName() {
  if constexpr (std::is_same_v<T, int> || std::is_same_v<T, int64_t> || std::is_same_v<T, size_t>) {
    return "integer";
  } else if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
    return "number";
  } else if constexpr (std::is_same_v<T, std::string>) {
    return "string";
  } else if constexpr (std::is_same_v<T, bool>) {
    return "boolean";
  } else {
    return "unknown type";
  }
}

/**
 * @brief Extracts a scalar from TOML, returning nullopt if not present.
 * @throws ValidationError if the field exists but has the wrong type.
 */
template <typename T>
std::optional<T> extractToml(const toml::table& config, const std::string& key) {
  if (!config.contains(key)) {
    return std::nullopt;
  }
  auto val = config[key].template value<T>();
  if (!val) {
    throw validators::ValidationError(key, "wrong type (expected " + getTypeName<T>() + ")");
  }
  return val;
}

/**
 * @brief Extracts a vector from TOML, returning nullopt if not present.
 * @throws ValidationError if the field exists but is not an array or has wrong element types.
 */
template <typename T>
std::optional<std::vector<T>> extractTomlVector(const toml::table& config, const std::string& key) {
  if (!config.contains(key)) {
    return std::nullopt;
  }
  auto* arr = config[key].as_array();
  if (!arr) {
    throw validators::ValidationError(key, "wrong type (expected array)");
  }
  std::vector<T> result;
  for (size_t i = 0; i < arr->size(); ++i) {
    auto val = arr->at(i).template value<T>();
    if (!val) {
      std::ostringstream oss;
      oss << "element [" << i << "] wrong type (expected " << getTypeName<T>() << ")";
      throw validators::ValidationError(key, oss.str());
    }
    result.push_back(*val);
  }
  return result;
}

/**
 * @brief Extracts a field that can be scalar (broadcasts to vector) or array.
 * Returns nullopt if not present. Expands scalar to vector of expected_size.
 * @throws ValidationError if the field has wrong type or array has wrong size.
 */
template <typename T>
std::optional<std::vector<T>> extractScalarOrVector(const toml::table& config, const std::string& key, size_t expected_size) {
  if (!config.contains(key)) {
    return std::nullopt;
  }
  if (auto* arr = config[key].as_array()) {
    if (arr->size() != expected_size) {
      std::ostringstream oss;
      oss << "array must have exactly " << expected_size << " elements, got " << arr->size();
      throw validators::ValidationError(key, oss.str());
    }
    std::vector<T> result;
    for (size_t i = 0; i < arr->size(); ++i) {
      auto val = arr->at(i).template value<T>();
      if (!val) {
        std::ostringstream oss;
        oss << "element [" << i << "] wrong type (expected " << getTypeName<T>() << ")";
        throw validators::ValidationError(key, oss.str());
      }
      result.push_back(*val);
    }
    return result;
  }
  if (auto val = config[key].template value<T>()) {
    return std::vector<T>(expected_size, *val);
  }
  throw validators::ValidationError(key, "must be either a scalar value or an array of " + getTypeName<T>());
}

/**
 * @brief Extracts an optional scalar value from a TOML node.
 */
template <typename T, typename NodeType>
std::optional<T> getOptional(const toml::node_view<NodeType>& node) {
  return node.template value<T>();
}

/**
 * @brief Extracts an optional vector from a TOML node.
 */
template <typename T, typename NodeType>
std::optional<std::vector<T>> getOptionalVector(const toml::node_view<NodeType>& node) {
  auto* arr = node.as_array();
  if (!arr) {
    return std::nullopt;
  }
  std::vector<T> result;
  for (size_t i = 0; i < arr->size(); ++i) {
    auto val = arr->at(i).template value<T>();
    if (!val) {
      return std::nullopt;
    }
    result.push_back(*val);
  }
  return result;
}

// Parsing helpers

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
        auto subsys_array_opt = extractTomlVector<size_t>(*elem_table, subsystem_key);
        auto subsys_array = validators::vectorField<size_t>(subsys_array_opt, subsystem_key).hasLength(2).value();
        size_t i = subsys_array[0];
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
        auto index_opt = extractToml<size_t>(*elem_table, subsystem_key);
        index = validators::field<size_t>(index_opt, subsystem_key).lessThan(num_subsystems).value();
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

// TOML parsing helpers for structured settings

ControlParameterizationSettings parseControlParamSpecsToml(const toml::table& param_table) {
  auto type_str_opt = extractToml<std::string>(param_table, "type");
  auto type_str = validators::field<std::string>(type_str_opt, "parameterization.type").value();

  auto type_enum = parseEnum(type_str, CONTROL_TYPE_MAP);
  if (!type_enum.has_value()) {
    throw validators::ValidationError("parameterization.type", "unknown control parameterization type: " + type_str);
  }

  ControlParameterizationSettings param;
  param.type = type_enum.value();
  param.nspline = extractToml<size_t>(param_table, "num");
  param.scaling = extractToml<double>(param_table, "scaling");
  param.tstart = extractToml<double>(param_table, "tstart");
  param.tstop = extractToml<double>(param_table, "tstop");
  return param;
}

ControlInitializationSettings parseControlInitSpecsToml(const toml::table& init_table) {
  auto type_str_opt = extractToml<std::string>(init_table, "type");
  auto type_str = validators::field<std::string>(type_str_opt, "initialization.type").value();

  auto type_enum = parseEnum(type_str, CONTROL_INITIALIZATION_TYPE_MAP);
  if (!type_enum.has_value()) {
    throw validators::ValidationError("initialization.type", "unknown control initialization type: " + type_str);
  }

  ControlInitializationSettings init;
  init.type = type_enum.value();
  init.filename = extractToml<std::string>(init_table, "filename");
  init.amplitude = extractToml<double>(init_table, "amplitude");
  init.phase = extractToml<double>(init_table, "phase");
  return init;
}

OptimTargetSettings parseOptimTargetToml(const toml::table& toml, size_t num_osc) {
  OptimTargetSettings target;

  if (toml.contains("target")) {
    if (toml["target"].as_table() == nullptr) {
      throw validators::ValidationError("target", "must be a table");
    }
    const auto* target_table = toml["target"].as_table();

    auto type_str_opt = extractToml<std::string>(*target_table, "type");
    auto type_str = validators::field<std::string>(type_str_opt, "target.type").valueOr("none");

    auto type_opt = parseEnum(type_str, TARGET_TYPE_MAP);
    if (!type_opt.has_value()) {
      throw validators::ValidationError("target.type", "unknown optim_target type: " + type_str);
    }
    target.type = type_opt.value();

    auto gate_type_str_opt = extractToml<std::string>(*target_table, "gate_type");
    if (gate_type_str_opt.has_value()) {
      target.gate_type = parseEnum(*gate_type_str_opt, GATE_TYPE_MAP);
      if (!target.gate_type.has_value()) {
        throw validators::ValidationError("target.gate_type", "unknown gate type: " + *gate_type_str_opt);
      }
    }

    target.filename = extractToml<std::string>(*target_table, "filename");
    target.levels = getOptionalVector<size_t>((*target_table)["levels"]);
    target.gate_rot_freq = extractScalarOrVector<double>(*target_table, "gate_rot_freq", num_osc);
  }

  return target;
}

/**
 * @brief Extracts configuration from TOML into a ConfigInput struct.
 *
 * This function handles all TOML-specific extraction:
 * - Scalar-to-vector expansion based on num_osc
 * - String-to-enum conversion
 * - Per-subsystem settings parsing
 * - Nested table handling
 *
 * Validation of values happens in the ConfigInput constructor.
 */
ConfigInput extractConfigInput(const toml::table& toml, const MPILogger& logger) {
  try {
    ConfigInput input;

    // Get section tables - only [system] is required
    const auto* system_table = toml["system"].as_table();
    if (system_table == nullptr) {
      throw validators::ValidationError("system", "table is required");
    }

    // Other tables are optional - use empty table as fallback
    const auto control_table = toml["control"].is_table() ? *toml["control"].as_table() : toml::table{};
    const auto optimization_table = toml["optimization"].is_table() ? *toml["optimization"].as_table() : toml::table{};
    const auto output_table = toml["output"].is_table() ? *toml["output"].as_table() : toml::table{};
    const auto solver_table = toml["solver"].is_table() ? *toml["solver"].as_table() : toml::table{};

    // Extract nlevels first - needed for scalar-to-vector expansion
    input.nlevels = extractTomlVector<size_t>(*system_table, "nlevels");
    if (!input.nlevels.has_value() || input.nlevels->empty()) {
      throw validators::ValidationError("nlevels", "field not found or empty");
    }
    size_t num_osc = input.nlevels->size();
    size_t num_pairs = (num_osc - 1) * num_osc / 2;

    // System parameters
    input.nessential = extractScalarOrVector<size_t>(*system_table, "nessential", num_osc);
    input.ntime = extractToml<size_t>(*system_table, "ntime");
    input.dt = extractToml<double>(*system_table, "dt");
    input.transfreq = extractScalarOrVector<double>(*system_table, "transfreq", num_osc);
    input.selfkerr = extractScalarOrVector<double>(*system_table, "selfkerr", num_osc);

    // Parse crosskerr coupling
    if (system_table->contains("crosskerr")) {
      if ((*system_table)["crosskerr"].is_value()) {
        auto single_val = extractToml<double>(*system_table, "crosskerr");
        if (single_val.has_value()) {
          input.crosskerr = std::vector<double>(num_pairs, *single_val);
        }
      } else {
        auto parseFunc = [](const toml::table& t) {
          auto val_opt = extractToml<double>(t, "value");
          return validators::field<double>(val_opt, "value").value();
        };
        input.crosskerr = parsePerSubsystemSettings<double>(*system_table, "crosskerr", num_osc, ConfigDefaults::CROSSKERR, parseFunc, logger);
      }
    }

    // Parse Jkl coupling
    if (system_table->contains("Jkl")) {
      if ((*system_table)["Jkl"].is_value()) {
        auto single_val = extractToml<double>(*system_table, "Jkl");
        if (single_val.has_value()) {
          input.Jkl = std::vector<double>(num_pairs, *single_val);
        }
      } else {
        auto parseFunc = [](const toml::table& t) {
          auto val_opt = extractToml<double>(t, "value");
          return validators::field<double>(val_opt, "value").value();
        };
        input.Jkl = parsePerSubsystemSettings<double>(*system_table, "Jkl", num_osc, ConfigDefaults::JKL, parseFunc, logger);
      }
    }

    input.rotfreq = extractScalarOrVector<double>(*system_table, "rotfreq", num_osc);

    // Parse decoherence settings
    if (system_table->contains("decoherence")) {
      auto* decoherence_table = (*system_table)["decoherence"].as_table();
      if (!decoherence_table) {
        throw validators::ValidationError("decoherence", "must be a table");
      }
      auto type_str = extractToml<std::string>(*decoherence_table, "type").value_or("none");
      input.decoherence_type = parseEnum(type_str, DECOHERENCE_TYPE_MAP, ConfigDefaults::DECOHERENCE_TYPE);
      input.decay_time = extractScalarOrVector<double>(*decoherence_table, "decay_time", num_osc);
      input.dephase_time = extractScalarOrVector<double>(*decoherence_table, "dephase_time", num_osc);
    }

    // Parse initial condition
    if (system_table->contains("initial_condition")) {
      auto* init_table = (*system_table)["initial_condition"].as_table();
      if (!init_table) {
        throw validators::ValidationError("initial_condition", "must be a table");
      }
      auto type_str = extractToml<std::string>(*init_table, "type");
      if (!type_str.has_value()) {
        throw validators::ValidationError("initial_condition.type", "field not found");
      }
      auto type_enum = parseEnum(*type_str, INITCOND_TYPE_MAP);
      if (!type_enum.has_value()) {
        throw validators::ValidationError("initial_condition.type", "unknown type: " + *type_str);
      }
      InitialConditionSettings init_cond;
      init_cond.type = type_enum.value();
      init_cond.levels = getOptionalVector<size_t>((*init_table)["levels"]);
      init_cond.filename = getOptional<std::string>((*init_table)["filename"]);
      init_cond.subsystem = getOptionalVector<size_t>((*init_table)["subsystem"]);
      input.initial_condition = init_cond;
    }

    // Hamiltonian files (inherently optional)
    input.hamiltonian_file_Hsys = getOptional<std::string>((*system_table)["hamiltonian_file_Hsys"]);
    input.hamiltonian_file_Hc = getOptional<std::string>((*system_table)["hamiltonian_file_Hc"]);

    // Control parameters
    input.control_zero_boundary_condition = extractToml<bool>(control_table, "zero_boundary_condition");

    if (control_table.contains("parameterization")) {
      input.control_parameterizations = parsePerSubsystemSettings<ControlParameterizationSettings>(
          control_table, "parameterization", num_osc, ControlParameterizationSettings{}, parseControlParamSpecsToml, logger);
    }

    if (control_table.contains("initialization")) {
      input.control_initializations = parsePerSubsystemSettings<ControlInitializationSettings>(
          control_table, "initialization", num_osc, ControlInitializationSettings{}, parseControlInitSpecsToml, logger);
    }

    input.control_amplitude_bounds = extractScalarOrVector<double>(control_table, "amplitude_bound", num_osc);

    // Parse carrier frequencies
    if (control_table.contains("carrier_frequency")) {
      auto* carrier_freq_array = control_table["carrier_frequency"].as_array();
      if (carrier_freq_array && !carrier_freq_array->empty() && !carrier_freq_array->front().is_table()) {
        // Direct array format: carrier_frequency = [1.0, 2.0]
        auto values = extractTomlVector<double>(control_table, "carrier_frequency");
        if (values.has_value()) {
          input.carrier_frequencies = std::vector<std::vector<double>>(num_osc, *values);
        }
      } else {
        // Table or array of tables format
        auto parseFunc = [](const toml::table& t) {
          auto val_opt = extractTomlVector<double>(t, "value");
          return validators::vectorField<double>(val_opt, "value").value();
        };
        std::vector<double> default_carrier_freq = {ConfigDefaults::CARRIER_FREQ};
        input.carrier_frequencies = parsePerSubsystemSettings<std::vector<double>>(
            control_table, "carrier_frequency", num_osc, default_carrier_freq, parseFunc, logger);
      }
    }

    // Optimization parameters
    if (optimization_table.contains("target")) {
      input.optim_target = parseOptimTargetToml(optimization_table, num_osc);
    }

    auto objective_str = optimization_table["objective"].value<std::string>();
    if (objective_str.has_value()) {
      input.optim_objective = parseEnum(*objective_str, OBJECTIVE_TYPE_MAP, ConfigDefaults::OPTIM_OBJECTIVE);
    }

    input.optim_weights = extractTomlVector<double>(optimization_table, "weights");
    input.optim_maxiter = extractToml<size_t>(optimization_table, "maxiter");

    // Tolerance table
    if (optimization_table.contains("tolerance")) {
      const auto* tol_table = optimization_table["tolerance"].as_table();
      if (tol_table) {
        input.optim_tol_grad_abs = extractToml<double>(*tol_table, "grad_abs");
        input.optim_tol_grad_rel = extractToml<double>(*tol_table, "grad_rel");
        input.optim_tol_finalcost = extractToml<double>(*tol_table, "final_cost");
        input.optim_tol_infidelity = extractToml<double>(*tol_table, "infidelity");
      }
    }

    // Tikhonov table
    if (optimization_table.contains("tikhonov")) {
      const auto* regul_table = optimization_table["tikhonov"].as_table();
      if (regul_table) {
        input.optim_tikhonov_coeff = extractToml<double>(*regul_table, "coeff");
        input.optim_tikhonov_use_x0 = extractToml<bool>(*regul_table, "use_x0");
      }
    }

    // Penalty table
    if (optimization_table.contains("penalty")) {
      const auto* penalty_table = optimization_table["penalty"].as_table();
      if (penalty_table) {
        input.optim_penalty_leakage = extractToml<double>(*penalty_table, "leakage");
        input.optim_penalty_weightedcost = extractToml<double>(*penalty_table, "weightedcost");
        input.optim_penalty_weightedcost_width = extractToml<double>(*penalty_table, "weightedcost_width");
        input.optim_penalty_dpdm = extractToml<double>(*penalty_table, "dpdm");
        input.optim_penalty_energy = extractToml<double>(*penalty_table, "energy");
        input.optim_penalty_variation = extractToml<double>(*penalty_table, "variation");
      }
    }

    // Output parameters
    input.output_directory = extractToml<std::string>(output_table, "directory");
    input.output_timestep_stride = extractToml<size_t>(output_table, "timestep_stride");
    input.output_optimization_stride = extractToml<size_t>(output_table, "optimization_stride");

    // Parse observables (string array to enum vector)
    if (auto observables_array = output_table["observables"].as_array()) {
      std::vector<OutputType> observables;
      for (auto&& elem : *observables_array) {
        if (auto str = elem.value<std::string>()) {
          auto enum_val = parseEnum(*str, OUTPUT_TYPE_MAP);
          if (!enum_val.has_value()) {
            throw validators::ValidationError("observables", "unknown output type: " + *str);
          }
          observables.push_back(enum_val.value());
        } else {
          throw validators::ValidationError("observables", "array must contain strings");
        }
      }
      input.output_observables = observables;
    }

    // Solver parameters
    auto runtype_str = solver_table["runtype"].value<std::string>();
    if (runtype_str.has_value()) {
      input.runtype = parseEnum(*runtype_str, RUN_TYPE_MAP, ConfigDefaults::RUNTYPE);
    }

    input.usematfree = extractToml<bool>(solver_table, "usematfree");

    // Linearsolver table
    if (solver_table.contains("linearsolver")) {
      const auto* ls_table = solver_table["linearsolver"].as_table();
      if (ls_table) {
        auto type_str = extractToml<std::string>(*ls_table, "type");
        if (type_str.has_value()) {
          input.linearsolver_type = parseEnum(*type_str, LINEAR_SOLVER_TYPE_MAP, ConfigDefaults::LINEARSOLVER_TYPE);
        }
        input.linearsolver_maxiter = extractToml<size_t>(*ls_table, "maxiter");
      }
    }

    auto timestepper_str = solver_table["timestepper"].value<std::string>();
    if (timestepper_str.has_value()) {
      input.timestepper_type = parseEnum(*timestepper_str, TIME_STEPPER_TYPE_MAP, ConfigDefaults::TIMESTEPPER_TYPE);
    }

    input.rand_seed = extractToml<int>(solver_table, "rand_seed");

    return input;
  } catch (const validators::ValidationError& e) {
    logger.exitWithError(std::string(e.what()));
  }
}

} // namespace

Config::Config(const toml::table& toml, bool quiet_mode)
    : Config(extractConfigInput(toml, MPILogger(quiet_mode)), quiet_mode) {}

Config::Config(const ConfigInput& input, bool quiet_mode) : logger(MPILogger(quiet_mode)) {
  try {
    // System parameters
    data.nlevels = validators::vectorField<size_t>(input.nlevels, "nlevels").minLength(1).positive().value();
    size_t num_osc = data.nlevels.size();

    data.nessential = validators::vectorField<size_t>(input.nessential, "nessential").valueOr(data.nlevels);

    data.ntime = validators::field<size_t>(input.ntime, "ntime").positive().value();

    data.dt = validators::field<double>(input.dt, "dt").positive().value();

    data.transfreq = validators::vectorField<double>(input.transfreq, "transfreq").hasLength(num_osc).value();

    data.selfkerr = validators::vectorField<double>(input.selfkerr, "selfkerr").valueOr(std::vector<double>(num_osc, ConfigDefaults::SELFKERR));

    size_t num_pairs = (num_osc - 1) * num_osc / 2;
    data.crosskerr = validators::vectorField<double>(input.crosskerr, "crosskerr").valueOr(std::vector<double>(num_pairs, ConfigDefaults::CROSSKERR));

    data.Jkl = validators::vectorField<double>(input.Jkl, "Jkl").valueOr(std::vector<double>(num_pairs, ConfigDefaults::JKL));

    data.rotfreq = validators::vectorField<double>(input.rotfreq, "rotfreq").valueOr(std::vector<double>(num_osc, ConfigDefaults::ROTFREQ));

    data.decoherence_type = input.decoherence_type.value_or(ConfigDefaults::DECOHERENCE_TYPE);

    data.decay_time = validators::vectorField<double>(input.decay_time, "decay_time").valueOr(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));

    data.dephase_time = validators::vectorField<double>(input.dephase_time, "dephase_time").valueOr(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));

    if (!input.initial_condition.has_value()) {
      throw validators::ValidationError("initial_condition", "field not found");
    }
    data.initial_condition = input.initial_condition.value();

    // Inherently optional fields
    data.hamiltonian_file_Hsys = input.hamiltonian_file_Hsys;
    data.hamiltonian_file_Hc = input.hamiltonian_file_Hc;

    // Control parameters
    data.control_zero_boundary_condition = input.control_zero_boundary_condition.value_or(ConfigDefaults::CONTROL_ZERO_BOUNDARY_CONDITION);

    ControlParameterizationSettings default_param;
    data.control_parameterizations = input.control_parameterizations.value_or(std::vector<ControlParameterizationSettings>(num_osc, default_param));

    ControlInitializationSettings default_init;
    data.control_initializations = input.control_initializations.value_or(std::vector<ControlInitializationSettings>(num_osc, default_init));

    data.control_amplitude_bounds = validators::vectorField<double>(input.control_amplitude_bounds, "control_amplitude_bounds").valueOr(std::vector<double>(num_osc, ConfigDefaults::CONTROL_AMPLITUDE_BOUND));

    std::vector<double> default_carrier_freq = {ConfigDefaults::CARRIER_FREQ};
    data.carrier_frequencies = input.carrier_frequencies.value_or(std::vector<std::vector<double>>(num_osc, default_carrier_freq));

    // Optimization parameters
    OptimTargetSettings default_target;
    data.optim_target = input.optim_target.value_or(default_target);

    data.optim_objective = input.optim_objective.value_or(ConfigDefaults::OPTIM_OBJECTIVE);

    data.optim_weights = validators::vectorField<double>(input.optim_weights, "optim_weights").valueOr(std::vector<double>{ConfigDefaults::OPTIM_WEIGHT});

    data.optim_tol_grad_abs = validators::field<double>(input.optim_tol_grad_abs, "optim_tol_grad_abs").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_ABS);

    data.optim_tol_grad_rel = validators::field<double>(input.optim_tol_grad_rel, "optim_tol_grad_rel").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_REL);

    data.optim_tol_finalcost = validators::field<double>(input.optim_tol_finalcost, "optim_tol_finalcost").positive().valueOr(ConfigDefaults::OPTIM_TOL_FINALCOST);

    data.optim_tol_infidelity = validators::field<double>(input.optim_tol_infidelity, "optim_tol_infidelity").positive().valueOr(ConfigDefaults::OPTIM_TOL_INFIDELITY);

    data.optim_maxiter = input.optim_maxiter.value_or(ConfigDefaults::OPTIM_MAXITER);

    data.optim_tikhonov_coeff = validators::field<double>(input.optim_tikhonov_coeff, "optim_tikhonov_coeff").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_TIKHONOV_COEFF);

    data.optim_tikhonov_use_x0 = input.optim_tikhonov_use_x0.value_or(ConfigDefaults::OPTIM_TIKHONOV_USE_X0);

    data.optim_penalty_leakage = validators::field<double>(input.optim_penalty_leakage, "optim_penalty_leakage").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);

    data.optim_penalty_weightedcost = validators::field<double>(input.optim_penalty_weightedcost, "optim_penalty_weightedcost").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);

    data.optim_penalty_weightedcost_width = validators::field<double>(input.optim_penalty_weightedcost_width, "optim_penalty_weightedcost_width").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);

    data.optim_penalty_dpdm = validators::field<double>(input.optim_penalty_dpdm, "optim_penalty_dpdm").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_DPDM);

    data.optim_penalty_energy = validators::field<double>(input.optim_penalty_energy, "optim_penalty_energy").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_ENERGY);

    data.optim_penalty_variation = validators::field<double>(input.optim_penalty_variation, "optim_penalty_variation").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY_VARIATION);

    // Output parameters
    data.output_directory = input.output_directory.value_or(ConfigDefaults::OUTPUT_DIRECTORY);

    data.output_observables = input.output_observables.value_or(std::vector<OutputType>{});

    data.output_timestep_stride = input.output_timestep_stride.value_or(ConfigDefaults::OUTPUT_TIMESTEP_STRIDE);

    data.output_optimization_stride = input.output_optimization_stride.value_or(ConfigDefaults::OUTPUT_OPTIMIZATION_STRIDE);

    // Solver parameters
    data.runtype = input.runtype.value_or(ConfigDefaults::RUNTYPE);

    data.usematfree = input.usematfree.value_or(ConfigDefaults::USEMATFREE);

    data.linearsolver_type = input.linearsolver_type.value_or(ConfigDefaults::LINEARSOLVER_TYPE);

    data.linearsolver_maxiter = validators::field<size_t>(input.linearsolver_maxiter, "linearsolver_maxiter").positive().valueOr(ConfigDefaults::LINEARSOLVER_MAXITER);

    data.timestepper_type = input.timestepper_type.value_or(ConfigDefaults::TIMESTEPPER_TYPE);

    int rand_seed_val = input.rand_seed.value_or(ConfigDefaults::RAND_SEED);
    setRandSeed(rand_seed_val);

  } catch (const validators::ValidationError& e) {
    logger.exitWithError(std::string(e.what()));
  }

  // Finalize interdependent settings, then validate
  finalize();
  validate();
}

Config::Config(const ParsedConfigData& settings, bool quiet_mode) : logger(MPILogger(quiet_mode)) {

  if (!settings.nlevels.has_value()) {
    logger.exitWithError("nlevels cannot be empty");
  }
  data.nlevels = settings.nlevels.value();
  size_t num_osc = data.nlevels.size();

  data.nessential = settings.nessential.value_or(data.nlevels);
  copyLast(data.nessential, num_osc);

  if (!settings.ntime.has_value()) {
    logger.exitWithError("ntime cannot be empty");
  }
  data.ntime = settings.ntime.value();
  if (data.ntime <= 0) {
    logger.exitWithError("ntime must be positive, got " + std::to_string(data.ntime));
  }

  if (!settings.dt.has_value()) {
    logger.exitWithError("dt cannot be empty");
  }
  data.dt = settings.dt.value();
  if (data.dt <= 0) {
    logger.exitWithError("dt must be positive, got " + std::to_string(data.dt));
  }

  if (!settings.transfreq.has_value()) {
    logger.exitWithError("transfreq cannot be empty");
  }

  data.transfreq = settings.transfreq.value();
  copyLast(data.transfreq, num_osc);

  data.selfkerr = settings.selfkerr.value_or(std::vector<double>(num_osc, ConfigDefaults::SELFKERR));
  copyLast(data.selfkerr, num_osc);

  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;
  data.crosskerr = settings.crosskerr.value_or(std::vector<double>(num_pairs_osc, ConfigDefaults::CROSSKERR));
  copyLast(data.crosskerr, num_pairs_osc);
  data.crosskerr.resize(num_pairs_osc);  // Truncate if larger than expected

  data.Jkl = settings.Jkl.value_or(std::vector<double>(num_pairs_osc, ConfigDefaults::JKL));
  copyLast(data.Jkl, num_pairs_osc);
  data.Jkl.resize(num_pairs_osc);  // Truncate if larger than expected

  data.rotfreq = settings.rotfreq.value_or(std::vector<double>(num_osc, ConfigDefaults::ROTFREQ));
  copyLast(data.rotfreq, num_osc);

  data.hamiltonian_file_Hsys = settings.hamiltonian_file_Hsys;
  data.hamiltonian_file_Hc = settings.hamiltonian_file_Hc;

  data.decoherence_type = settings.decoherence_type.value_or(ConfigDefaults::DECOHERENCE_TYPE);

  data.decay_time = settings.decay_time.value_or(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
  copyLast(data.decay_time, num_osc);

  data.dephase_time = settings.dephase_time.value_or(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
  copyLast(data.dephase_time, num_osc);

  if (!settings.initialcondition.has_value()) {
    logger.exitWithError("initialcondition cannot be empty");
  }
  data.initial_condition.type = settings.initialcondition.value().type;
  data.initial_condition.filename = settings.initialcondition.value().filename;
  data.initial_condition.levels = settings.initialcondition.value().levels;
  data.initial_condition.subsystem = settings.initialcondition.value().subsystem;

  // Control and optimization parameters
  data.control_zero_boundary_condition = settings.control_zero_boundary_condition.value_or(ConfigDefaults::CONTROL_ZERO_BOUNDARY_CONDITION);

  data.control_parameterizations = parseControlParameterizationsCfg(settings.indexed_control_parameterizations);

  // Control initialization
  if (settings.indexed_control_init.has_value()) {
    auto init_map = settings.indexed_control_init.value();
    // First check for global file initialization and populate to all oscillators if present
    if (init_map.find(0) != init_map.end() && init_map[0].filename.has_value()) {
      data.control_initializations.resize(num_osc);
      std::string control_initialization_file = init_map[0].filename.value();
      for (size_t i = 0; i < num_osc; i++) {
        data.control_initializations[i] = ControlInitializationSettings{ControlInitializationType::FILE, std::nullopt, std::nullopt, control_initialization_file};
      }
    } else {
      data.control_initializations = parseControlInitializationsCfg(settings.indexed_control_init);
    }
  } else {
    // Initialize with defaults when no control initialization is provided
    data.control_initializations.resize(num_osc);
    for (size_t i = 0; i < num_osc; i++) {
      data.control_initializations[i] = ControlInitializationSettings{
        ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt};
    }
  }

  // Parse control_amplitude_bounds from CFG format (returns vector of vectors, but we need vector)
  auto control_amplitude_bounds_cfg = parseOscillatorSettingsCfg<double>(settings.indexed_control_amplitude_bounds, data.control_parameterizations.size(), {ConfigDefaults::CONTROL_AMPLITUDE_BOUND});
  data.control_amplitude_bounds.resize(control_amplitude_bounds_cfg.size());
  for (size_t i = 0; i < control_amplitude_bounds_cfg.size(); ++i) {
    data.control_amplitude_bounds[i] = control_amplitude_bounds_cfg[i].empty() ? ConfigDefaults::CONTROL_AMPLITUDE_BOUND : control_amplitude_bounds_cfg[i][0];
  }

  data.carrier_frequencies.resize(num_osc);
  data.carrier_frequencies = parseOscillatorSettingsCfg<double>(settings.indexed_carrier_frequencies, num_osc, {ConfigDefaults::CARRIER_FREQ});

  if (settings.optim_target.has_value()) {
    data.optim_target = settings.optim_target.value();
    data.optim_target.gate_rot_freq = settings.gate_rot_freq.value_or(std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));
  } else {
    // No optim_target specified, use default (no target)
    OptimTargetSettings default_target;
    data.optim_target = default_target;
  }

  data.optim_objective = settings.optim_objective.value_or(ConfigDefaults::OPTIM_OBJECTIVE);

  data.optim_weights = settings.optim_weights.value_or(std::vector<double>{ConfigDefaults::OPTIM_WEIGHT});

  data.optim_tol_grad_abs = settings.optim_tol_grad_abs.value_or(ConfigDefaults::OPTIM_TOL_GRAD_ABS);
  data.optim_tol_grad_rel = settings.optim_tol_grad_rel.value_or(ConfigDefaults::OPTIM_TOL_GRAD_REL);
  data.optim_tol_finalcost = settings.optim_tol_finalcost.value_or(ConfigDefaults::OPTIM_TOL_FINALCOST);
  data.optim_tol_infidelity = settings.optim_tol_infidelity.value_or(ConfigDefaults::OPTIM_TOL_INFIDELITY);
  data.optim_maxiter = settings.optim_maxiter.value_or(ConfigDefaults::OPTIM_MAXITER);

  data.optim_tikhonov_coeff = settings.optim_regul.value_or(ConfigDefaults::OPTIM_TIKHONOV_COEFF);
  data.optim_tikhonov_use_x0 = settings.optim_regul_tik0.value_or(ConfigDefaults::OPTIM_TIKHONOV_USE_X0);
  if (settings.optim_regul_interpolate.has_value()) {
    // Handle deprecated optim_regul_interpolate logic
    data.optim_tikhonov_use_x0 = settings.optim_regul_interpolate.value();
    logger.log("# Warning: 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
  }

  data.optim_penalty_leakage = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
  data.optim_penalty_weightedcost = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
  data.optim_penalty_weightedcost_width = settings.optim_penalty_param.value_or(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
  data.optim_penalty_dpdm = settings.optim_penalty_dpdm.value_or(ConfigDefaults::OPTIM_PENALTY_DPDM);
  data.optim_penalty_energy = settings.optim_penalty_energy.value_or(ConfigDefaults::OPTIM_PENALTY_ENERGY);
  data.optim_penalty_variation = settings.optim_penalty_variation.value_or(ConfigDefaults::OPTIM_PENALTY_VARIATION);

  // Output parameters
  data.output_directory = settings.datadir.value_or(ConfigDefaults::OUTPUT_DIRECTORY);

  // Convert old per-oscillator output to global output_observables (apply to all oscillators)
  auto indexed_output_vec = parseOscillatorSettingsCfg<OutputType>(settings.indexed_output, num_osc);
  data.output_observables.clear();
  // Collect unique output types from all oscillators
  std::set<OutputType> unique_types;
  for (const auto& osc_output : indexed_output_vec) {
    for (const auto& type : osc_output) {
      unique_types.insert(type);
    }
  }
  // Convert set to vector
  data.output_observables.assign(unique_types.begin(), unique_types.end());

  data.output_timestep_stride = settings.output_timestep_stride.value_or(ConfigDefaults::OUTPUT_TIMESTEP_STRIDE);
  data.output_optimization_stride = settings.output_optimization_stride.value_or(ConfigDefaults::OUTPUT_OPTIMIZATION_STRIDE);
  data.runtype = settings.runtype.value_or(ConfigDefaults::RUNTYPE);
  data.usematfree = settings.usematfree.value_or(ConfigDefaults::USEMATFREE);
  data.linearsolver_type = settings.linearsolver_type.value_or(ConfigDefaults::LINEARSOLVER_TYPE);
  data.linearsolver_maxiter = settings.linearsolver_maxiter.value_or(ConfigDefaults::LINEARSOLVER_MAXITER);
  data.timestepper_type = settings.timestepper_type.value_or(ConfigDefaults::TIMESTEPPER_TYPE);
  setRandSeed(settings.rand_seed.value_or(ConfigDefaults::RAND_SEED));

  // Finalize interdependent settings, then validate
  finalize();
  validate();
}

Config Config::fromFile(const std::string& filename, bool quiet_mode) {
  if (hasSuffix(filename, ".toml")) {
    return Config::fromToml(filename, quiet_mode);
  } else {
    // TODO cfg: delete this when .cfg format is removed.
    MPILogger logger = MPILogger(quiet_mode);
    logger.log(
        "# Warning: Config file does not have .toml extension. "
        "The deprecated .cfg format will be removed in future versions.\n");
    return Config::fromCfg(filename, quiet_mode);
  }
}

Config Config::fromToml(const std::string& filename, bool quiet_mode) {
  toml::table toml = toml::parse_file(filename);
  return Config(toml, quiet_mode);
}

Config Config::fromTomlString(const std::string& toml_content, bool quiet_mode) {
  toml::table toml = toml::parse(toml_content);
  return Config(toml, quiet_mode);
}

Config Config::fromCfg(const std::string& filename, bool quiet_mode) {
  MPILogger logger = MPILogger(quiet_mode);
  CfgParser parser(logger);
  ParsedConfigData settings = parser.parseFile(filename);
  return Config(settings, quiet_mode);
}

Config Config::fromCfgString(const std::string& cfg_content, bool quiet_mode) {
  MPILogger logger = MPILogger(quiet_mode);
  CfgParser parser(logger);
  ParsedConfigData settings = parser.parseString(cfg_content);
  return Config(settings, quiet_mode);
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
  log << "nlevels = " << printVector(data.nlevels) << "\n";
  log << "nessential = " << printVector(data.nessential) << "\n";
  log << "ntime = " << data.ntime << "\n";
  log << "dt = " << data.dt << "\n";
  log << "transfreq = " << printVector(data.transfreq) << "\n";
  log << "selfkerr = " << printVector(data.selfkerr) << "\n";
  log << "crosskerr = " << toStringCoupling(data.crosskerr, data.nlevels.size()) << "\n";
  log << "Jkl = " << toStringCoupling(data.Jkl, data.nlevels.size()) << "\n";
  log << "rotfreq = " << printVector(data.rotfreq) << "\n";
  log << "decoherence = {\n";
  log << "  type = \"" << enumToString(data.decoherence_type, DECOHERENCE_TYPE_MAP) << "\",\n";
  log << "  decay_time = " << printVector(data.decay_time) << ",\n";
  log << "  dephase_time = " << printVector(data.dephase_time) << "\n";
  log << "}\n";
  log << "initial_condition = " << toString(data.initial_condition) << "\n";
  if (data.hamiltonian_file_Hsys.has_value()) {
    log << "hamiltonian_file_Hsys = \"" << data.hamiltonian_file_Hsys.value() << "\"\n";
  }
  if (data.hamiltonian_file_Hc.has_value()) {
    log << "hamiltonian_file_Hc = \"" << data.hamiltonian_file_Hc.value() << "\"\n";
  }

  log << "\n";
  log << "[control]\n";

  log << "parameterization = " << toString(data.control_parameterizations) << "\n";
  log << "carrier_frequency = " << toString(data.carrier_frequencies) << "\n";
  log << "initialization = " << toString(data.control_initializations) << "\n";
  log << "amplitude_bound = " << toString(data.control_amplitude_bounds) << "\n";
  log << "zero_boundary_condition = " << (data.control_zero_boundary_condition ? "true" : "false") << "\n";

  log << "\n";
  log << "[optimization]\n";

  log << "target = " << toString(data.optim_target) << "\n";
  log << "objective = \"" << enumToString(data.optim_objective, OBJECTIVE_TYPE_MAP) << "\"\n";
  bool uniform_weights = std::adjacent_find(data.optim_weights.begin(), data.optim_weights.end(), std::not_equal_to<double>{}) == data.optim_weights.end();
  if (!uniform_weights) {
    log << "weights = " << printVector(data.optim_weights) << "\n";
  }
  log << "tolerance = { grad_abs = " << data.optim_tol_grad_abs
      << ", grad_rel = " << data.optim_tol_grad_rel
      << ", final_cost = " << data.optim_tol_finalcost
      << ", infidelity = " << data.optim_tol_infidelity << " }\n";
  log << "maxiter = " << data.optim_maxiter << "\n";
  log << "tikhonov = { coeff = " << data.optim_tikhonov_coeff
      << ", use_x0 = " << (data.optim_tikhonov_use_x0 ? "true" : "false") << " }\n";
  log << "penalty = { leakage = " << data.optim_penalty_leakage
      << ", energy = " << data.optim_penalty_energy
      << ", dpdm = " << data.optim_penalty_dpdm
      << ", variation = " << data.optim_penalty_variation
      << ", weightedcost = " << data.optim_penalty_weightedcost
      << ", weightedcost_width = " << data.optim_penalty_weightedcost_width << " }\n";

  log << "\n";
  log << "[output]\n";

  log << "directory = \"" << data.output_directory << "\"\n";
  log << "observables = [";
  for (size_t j = 0; j < data.output_observables.size(); ++j) {
    log << "\"" << enumToString(data.output_observables[j], OUTPUT_TYPE_MAP) << "\"";
    if (j < data.output_observables.size() - 1) log << ", ";
  }
  log << "]\n";
  log << "timestep_stride = " << data.output_timestep_stride << "\n";
  log << "optimization_stride = " << data.output_optimization_stride << "\n";

  log << "\n";
  log << "[solver]\n";

  log << "runtype = \"" << enumToString(data.runtype, RUN_TYPE_MAP) << "\"\n";
  log << "usematfree = " << (data.usematfree ? "true" : "false") << "\n";
  log << "linearsolver = { type = \"" << enumToString(data.linearsolver_type, LINEAR_SOLVER_TYPE_MAP) << "\", maxiter = " << data.linearsolver_maxiter << " }\n";
  log << "timestepper = \"" << enumToString(data.timestepper_type, TIME_STEPPER_TYPE_MAP) << "\"\n";
  log << "rand_seed = " << data.rand_seed << "\n";
}

void Config::finalize() {
  // Hamiltonian file + matrix-free compatibility check
  if ((data.hamiltonian_file_Hsys.has_value() || data.hamiltonian_file_Hc.has_value()) && data.usematfree) {
    logger.log(
        "# Warning: Matrix-free solver cannot be used when Hamiltonian is read from file. Switching to sparse-matrix "
        "version.\n");
    data.usematfree = false;
  }

  if (data.usematfree && data.nlevels.size() > 5) {
    logger.log(
        "Warning: Matrix free solver is only implemented for systems with 2, 3, 4, or 5 oscillators."
        "Switching to sparse-matrix solver now.\n");
    data.usematfree = false;
  }

  // DIAGONAL and BASIS initial conditions in the Schroedinger case are the same. Overwrite it to DIAGONAL
  if (data.decoherence_type == DecoherenceType::NONE && data.initial_condition.type == InitialConditionType::BASIS) {
    data.initial_condition.type = InitialConditionType::DIAGONAL;
  }

  // For BASIS, ENSEMBLE, and DIAGONAL, or default to all oscillators IDs
  if (data.initial_condition.type == InitialConditionType::BASIS || data.initial_condition.type == InitialConditionType::ENSEMBLE || data.initial_condition.type == InitialConditionType::DIAGONAL) {
    if (!data.initial_condition.subsystem.has_value()) {
      data.initial_condition.subsystem = std::vector<size_t>(data.nlevels.size());
      for (size_t i = 0; i < data.nlevels.size(); i++) {
        data.initial_condition.subsystem->at(i) = i;
      }
    }
  }

  // Compute number of initial conditions
  n_initial_conditions = computeNumInitialConditions(data.initial_condition, data.nlevels, data.nessential, data.decoherence_type);

  // overwrite decay or dephase times with zeros, if the decoherence type is only one of them, or none.
  if (data.decoherence_type == DecoherenceType::DECAY) {
    std::fill(data.dephase_time.begin(), data.dephase_time.end(), 0);
  } else if (data.decoherence_type == DecoherenceType::DEPHASE) {
    std::fill(data.decay_time.begin(), data.decay_time.end(), 0);
  } else if (data.decoherence_type == DecoherenceType::NONE) {
    std::fill(data.decay_time.begin(), data.decay_time.end(), 0);
    std::fill(data.dephase_time.begin(), data.dephase_time.end(), 0);
  }

  // Scale optimization weights such that they sum up to one
  // If a single value was provided, replicate it for all initial conditions
  if (data.optim_weights.size() == 1) {
    // TODO remove this when removing cfg format
    copyLast(data.optim_weights, n_initial_conditions);
  } else if (data.optim_weights.size() != n_initial_conditions) {
    logger.exitWithError("optim_weights vector has length " + std::to_string(data.optim_weights.size()) + " but must have length " + std::to_string(n_initial_conditions) + " (number of initial conditions)");
  }
  // Scale the weights so that they sum up to 1
  double scaleweights = 0.0;
  for (size_t i = 0; i < data.optim_weights.size(); i++) scaleweights += data.optim_weights[i];
  if (scaleweights == 0.0) {
    logger.exitWithError("optim_weights sum to zero; at least one weight must be positive");
  }
  for (size_t i = 0; i < data.optim_weights.size(); i++) data.optim_weights[i] = data.optim_weights[i] / scaleweights;

  // Set weightedcost width to zero if weightedcost penalty is zero
  if (data.optim_penalty_weightedcost == 0.0) {
    data.optim_penalty_weightedcost_width = 0.0;
  }

  // Set control variation penalty to zero if not using 2nd order Bspline parameterization
  for (size_t i = 0; i < data.control_parameterizations.size(); i++) {
    if (data.control_parameterizations[i].type != ControlType::BSPLINE0) {
      data.optim_penalty_variation = 0.0;
      break;
    }
  }

  // Unset control initialization phase paremeter, unless BSPLINEAMP parameterization is used
  for (size_t i = 0; i < data.control_initializations.size(); i++) {
    if (data.control_parameterizations[i].type != ControlType::BSPLINEAMP) {
      data.control_initializations[i].phase = std::nullopt;
    }
  }

  // Apply defaults for control initialization amplitude
  for (size_t i = 0; i < data.control_initializations.size(); i++) {
    if (data.control_initializations[i].type != ControlInitializationType::FILE) {
      if (!data.control_initializations[i].amplitude.has_value()) {
        data.control_initializations[i].amplitude = ConfigDefaults::CONTROL_INIT_AMPLITUDE;
      }
    }
  }

  // Apply defaults and transformations for optimization target
  if (data.optim_target.type == TargetType::GATE) {
    // Prioritize gate from file: if filename is set, use FILE gate type
    if (data.optim_target.filename.has_value()) {
      data.optim_target.gate_type = GateType::FILE;
    }
    // Apply default gate_rot_freq if not set
    if (!data.optim_target.gate_rot_freq.has_value()) {
      data.optim_target.gate_rot_freq = std::vector<double>(data.nlevels.size(), ConfigDefaults::GATE_ROT_FREQ);
    }
  } else if (data.optim_target.type == TargetType::STATE) {
    // Prioritize state from file: if filename is set, clear levels
    if (data.optim_target.filename.has_value()) {
      data.optim_target.levels = std::nullopt;
    }
  }
}

void Config::validate() const {

  // Validate essential levels don't exceed total levels
  if (data.nessential.size() != data.nlevels.size()) {
    logger.exitWithError("nessential size must match nlevels size");
  }
  for (size_t i = 0; i < data.nlevels.size(); i++) {
    if (data.nessential[i] > data.nlevels[i]) {
      logger.exitWithError("nessential[" + std::to_string(i) + "] = " + std::to_string(data.nessential[i]) + " cannot exceed nlevels[" + std::to_string(i) + "] = " + std::to_string(data.nlevels[i]));
    }
  }

  // Validate control parameterization settings
  for (size_t i = 0; i < data.control_parameterizations.size(); i++) {
    const auto& param = data.control_parameterizations[i];
    if (param.type == ControlType::BSPLINE ||
        param.type == ControlType::BSPLINE0 ||
        param.type == ControlType::BSPLINEAMP) {
      if (!param.nspline.has_value()) {
        logger.exitWithError("control parameterization[" + std::to_string(i) + "] of type BSPLINE/BSPLINE0/BSPLINEAMP requires 'num'");
      }
    }
    if (param.type == ControlType::BSPLINEAMP) {
      if (!param.scaling.has_value()) {
        logger.exitWithError("control parameterization[" + std::to_string(i) + "] of type BSPLINEAMP requires 'scaling'");
      }
    }
  }

  // Validate control initialization settings
  for (size_t i = 0; i < data.control_initializations.size(); i++) {
    const auto& init = data.control_initializations[i];
    if (init.type == ControlInitializationType::FILE) {
      if (!init.filename.has_value()) {
        logger.exitWithError("control initialization[" + std::to_string(i) + "] of type FILE requires 'filename'");
      }
    }
  }

  // Validate optimization target settings
  if (data.optim_target.type == TargetType::GATE) {
    if (!data.optim_target.gate_type.has_value() && !data.optim_target.filename.has_value()) {
      logger.exitWithError("optimization target of type GATE requires 'gate_type' or 'filename'");
    }
  } else if (data.optim_target.type == TargetType::STATE) {
    if (!data.optim_target.levels.has_value() && !data.optim_target.filename.has_value()) {
      logger.exitWithError("optimization target of type STATE requires 'levels' or 'filename'");
    }
    // Validate levels size and values if provided
    if (data.optim_target.levels.has_value()) {
      if (data.optim_target.levels->size() != data.nlevels.size()) {
        logger.exitWithError("optimization target levels size (" + std::to_string(data.optim_target.levels->size()) + ") must match number of oscillators (" + std::to_string(data.nlevels.size()) + ")");
      }
      for (size_t i = 0; i < data.optim_target.levels->size(); i++) {
        if (data.optim_target.levels->at(i) >= data.nlevels[i]) {
          logger.exitWithError("optimization target levels[" + std::to_string(i) + "] = " + std::to_string(data.optim_target.levels->at(i)) + " exceeds number of modeled levels (" + std::to_string(data.nlevels[i]) + ")");
        }
      }
    }
  }

  /* Sanity check for Schrodinger solver initial conditions */
  if (data.decoherence_type == DecoherenceType::NONE) {
    if (data.initial_condition.type == InitialConditionType::ENSEMBLE ||
        data.initial_condition.type == InitialConditionType::THREESTATES ||
        data.initial_condition.type == InitialConditionType::NPLUSONE) {
      logger.exitWithError(
          "\n\n ERROR for initial condition setting: \n When running Schroedingers solver,"
          " the initial condition needs to be either 'state' or 'file' or 'diagonal' or "
          "'basis'."
          " Note that 'diagonal' and 'basis' in the Schroedinger case are the same (all unit vectors).\n\n");
    }
  }

  // Validate control bounds are positive
  for (size_t i = 0; i < data.control_amplitude_bounds.size(); i++) {
    if (data.control_amplitude_bounds[i] <= 0.0) {
      logger.exitWithError("control_amplitude_bounds[" + std::to_string(i) + "] must be positive");
    }
  }

  // Validate initial condition settings
  if (data.initial_condition.type == InitialConditionType::FROMFILE) {
    if (!data.initial_condition.filename.has_value()) {
      logger.exitWithError("initialcondition of type FROMFILE must have a filename");
    }
  }
  if (data.initial_condition.type == InitialConditionType::PRODUCT_STATE) {
    if (!data.initial_condition.levels.has_value()) {
      logger.exitWithError("initialcondition of type PRODUCT_STATE must have 'levels'");
    }
    if (data.initial_condition.levels->size() != data.nlevels.size()) {
      logger.exitWithError("initialcondition of type PRODUCT_STATE must have exactly " + std::to_string(data.nlevels.size()) + " parameters, got " + std::to_string(data.initial_condition.levels->size()));
    }
    for (size_t k = 0; k < data.initial_condition.levels->size(); k++) {
      if (data.initial_condition.levels->at(k) >= data.nlevels[k]) {
        logger.exitWithError("ERROR in config setting. The requested product state initialization " + std::to_string(data.initial_condition.levels->at(k)) + " exceeds the number of allowed levels for that oscillator (" + std::to_string(data.nlevels[k]) + ").\n");
      }
    }
  }
  if (data.initial_condition.type == InitialConditionType::BASIS ||
      data.initial_condition.type == InitialConditionType::DIAGONAL ||
      data.initial_condition.type == InitialConditionType::ENSEMBLE) {
    if (!data.initial_condition.subsystem.has_value()) {
      logger.exitWithError("initialcondition of type BASIS, DIAGONAL, or ENSEMBLE must have 'subsystem'");
    }
    if (data.initial_condition.subsystem->back() >= data.nlevels.size()) {
      logger.exitWithError("Last element in initialcondition params exceeds number of oscillators");
    }
    for (size_t i = 1; i < data.initial_condition.subsystem->size() - 1; i++) {
      if (data.initial_condition.subsystem->at(i) + 1 != data.initial_condition.subsystem->at(i + 1)) {
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
  data.rand_seed = rand_seed_;
  if (data.rand_seed < 0) {
    std::random_device rd;
    data.rand_seed = rd(); // random non-reproducable seed
  }
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
    return std::vector<ControlParameterizationSettings>(data.nlevels.size(), default_parameterization);
  }

  // Otherwise, parse specified parameterizations for each oscillator
  auto parsed_parameterizations = std::vector<ControlParameterizationSettings>(data.nlevels.size(), default_parameterization);
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

  std::vector<ControlInitializationSettings> control_initializations(data.nlevels.size(), default_init);

  if (init_configs.has_value()) {
    for (size_t i = 0; i < data.nlevels.size(); i++) {
      if (init_configs->find(static_cast<int>(i)) != init_configs->end()) {
        control_initializations[i] = init_configs->at(static_cast<int>(i));
      }
    }
  }

  return control_initializations;
}
