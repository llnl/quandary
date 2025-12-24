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
}


Config::Config(const MPILogger& logger, const toml::table& toml) : logger(logger) {
  try {
    // General options
    nlevels = validators::vectorField<size_t>(toml, "nlevels").minLength(1).positive().value();

    size_t num_osc = nlevels.size();

    nessential = validators::vectorField<size_t>(toml, "nessential").minLength(1).positive().valueOr(nlevels);
    copyLast(nessential, num_osc);

    ntime = toml["ntime"].value_or(ConfigDefaults::NTIME);

    dt = toml["dt"].value_or(ConfigDefaults::DT);

    transfreq = validators::vectorField<double>(toml, "transfreq").minLength(1).value();
    copyLast(transfreq, num_osc);

    selfkerr = validators::vectorField<double>(toml, "selfkerr")
                   .minLength(1)
                   .valueOr(std::vector<double>(num_osc, ConfigDefaults::SELFKERR));
    copyLast(selfkerr, num_osc);

    crosskerr = parseCouplingParameters(toml, "crosskerr", num_osc, ConfigDefaults::CROSSKERR);
    Jkl = parseCouplingParameters(toml, "Jkl", num_osc, ConfigDefaults::JKL);

    rotfreq = validators::vectorField<double>(toml, "rotfreq").minLength(1).value();
    copyLast(rotfreq, num_osc);

    collapse_type = ConfigDefaults::COLLAPSE_TYPE_ENUM;
    decay_time = std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME);
    dephase_time = std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME);
    if (toml.contains("decoherence")) {
      auto* decoherence_table = toml["decoherence"].as_table();
      if (decoherence_table) {
        auto type_str = validators::field<std::string>(*decoherence_table, "type")
                            .valueOr("none");
        collapse_type = parseEnum(type_str, LINDBLAD_TYPE_MAP, ConfigDefaults::COLLAPSE_TYPE_ENUM);
        decay_time = validators::vectorField<double>(*decoherence_table, "decay_time")
                         .valueOr(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
        dephase_time = validators::vectorField<double>(*decoherence_table, "dephase_time")
                           .valueOr(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
      }
    }

    auto init_cond_table = validators::getRequiredTable(toml, "initial_condition");
    auto type_opt = parseEnum(validators::field<std::string>(init_cond_table, "type").value(), INITCOND_TYPE_MAP);
    std::optional<std::vector<size_t>> levels = validators::getOptionalVector<size_t>(init_cond_table["levels"]);
    std::optional<std::vector<size_t>> osc_IDs = validators::getOptionalVector<size_t>(init_cond_table["oscIDs"]);
    std::optional<std::string> filename = init_cond_table["filename"].value<std::string>();
    initial_condition = parseInitialCondition(type_opt, filename, levels, osc_IDs);
    n_initial_conditions = computeNumInitialConditions();

    apply_pipulse = std::vector<std::vector<PiPulseSegment>>(nlevels.size());
    auto apply_pipulse_array_of_tables = validators::getArrayOfTables(toml, "apply_pipulse");
    for (auto& elem : apply_pipulse_array_of_tables) {
      auto pipulse_table = *elem.as_table();

      // Validate allowed keys in apply_pipulse table
      const std::string tstart_key = "tstart";
      const std::string tstop_key = "tstop";
      const std::string amp_key = "amp";
      const std::set<std::string> allowed_keys = {OSC_ID_KEY, tstart_key, tstop_key, amp_key};
      validateTableKeys(pipulse_table, allowed_keys, "apply_pipulse");

      size_t oscilID = validators::field<size_t>(pipulse_table, OSC_ID_KEY).value();
      double tstart = validators::field<double>(pipulse_table, tstart_key).value();
      double tstop = validators::field<double>(pipulse_table, tstop_key).value();
      double amp = validators::field<double>(pipulse_table, amp_key).value();

      addPiPulseSegment(apply_pipulse, oscilID, tstart, tstop, amp);
    }

    hamiltonian_file_Hsys = toml["hamiltonian_file_Hsys"].value<std::string>();
    hamiltonian_file_Hc = toml["hamiltonian_file_Hc"].value<std::string>();

    // Optimization options
    control_enforceBC = toml["control_enforceBC"].value_or(ConfigDefaults::CONTROL_ENFORCE_BC);

    // Parse control parameterizations
    auto control_array = validators::getArrayOfTables(toml, "control_parameterizations");
    control_parameterizations = parseControlParameterizations(control_array, num_osc);

    // Parse control initialization: table with optional default and per-subsystem overrides
    if (toml.contains("control_initialization")) {
      if (toml["control_initialization"].is_table()) {
        auto* control_init_table = toml["control_initialization"].as_table();
        control_initializations = parseControlInitializations(*control_init_table, num_osc);
      } else {
        logger.exitWithError("control_initialization must be a table");
      }
    } else {
      // No control_initialization specified, use defaults
      ControlInitialization default_init = ControlInitialization{ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt}; 
      control_initializations.resize(num_osc, default_init);
    }

    // Parse control bounds
    auto control_bounds_array = validators::getArrayOfTables(toml, "control_bounds");

    // Allowed keys for control_bounds table
    const std::string values_key = "values";
    const std::set<std::string> control_bounds_allowed_keys = {OSC_ID_KEY, values_key};

    control_bounds = parseOscillatorSettings<double>(control_bounds_array, num_osc, {ConfigDefaults::CONTROL_BOUND},
                                                     values_key, control_bounds_allowed_keys, "control_bounds");

    // Parse carrier frequencies
    auto carrier_freq_array = validators::getArrayOfTables(toml, "carrier_frequency");

    // Allowed keys for carrier_frequency table
    const std::string carrier_values_key = "values";
    const std::set<std::string> carrier_freq_allowed_keys = {OSC_ID_KEY, carrier_values_key};

    carrier_frequencies =
        parseOscillatorSettings<double>(carrier_freq_array, num_osc, {ConfigDefaults::CARRIER_FREQ}, carrier_values_key,
                                        carrier_freq_allowed_keys, "carrier_frequency");

    // optim_target
    if (toml.contains("optim_target")) {
      auto target_table = *toml["optim_target"].as_table();
      std::string type_str = validators::field<std::string>(target_table, "target_type").value();
      auto target_type_opt = parseEnum(type_str, TARGET_TYPE_MAP);

      if (!target_type_opt.has_value()) {
        logger.exitWithError("Unknown optimization target type: " + type_str);
      }

      std::optional<GateType> gate_type_opt;
      std::optional<std::string> gate_type_str = target_table["gate_type"].value<std::string>();
      if (gate_type_str.has_value()) {
        gate_type_opt = parseEnum(gate_type_str.value(), GATE_TYPE_MAP);
      }

      auto gate_file_opt = target_table["gate_file"].value<std::string>();
      auto target_levels = validators::getOptionalVector<size_t>(target_table["levels"]);
      auto target_filename = target_table["filename"].value<std::string>();

      optim_target =
          parseOptimTarget(target_type_opt.value(), gate_type_opt, gate_file_opt, target_levels, target_filename);
    } else {
      optim_target =
          parseOptimTarget(ConfigDefaults::OPTIM_TARGET, std::nullopt, std::nullopt, std::nullopt, std::nullopt);
    }

    gate_rot_freq = validators::vectorField<double>(toml, "gate_rot_freq")
                        .valueOr(std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));
    copyLast(gate_rot_freq, num_osc);

    optim_objective =
        parseEnum(toml["optim_objective"].value<std::string>(), OBJECTIVE_TYPE_MAP, ConfigDefaults::OPTIM_OBJECTIVE);

    std::optional<std::vector<double>> optim_weights_opt = validators::getOptionalVector<double>(toml["optim_weights"]);
    optim_weights = optim_weights_opt.value_or(std::vector<double>{ConfigDefaults::OPTIM_WEIGHT});
    copyLast(optim_weights, n_initial_conditions);
    optim_weights.resize(n_initial_conditions);

    // Optimization tolerances: use nested table `optim_tolerance` only; default all if missing
    if (toml.contains("optim_tolerance")) {
      auto* tol_table = toml["optim_tolerance"].as_table();
      if (!tol_table) {
        logger.exitWithError("optim_tolerance must be a table");
      }
      optim_tol_grad_abs = validators::field<double>(*tol_table, "grad_abs").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_ABS);
      optim_tol_grad_rel = validators::field<double>(*tol_table, "grad_rel").positive().valueOr(ConfigDefaults::OPTIM_TOL_GRAD_REL);
      optim_tol_finalcost = validators::field<double>(*tol_table, "final_cost").positive().valueOr(ConfigDefaults::OPTIM_TOL_FINALCOST);
      optim_tol_infidelity = validators::field<double>(*tol_table, "infidelity").positive().valueOr(ConfigDefaults::OPTIM_TOL_INFIDELITY);
    } else {
      optim_tol_grad_abs = ConfigDefaults::OPTIM_TOL_GRAD_ABS;
      optim_tol_grad_rel = ConfigDefaults::OPTIM_TOL_GRAD_REL;
      optim_tol_finalcost = ConfigDefaults::OPTIM_TOL_FINALCOST;
      optim_tol_infidelity = ConfigDefaults::OPTIM_TOL_INFIDELITY;
    }
    optim_maxiter = validators::field<size_t>(toml, "optim_maxiter").positive().valueOr(ConfigDefaults::OPTIM_MAXITER);

    // Parse optim_tikhonov inline table
    if (toml.contains("optim_tikhonov")) {
      auto regul_table = toml["optim_tikhonov"].as_table();
      if (!regul_table) {
        logger.exitWithError("optim_tikhonov must be a table");
      }
      optim_tikhonov_coeff= validators::field<double>(*regul_table, "coeff")
          .greaterThanEqual(0.0)
          .valueOr(ConfigDefaults::OPTIM_TIKHONOV_COEFF);
      optim_tikhonov_use_x0 = validators::field<bool>(*regul_table, "use_x0")
          .valueOr(ConfigDefaults::OPTIM_TIKHONOV_USE_X0);
    } else {
      optim_tikhonov_coeff = ConfigDefaults::OPTIM_TIKHONOV_COEFF;
      optim_tikhonov_use_x0 = ConfigDefaults::OPTIM_TIKHONOV_USE_X0;
    }

    // Parse optim_penalty inline table
    if (toml.contains("optim_penalty")) {
      auto penalty_table = toml["optim_penalty"].as_table();
      if (!penalty_table) {
        logger.exitWithError("optim_penalty must be a table");
      }
      optim_penalty_leakage = validators::field<double>(*penalty_table, "leakage")
          .greaterThanEqual(0.0)
          .valueOr(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
      optim_penalty_weightedcost = validators::field<double>(*penalty_table, "weightedcost")
          .greaterThanEqual(0.0)
          .valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
      optim_penalty_weightedcost_width = validators::field<double>(*penalty_table, "weightedcost_width")
          .greaterThanEqual(0.0)
          .valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
      optim_penalty_dpdm = validators::field<double>(*penalty_table, "dpdm")
                               .greaterThanEqual(0.0)
                               .valueOr(ConfigDefaults::OPTIM_PENALTY_DPDM);
      optim_penalty_energy = validators::field<double>(*penalty_table, "energy")
                                 .greaterThanEqual(0.0)
                                 .valueOr(ConfigDefaults::OPTIM_PENALTY_ENERGY);
      optim_penalty_variation = validators::field<double>(*penalty_table, "variation")
                                    .greaterThanEqual(0.0)
                                    .valueOr(ConfigDefaults::OPTIM_PENALTY_VARIATION);
    } else {
      // Set defaults of optim_penalty table is not given 
      optim_penalty_leakage = ConfigDefaults::OPTIM_PENALTY_LEAKAGE;
      optim_penalty_weightedcost = ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST;
      optim_penalty_weightedcost_width = ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH;
      optim_penalty_dpdm = ConfigDefaults::OPTIM_PENALTY_DPDM;
      optim_penalty_energy = ConfigDefaults::OPTIM_PENALTY_ENERGY;
      optim_penalty_variation = ConfigDefaults::OPTIM_PENALTY_VARIATION;
    }

    datadir = toml["datadir"].value_or(ConfigDefaults::DATADIR);

    // Parse output_type as an array of strings
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

    output_frequency = toml["output_frequency"].value_or(ConfigDefaults::OUTPUT_FREQUENCY);
    optim_monitor_frequency = toml["optim_monitor_frequency"].value_or(ConfigDefaults::OPTIM_MONITOR_FREQUENCY);

    runtype = parseEnum(toml["runtype"].value<std::string>(), RUN_TYPE_MAP, ConfigDefaults::RUNTYPE);

    usematfree = toml["usematfree"].value_or(ConfigDefaults::USEMATFREE);

    linearsolver_type = parseEnum(toml["linearsolver_type"].value<std::string>(), LINEAR_SOLVER_TYPE_MAP,
                                  ConfigDefaults::LINEARSOLVER_TYPE);

    linearsolver_maxiter = toml["linearsolver_maxiter"].value_or(ConfigDefaults::LINEARSOLVER_MAXITER);

    timestepper_type =
        parseEnum(toml["timestepper"].value<std::string>(), TIME_STEPPER_TYPE_MAP, ConfigDefaults::TIMESTEPPER_TYPE);

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
  size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;

  nessential = settings.nessential.value_or(nlevels);
  copyLast(nessential, num_osc);

  ntime = settings.ntime.value_or(ConfigDefaults::NTIME);

  dt = settings.dt.value_or(ConfigDefaults::DT);

  if (!settings.transfreq.has_value()) {
    logger.exitWithError("transfreq cannot be empty");
  }
  transfreq = settings.transfreq.value();
  copyLast(transfreq, num_osc);

  selfkerr = settings.selfkerr.value_or(std::vector<double>(num_osc, ConfigDefaults::SELFKERR));
  copyLast(selfkerr, num_osc);

  crosskerr = settings.crosskerr.value_or(std::vector<double>(num_pairs_osc, ConfigDefaults::CROSSKERR));
  copyLast(crosskerr, num_pairs_osc);

  Jkl = settings.Jkl.value_or(std::vector<double>(num_pairs_osc, ConfigDefaults::JKL));
  copyLast(Jkl, num_pairs_osc);

  if (!settings.rotfreq.has_value()) {
    logger.exitWithError("rotfreq cannot be empty");
  }
  rotfreq = settings.rotfreq.value();
  copyLast(rotfreq, num_osc);

  collapse_type = settings.collapse_type.value_or(ConfigDefaults::COLLAPSE_TYPE_ENUM);

  decay_time = settings.decay_time.value_or(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
  copyLast(decay_time, num_osc);

  dephase_time = settings.dephase_time.value_or(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
  copyLast(dephase_time, num_osc);

  if (!settings.initialcondition.has_value()) {
    logger.exitWithError("initialcondition cannot be empty");
  }
  initial_condition =
      parseInitialCondition(settings.initialcondition.value().type, settings.initialcondition.value().filename,
                            settings.initialcondition.value().levels, settings.initialcondition.value().osc_IDs);
  n_initial_conditions = computeNumInitialConditions();

  apply_pipulse = std::vector<std::vector<PiPulseSegment>>(nlevels.size());
  if (settings.apply_pipulse.has_value()) {
    for (const auto& pulse_config : settings.apply_pipulse.value()) {
      addPiPulseSegment(apply_pipulse, pulse_config.oscil_id, pulse_config.tstart, pulse_config.tstop,
                        pulse_config.amp);
    }
  }

  hamiltonian_file_Hsys = settings.hamiltonian_file_Hsys;
  hamiltonian_file_Hc = settings.hamiltonian_file_Hc;

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
        control_initializations[i] = ControlInitialization{ControlInitializationType::FILE, std::nullopt, std::nullopt, control_initialization_file};
      }
    } else {
      control_initializations = parseControlInitializationsCfg(settings.indexed_control_init);
    }
  } else {
    // Initialize with defaults when no control initialization is provided
    control_initializations.resize(num_osc);
    for (size_t i = 0; i < num_osc; i++) {
      control_initializations[i] = ControlInitialization{
        ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt};
    }
  }

  control_bounds = parseOscillatorSettingsCfg<double>(settings.indexed_control_bounds, control_parameterizations.size(),
                                                      {ConfigDefaults::CONTROL_BOUND});

  carrier_frequencies.resize(num_osc);
  carrier_frequencies = parseOscillatorSettingsCfg<double>(settings.indexed_carrier_frequencies, num_osc, {ConfigDefaults::CARRIER_FREQ});
  if (settings.optim_target.has_value()) {
    const OptimTargetSettings& target_config = settings.optim_target.value();
    optim_target = parseOptimTarget(target_config.type, target_config.gate_type, target_config.gate_file,
                                    target_config.levels, target_config.file);
  } else {
    optim_target =
        parseOptimTarget(ConfigDefaults::OPTIM_TARGET, std::nullopt, std::nullopt, std::nullopt, std::nullopt);
  }

  gate_rot_freq = settings.gate_rot_freq.value_or(std::vector<double>(num_osc, ConfigDefaults::GATE_ROT_FREQ));
  copyLast(gate_rot_freq, num_osc);

  optim_objective = settings.optim_objective.value_or(ConfigDefaults::OPTIM_OBJECTIVE);

  optim_weights = settings.optim_weights.value_or(std::vector<double>{ConfigDefaults::OPTIM_WEIGHT});
  copyLast(optim_weights, n_initial_conditions);
  optim_weights.resize(n_initial_conditions);

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
  
  output_frequency = settings.output_frequency.value_or(ConfigDefaults::OUTPUT_FREQUENCY);
  optim_monitor_frequency = settings.optim_monitor_frequency.value_or(ConfigDefaults::OPTIM_MONITOR_FREQUENCY);
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

std::vector<double> Config::parseCouplingParameters(const toml::table& toml, const std::string& key, size_t num_osc, double default_value) const {
  size_t num_pairs = (num_osc - 1) * num_osc / 2;
  std::vector<double> result(num_pairs, default_value);

  // If the key doesn't exist, return default vector, otherwise grab the table
  if (!toml.contains(key)) {
    return result;
  }
  auto* table = toml[key].as_table();

  // If it's not a table, fall back to the old vector format for backwards compatibility
  if (!table) {
    auto* arr = toml[key].as_array();
    if (arr) {
      for (size_t i = 0; i < arr->size() && i < result.size(); i++) {
        auto val = arr->at(i).value<double>();
        if (val) {
          result[i] = *val;
        }
      }
    }
    return result;
  }

  // Parse table format: keys are "i-j" pairs
  for (auto& [pair_key, value_node] : *table) { // This iterates over key-value pairs in the table

    std::string key_str(pair_key.str()); 
    size_t dash_pos = key_str.find('-'); 
    if (dash_pos == std::string::npos) {
      throw validators::ValidationError(key, "coupling key must be in format 'i-j' (e.g., '0-1')");
    }
    size_t i = std::stoul(key_str.substr(0, dash_pos));
    size_t j = std::stoul(key_str.substr(dash_pos + 1));

    // Ensure i < j and i,j < num_oscillators
    if (i >= num_osc || j >= num_osc) {
      throw validators::ValidationError(key, "oscillator index out of range for key '" + key_str + "'");
    }
    if (i > j) {
      std::swap(i, j);
    }

    // Convert (i,j) pair to linear index: (0,1), (0,2), ..., (0,n-1), (1,2), ..., (1,n-1), ..., (n-2,n-1)
    // Offset for first element i: i*num_osc - i - i*(i-1)/2
    // Index within pairs starting with i: j - i - 1
    size_t pair_index = i * num_osc - i - i * (i - 1) / 2 + (j - i - 1);
    auto coupling_val = value_node.value<double>();
    if (!coupling_val) {
      throw validators::ValidationError(key, "coupling value for key '" + key_str + "' must be a number");
    }
    result[pair_index] = *coupling_val;
  }

  return result;
}

std::string Config::printCouplingParameters(const std::vector<double>& couplings, size_t num_osc) const {
  if (couplings.empty()) return "{}";

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
  std::string result = "{ ";
  for (size_t i = 0; i < nonzero_couplings.size(); ++i) {
    auto [pair, value] = nonzero_couplings[i];
    auto [first, second] = pair;
    result += "\"" + std::to_string(first) + "-" + std::to_string(second) + "\" = " + std::to_string(value);
    if (i < nonzero_couplings.size() - 1) {
      result += ", ";
    }
  }
  result += " }";
  return result;
}

namespace {

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


std::string toString(const InitialCondition& initial_condition) {
  auto type_str = "type = \"" + enumToString(initial_condition.type, INITCOND_TYPE_MAP) + "\"";
  switch (initial_condition.type) {
    case InitialConditionType::FROMFILE:
      return "{" + type_str + ", filename = \"" + initial_condition.filename.value() + "\"}";
    case InitialConditionType::PURE: {
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
  auto type_str = "target_type = \"" + enumToString(optim_target.type, TARGET_TYPE_MAP) + "\"";
  switch (optim_target.type) {
    case TargetType::GATE: {
      std::string out = "{" + type_str;
      if (optim_target.gate_type.has_value()) {
        out += ", gate_type = \"" + enumToString(optim_target.gate_type.value(), GATE_TYPE_MAP) + "\"";
      }
      if (optim_target.gate_file.has_value() && !optim_target.gate_file.value().empty()) {
        out += ", gate_file = \"" + optim_target.gate_file.value() + "\"";
      }
      out += "}";
      return out;
    }
    case TargetType::PURE: {
      std::string out = "{" + type_str;
      if (optim_target.levels.has_value()) {
        out += ", levels = " + printVector(optim_target.levels.value());
      }
      out += "}";
      return out;
    }
    case TargetType::FROMFILE: {
      std::string out = "{" + type_str;
      if (optim_target.file.has_value()) {
        out += ", file = \"" + optim_target.file.value() + "\"";
      }
      out += "}";
      return out;
    }
  }
  return "unknown";
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
  log << "crosskerr = " << printCouplingParameters(crosskerr, nlevels.size()) << "\n";
  log << "Jkl = " << printCouplingParameters(Jkl, nlevels.size()) << "\n";
  log << "rotfreq = " << printVector(rotfreq) << "\n";
  log << "decoherence = {\n";
  log << "  type = \"" << enumToString(collapse_type, LINDBLAD_TYPE_MAP) << "\",\n";
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

  // Optimization parameters
  log << "control_enforceBC = " << (control_enforceBC ? "true" : "false") << "\n";
  log << "optim_target = " << toString(optim_target) << "\n";
  log << "gate_rot_freq = " << printVector(gate_rot_freq) << "\n";
  log << "optim_objective = \"" << enumToString(optim_objective, OBJECTIVE_TYPE_MAP) << "\"\n";
  log << "optim_weights = " << printVector(optim_weights) << "\n";
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

  // Output parameters
  log << "datadir = \"" << datadir << "\"\n";
  log << "output_frequency = " << output_frequency << "\n";
  log << "output_type = [";
  for (size_t j = 0; j < output_type.size(); ++j) {
    log << "\"" << enumToString(output_type[j], OUTPUT_TYPE_MAP) << "\"";
    if (j < output_type.size() - 1) log << ", ";
  }
  log << "]\n\n";
  log << "optim_monitor_frequency = " << optim_monitor_frequency << "\n";
  log << "runtype = \"" << enumToString(runtype, RUN_TYPE_MAP) << "\"\n";
  log << "usematfree = " << (usematfree ? "true" : "false") << "\n";
  log << "linearsolver_type = \"" << enumToString(linearsolver_type, LINEAR_SOLVER_TYPE_MAP) << "\"\n";
  log << "linearsolver_maxiter = " << linearsolver_maxiter << "\n";
  log << "timestepper = \"" << enumToString(timestepper_type, TIME_STEPPER_TYPE_MAP) << "\"\n";
  log << "rand_seed = " << rand_seed << "\n";

  // Control initialization 

  // Check if all subsystems have the same initialization, in which case a simple table is written instead of per-subsystem overrides
  bool all_same = true;
  const auto& first_init = control_initializations[0];
  for (size_t i = 1; i < control_initializations.size(); ++i) {
    if (control_initializations[i].type != first_init.type ||
        control_initializations[i].amplitude != first_init.amplitude ||
        control_initializations[i].phase != first_init.phase) {
      all_same = false;
      break;
    }
  }
  
  log << "control_initialization = {\n";
  if (all_same) {
    // All subsystems have the same initialization - output as simple table 
    log << "type = \"" << enumToString(first_init.type, CONTROL_INITIALIZATION_TYPE_MAP) << "\"";
    if (first_init.filename.has_value()) log << ", filename = \"" << first_init.filename.value() << "\"";
    if (first_init.amplitude.has_value()) log << ", amplitude = " << first_init.amplitude.value();
    if (first_init.phase.has_value())  log << ", phase = " << first_init.phase.value();
    log << " }\n\n";
  } else {
    // Subsystems have different initializations - output with per-subsystem overrides
    for (size_t i = 0; i < control_initializations.size(); ++i) {
      const auto& init = control_initializations[i];
      log << "  \"" << i << "\" = { ";
      log << "type = \"" << enumToString(init.type, CONTROL_INITIALIZATION_TYPE_MAP) << "\"";
      if (init.filename.has_value()) log << ", filename = \"" << init.filename.value() << "\"";
      if (init.amplitude.has_value()) log << ", amplitude = " << init.amplitude.value();
      if (init.phase.has_value()) log << ", phase = " << init.phase.value();
      log << " }";
      if (i < control_initializations.size() - 1) {
        log << ",";
      }
      log << "\n";
    }
    log << "}\n\n";
  }

  // Section 2: All array-of-tables at the end
  log << "\n";

  // Apply pi-pulse array of tables
  for (size_t i = 0; i < apply_pipulse.size(); ++i) {
    for (const auto& parameterization : apply_pipulse[i]) {
      log << "[[apply_pipulse]]\n";
      log << "oscID = " << i << "\n";
      log << "tstart = " << parameterization.tstart << "\n";
      log << "tstop = " << parameterization.tstop << "\n";
      log << "amp = " << parameterization.amp << "\n";
      log << "\n";
    }
  }

  // Control parameterizations as array of tables
  for (size_t i = 0; i < control_parameterizations.size(); ++i) {
      const auto& seg = control_parameterizations[i];
      log << "[[control_parameterizations]]\n";
      log << "oscID = " << i << "\n";
      log << "type = \"" << enumToString(seg.type, CONTROL_TYPE_MAP) << "\"\n";

      // Add parameterization-specific parameters
      if (seg.type == ControlType::BSPLINE || seg.type == ControlType::BSPLINE0) {
        log << "num = " << seg.nspline.value() << "\n";
        log << "tstart = " << seg.tstart.value() << "\n";
        log << "tstop = " << seg.tstop.value() << "\n";
      } else if (seg.type == ControlType::BSPLINEAMP) {
        log << "num = " << seg.nspline.value() << "\n";
        log << "scaling = " << seg.scaling.value() << "\n";
        log << "tstart = " << seg.tstart.value() << "\n";
        log << "tstop = " << seg.tstop.value() << "\n";
      }
      log << "\n";
  }

  // Control bounds as array of tables
  for (size_t i = 0; i < control_bounds.size(); ++i) {
    if (!control_bounds[i].empty()) {
      log << "[[control_bounds]]\n";
      log << "oscID = " << i << "\n";
      log << "values = " << printVector(control_bounds[i]) << "\n\n";
    }
  }

  // Carrier frequencies as array of tables
  for (size_t i = 0; i < carrier_frequencies.size(); ++i) {
    if (!carrier_frequencies[i].empty()) {
      log << "[[carrier_frequency]]\n";
      log << "oscID = " << i << "\n";
      log << "values = " << printVector(carrier_frequencies[i]) << "\n\n";
    }
  }
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

  if (collapse_type == LindbladType::NONE && initial_condition.type == InitialConditionType::BASIS) {
    // DIAGONAL and BASIS initial conditions in the Schroedinger case are the same. Overwrite it to DIAGONAL
    initial_condition.type = InitialConditionType::DIAGONAL;
  }

  // overwrite decay or dephase times with zeros, if the collapse type is only one of them.
  if (collapse_type == LindbladType::DECAY) {
    std::fill(dephase_time.begin(), dephase_time.end(), 0);
  } else if (collapse_type == LindbladType::DEPHASE) {
    std::fill(decay_time.begin(), decay_time.end(), 0);
  }

  // Scale optimization weights such that they sum up to one
  double scaleweights = 0.0;
  assert(n_initial_conditions == optim_weights.size());
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
  if (ntime <= 0) {
    logger.exitWithError("ntime must be positive, got " + std::to_string(ntime));
  }

  if (dt <= 0) {
    logger.exitWithError("dt must be positive, got " + std::to_string(dt));
  }

  // Validate essential levels don't exceed total levels
  if (nessential.size() != nlevels.size()) {
    logger.exitWithError("nessential size must match nlevels size");
  }

  for (size_t i = 0; i < nlevels.size(); i++) {
    if (nessential[i] > nlevels[i]) {
      logger.exitWithError("nessential[" + std::to_string(i) + "] = " + std::to_string(nessential[i]) +
                           " cannot exceed nlevels[" + std::to_string(i) + "] = " + std::to_string(nlevels[i]));
    }
  }

  /* Sanity check for Schrodinger solver initial conditions */
  if (collapse_type == LindbladType::NONE) {
    if (initial_condition.type == InitialConditionType::ENSEMBLE ||
        initial_condition.type == InitialConditionType::THREESTATES ||
        initial_condition.type == InitialConditionType::NPLUSONE) {
      logger.exitWithError(
          "\n\n ERROR for initial condition setting: \n When running Schroedingers solver"
          " (collapse_type == NONE), the initial condition needs to be either 'pure' or 'from file' or 'diagonal' or "
          "'basis'."
          " Note that 'diagonal' and 'basis' in the Schroedinger case are the same (all unit vectors).\n\n");
    }
  }
}

size_t Config::computeNumInitialConditions() const {
  size_t n_initial_conditions = 0;
  switch (initial_condition.type) {
    case InitialConditionType::FROMFILE:
    case InitialConditionType::PURE:
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
      if (!initial_condition.osc_IDs.has_value()) {
        logger.exitWithError("expected diagonal initial condition to have list of oscIDs");
      }
      n_initial_conditions = 1;
      for (size_t oscilID : initial_condition.osc_IDs.value()) {
        if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
      }
      break;
    case InitialConditionType::BASIS:
      /* Compute ninit = dim(subsystem defined by list of oscil IDs) */
      if (!initial_condition.osc_IDs.has_value()) {
        logger.exitWithError("expected diagonal initial condition to have list of oscIDs");
      }
      n_initial_conditions = 1;
      for (size_t oscilID : initial_condition.osc_IDs.value()) {
        if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
      }
      // if Schroedinger solver: ninit = N, do nothing.
      // else Lindblad solver: ninit = N^2
      if (collapse_type != LindbladType::NONE) {
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

void Config::validateTableKeys(const toml::table& table, const std::set<std::string>& allowed_keys,
                               const std::string& table_name) const {
  for (const auto& [key, _] : table) {
    if (allowed_keys.find(std::string(key.str())) == allowed_keys.end()) {
      logger.exitWithError("Unknown key '" + std::string(key.str()) + "' in " + table_name + ".");
    }
  }
}

template <typename T>
std::vector<std::vector<T>> Config::parseOscillatorSettings(const toml::array& array_of_tables, size_t num_entries,
                                                            std::vector<T> default_values,
                                                            const std::string& field_name,
                                                            const std::set<std::string>& allowed_keys,
                                                            const std::string& table_name) const {
  std::vector<std::vector<T>> result(num_entries, default_values);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    validateTableKeys(table, allowed_keys, table_name);

    std::vector<T> values = validators::vectorField<T>(table, field_name).value();

    if (auto osc_id_node = table[OSC_ID_KEY]) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, OSC_ID_KEY).lessThan(num_entries).value();
      result[osc_id] = values;
    } else {
      // Apply to ALL oscillators
      for (size_t i = 0; i < num_entries; i++) {
        result[i] = values;
      }
    }
  }

  return result;
}

InitialCondition Config::parseInitialCondition(std::optional<InitialConditionType> opt_type,
                                               const std::optional<std::string>& filename,
                                               const std::optional<std::vector<size_t>>& levels,
                                               const std::optional<std::vector<size_t>>& osc_IDs) const {
  if (!opt_type.has_value()) {
    logger.exitWithError("initial condition type not found.");
  }
  InitialConditionType type = opt_type.value();

  // If no params are given for BASIS, ENSEMBLE, or DIAGONAL, default to all oscillators
  auto init_cond_IDs = osc_IDs.value_or(std::vector<size_t>{});
  if (!osc_IDs.has_value() &&
      (type == InitialConditionType::BASIS || type == InitialConditionType::ENSEMBLE ||
       type == InitialConditionType::DIAGONAL)) {
    for (size_t i = 0; i < nlevels.size(); i++) {
      init_cond_IDs.push_back(i);
    }
  }

  InitialCondition result;
  result.type = type;

  switch (type) {
    case InitialConditionType::FROMFILE:
      if (!filename.has_value()) {
        logger.exitWithError("initialcondition of type FROMFILE must have a filename");
      }
      result.filename = filename.value();
      break;
    case InitialConditionType::PURE:
      if (!levels.has_value()) {
        logger.exitWithError("initialcondition of type PURE must have 'levels'");
      }
      if (levels.value().size() != nlevels.size()) {
        logger.exitWithError("initialcondition of type PURE must have exactly " + std::to_string(nlevels.size()) +
                             " parameters, got " + std::to_string(levels.value().size()));
      }
      for (size_t k = 0; k < levels.value().size(); k++) {
        if (levels.value()[k] >= nlevels[k]) {
          logger.exitWithError(
              "ERROR in config setting. The requested pure state initialization " + std::to_string(levels.value()[k]) +
              " exceeds the number of allowed levels for that oscillator (" + std::to_string(nlevels[k]) + ").\n");
        }
      }
      result.levels = levels.value();
      break;

    case InitialConditionType::BASIS:
    case InitialConditionType::DIAGONAL:
      result.osc_IDs = init_cond_IDs;
      break;

    case InitialConditionType::ENSEMBLE:
      if (init_cond_IDs.back() >= nlevels.size()) {
        logger.exitWithError("Last element in initialcondition params exceeds number of oscillators");
      }

      for (size_t i = 1; i < init_cond_IDs.size() - 1; i++) {
        if (init_cond_IDs[i] + 1 != init_cond_IDs[i + 1]) {
          logger.exitWithError("List of oscillators for ensemble initialization should be consecutive!\n");
        }
      }
      result.osc_IDs = init_cond_IDs;
      break;

    case InitialConditionType::THREESTATES:
    case InitialConditionType::NPLUSONE:
    case InitialConditionType::PERFORMANCE:
      // No additional fields needed for these types
      break;
  }

  return result;
}

void Config::addPiPulseSegment(std::vector<std::vector<PiPulseSegment>>& apply_pipulse, size_t oscilID, double tstart,
                               double tstop, double amp) const {
  if (oscilID < getNumOsc()) {
    PiPulseSegment parameterization = {tstart, tstop, amp};
    apply_pipulse[oscilID].push_back(parameterization);

    logger.log("Applying PiPulse to oscillator " + std::to_string(oscilID) + " in [" + std::to_string(tstart) + ", " +
               std::to_string(tstop) + "]: |p+iq|=" + std::to_string(amp) + "\n");

    // Set zero control for all other oscillators during this pipulse
    for (size_t i = 0; i < getNumOsc(); i++) {
      if (i != oscilID) {
        PiPulseSegment zero_parameterization = {tstart, tstop, 0.0};
        apply_pipulse[i].push_back(zero_parameterization);
      }
    }
  }
}

ControlParameterization Config::parseControlParameterization(const toml::table& table) const {
  ControlParameterization parameterization;

  std::string type_str = validators::field<std::string>(table, "type").value();
  std::optional<ControlType> type = parseEnum(type_str, CONTROL_TYPE_MAP);
  if (!type.has_value()) {
    logger.exitWithError("Unrecognized type '" + type_str + "' in control parameterization.");
  }
  parameterization.type = *type;

  const std::string type_key = "type";
  const std::string num_key = "num";
  const std::string scaling_key = "scaling";
  const std::string tstart_key = "tstart";
  const std::string tstop_key = "tstop";

  std::set<std::string> allowed_keys = {OSC_ID_KEY, type_key, num_key, scaling_key, tstart_key,tstop_key};
  validateTableKeys(table, allowed_keys, "control_parameterizations");

  switch (*type) {
    case ControlType::BSPLINE:
    case ControlType::BSPLINE0: {
      parameterization.nspline = validators::field<size_t>(table, num_key).value();
      parameterization.tstart = validators::field<double>(table, tstart_key).valueOr(ConfigDefaults::CONTROL_TSTART);
      parameterization.tstop = validators::field<double>(table, tstop_key).valueOr(getTotalTime());
      break;
    }
    case ControlType::BSPLINEAMP: {
      parameterization.nspline = validators::field<size_t>(table, num_key).value();
      parameterization.scaling = validators::field<double>(table, scaling_key).value();
      parameterization.tstart = validators::field<double>(table, tstart_key).valueOr(ConfigDefaults::CONTROL_TSTART);
      parameterization.tstop = validators::field<double>(table, tstop_key).valueOr(getTotalTime());
      break;
    }
    case ControlType::NONE:
      // logger.exitWithError("Unexpected control type " + type_str);
      // Do nothing, no parameters needed
      break;
  }

  return parameterization;
}

std::vector<ControlParameterization> Config::parseControlParameterizations(const toml::array& array_of_tables,size_t num_entries) const {
  ControlParameterization default_parameterization;
  default_parameterization.type = ConfigDefaults::CONTROL_TYPE;
  default_parameterization.nspline = ConfigDefaults::CONTROL_SPLINE_COUNT;
  default_parameterization.tstart = ConfigDefaults::CONTROL_TSTART;
  default_parameterization.tstop = getTotalTime();

  std::vector<ControlParameterization> result(num_entries, default_parameterization);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    ControlParameterization parameterization = parseControlParameterization(table);

    if (auto osc_id_node = table[OSC_ID_KEY]) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, OSC_ID_KEY).lessThan(num_entries).value();
      result[osc_id] = parameterization;
    } else {
      // Apply to ALL oscillators
      for (size_t i = 0; i < num_entries; i++) {
        result[i] = parameterization;
      }
    }
  }

  return result;
}

std::vector<ControlInitialization> Config::parseControlInitializations(const toml::table& table, size_t num_entries) const {

  const std::string type_key = "type";
  const std::string filename_key = "filename";
  const std::string amplitude_key = "amplitude";
  const std::string phase_key = "phase";
  const std::set<std::string> allowed_keys = {type_key, filename_key, amplitude_key, phase_key};

  // Helper function to parse a single initialization spec from a table
  auto parseInitSpec = [&](const toml::table& init_table, bool validate) -> ControlInitialization {
    if (validate) {
      validateTableKeys(init_table, allowed_keys, "control_initialization");
    }
    
    ControlInitialization init;

    // Parse type of initialization
    std::string type = validators::field<std::string>(init_table, type_key).value();
    auto type_enum = parseEnum(type, CONTROL_INITIALIZATION_TYPE_MAP);
    if (!type_enum.has_value()) {
      logger.exitWithError("Unknown control initialization type: " + type);
    }
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
  };

  ControlInitialization default_init = ControlInitialization{ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt}; 
  std::vector<ControlInitialization> result(num_entries, default_init);

  // First pass: Check if there is a default initialization (table with "type" key directly)
  if (table.contains(type_key)) {
    // This is a single default initialization that applies to all subsystems. Don't validate the main table since it may contain subsystem ID keys
    ControlInitialization global_default = parseInitSpec(table, false);
    for (size_t i = 0; i < num_entries; i++) {
      result[i] = global_default;
    }
  }

  // Second pass: Parse subsystem-specific overrides (with keys "0", "1", etc.)
  for (auto& [key, value] : table) {
    std::string key_str(key.str());

    // Skip keys that are part of the default specs
    if (allowed_keys.find(key_str) != allowed_keys.end()) {
      continue;
    }
    
    // Try to parse the key as a subsystem ID
    try {
      size_t subsys_id = std::stoul(key_str);
      if (subsys_id >= num_entries) {
        logger.exitWithError("control_initialization: subsystem ID " + key_str + " out of range (must be < " + std::to_string(num_entries) + ")");
      }

      // Clear any default for this subsystem and add the override setting
      if (value.is_table()) {
        auto* subsys_table = value.as_table();
        result[subsys_id] = parseInitSpec(*subsys_table, true);
      } else {
        logger.exitWithError("control_initialization: value for subsystem '" + key_str + "' must be a table");
      }
    } catch (const std::invalid_argument& e) {
      // Not a numeric key - might be a typo or invalid configuration
      logger.exitWithError("control_initialization: unexpected key '" + key_str + "'.");
    }
  }

  return result;
}

OptimTargetSettings Config::parseOptimTarget(TargetType type, const std::optional<GateType>& gate_type,
                                             const std::optional<std::string>& gate_file,
                                             const std::optional<std::vector<size_t>>& levels,
                                             const std::optional<std::string>& file) const {
  OptimTargetSettings target_settings;
  target_settings.type = type;

  switch (type) {
    case TargetType::GATE: {
      target_settings.gate_type = gate_type.value_or(ConfigDefaults::GATE_TYPE);
      target_settings.gate_file = gate_file.value_or("");
      break;
    }

    case TargetType::PURE: {
      if (levels.has_value() && !levels->empty()) {
        std::vector<size_t> pure_levels = levels.value();
        pure_levels.resize(nlevels.size(), nlevels.back());

        for (size_t i = 0; i < nlevels.size(); i++) {
          if (pure_levels[i] >= nlevels[i]) {
            logger.exitWithError(
                "ERROR in config setting. The requested pure state target |" + std::to_string(pure_levels[i]) +
                "> exceeds the number of modeled levels for that oscillator (" + std::to_string(nlevels[i]) + ").\n");
          }
        }
        target_settings.levels = pure_levels;
      } else {
        logger.log(
            "# Warning: You want to prepare a pure state, but didn't specify which one."
            " Taking default: ground-state |0...0> \n");
        target_settings.levels = std::vector<size_t>(nlevels.size(), 0);
      }
      break;
    }

    case TargetType::FROMFILE: {
      if (!file.has_value()) {
        logger.exitWithError("Optimization target of type FILE must have a filename");
      }
      target_settings.file = file.value();
      break;
    }
  }

  return target_settings;
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

std::vector<ControlParameterization> Config::parseControlParameterizationsCfg(const std::optional<std::map<int, ControlParameterizationData>>& parameterizations_map) const {
  // Set up a default parameterization for one oscillator
  ControlParameterization default_parameterization;
  default_parameterization.type = ConfigDefaults::CONTROL_TYPE;
  default_parameterization.nspline = ConfigDefaults::CONTROL_SPLINE_COUNT;
  default_parameterization.tstart = ConfigDefaults::CONTROL_TSTART;
  default_parameterization.tstop = getTotalTime();

  // Populate default if paramterization is not specified
  if (!parameterizations_map.has_value()) {
    return std::vector<ControlParameterization>(nlevels.size(), default_parameterization);
  }

  // Otherwise, parse specified parameterizations for each oscillator
  auto parsed_parameterizations = std::vector<ControlParameterization>(nlevels.size(), default_parameterization);
  for (size_t i = 0; i < parsed_parameterizations.size(); i++) {
    if (parameterizations_map.value().find(static_cast<int>(i)) != parameterizations_map.value().end()) {

      // auto parameterization = parseControlParameterizationCfg(parameterizations_map.value().at(i));
      auto oscil_config = parameterizations_map.value().at(static_cast<int>(i));
      const auto& params = oscil_config.parameters;

      // Create and store the parameterization
      ControlParameterization parameterization;
      parameterization.type = oscil_config.control_type;
      if (oscil_config.control_type == ControlType::BSPLINE || oscil_config.control_type == ControlType::BSPLINE0) {
        assert(params.size() >= 1); // nspline is required, should be validated in CfgParser
        parameterization.nspline = static_cast<size_t>(params[0]);
        parameterization.tstart = params.size() > 1 ? params[1] : ConfigDefaults::CONTROL_TSTART;
        parameterization.tstop = params.size() > 2 ? params[2] : getTotalTime();
      } else if (oscil_config.control_type == ControlType::BSPLINEAMP) {
        assert(params.size() >= 2); // nspline and scaling are required, should be validated in CfgParser
        parameterization.nspline = static_cast<size_t>(params[0]);
        parameterization.scaling = static_cast<double>(params[1]);
        parameterization.tstart = params.size() > 2 ? params[2] : ConfigDefaults::CONTROL_TSTART;
        parameterization.tstop = params.size() > 3 ? params[3] : getTotalTime();
      } 
      parsed_parameterizations[i] = parameterization;
    }
  }
  return parsed_parameterizations;
}


std::vector<ControlInitialization> Config::parseControlInitializationsCfg(const std::optional<std::map<int, ControlInitialization>>& init_configs) const {

  ControlInitialization default_init = ControlInitialization{ConfigDefaults::CONTROL_INIT_TYPE, ConfigDefaults::CONTROL_INIT_AMPLITUDE, std::nullopt, std::nullopt};

  std::vector<ControlInitialization> control_initializations(nlevels.size(), default_init);

  if (init_configs.has_value()) {
    for (size_t i = 0; i < nlevels.size(); i++) {
      if (init_configs->find(static_cast<int>(i)) != init_configs->end()) {
        control_initializations[i] = init_configs->at(static_cast<int>(i));
      }
    }
  }
  
  return control_initializations;
}
