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

    // Parse control segments
    auto control_seg_array = validators::getArrayOfTables(toml, "control_segments");
    control_segments = parseControlSegments(control_seg_array, num_osc);

    // Parse control initialization
    auto control_init_array = validators::getArrayOfTables(toml, "control_initialization");
    control_initializations = parseControlInitializations(control_init_array, num_osc, control_initialization_file);

    // Parse control bounds
    auto control_bounds_array = validators::getArrayOfTables(toml, "control_bounds");

    // Allowed keys for control_bounds table
    const std::string values_key = "values";
    const std::set<std::string> control_bounds_allowed_keys = {OSC_ID_KEY, values_key};

    control_bounds = parseOscillatorSettings<double>(control_bounds_array, num_osc, {ConfigDefaults::CONTROL_BOUND},
                                                     values_key, control_bounds_allowed_keys, "control_bounds");
    // Extend bounds to match number of control segments
    for (size_t i = 0; i < control_bounds.size(); i++) {
      copyLast(control_bounds[i], control_segments[i].size());
    }

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
    optim_regul =
        validators::field<double>(toml, "optim_regul").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_REGUL);

    optim_penalty_leakage = validators::field<double>(toml, "optim_penalty_leakage")
        .greaterThanEqual(0.0)
        .valueOr(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
    optim_penalty_weightedcost = validators::field<double>(toml, "optim_penalty_weightedcost")
        .greaterThanEqual(0.0)
        .valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
    optim_penalty_weightedcost_width = validators::field<double>(toml, "optim_penalty_weightedcost_width")
        .greaterThanEqual(0.0)
        .valueOr(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
    optim_penalty_dpdm = validators::field<double>(toml, "optim_penalty_dpdm")
                             .greaterThanEqual(0.0)
                             .valueOr(ConfigDefaults::OPTIM_PENALTY_DPDM);
    optim_penalty_energy = validators::field<double>(toml, "optim_penalty_energy")
                               .greaterThanEqual(0.0)
                               .valueOr(ConfigDefaults::OPTIM_PENALTY_ENERGY);
    optim_penalty_variation = validators::field<double>(toml, "optim_penalty_variation")
                                  .greaterThanEqual(0.0)
                                  .valueOr(ConfigDefaults::OPTIM_PENALTY_VARIATION);

    if (!toml.contains("optim_regul_tik0") && toml.contains("optim_regul_interpolate")) {
      // Handle deprecated optim_regul_interpolate logic
      optim_regul_tik0 = validators::field<bool>(toml, "optim_regul_interpolate").value();
      logger.log("# Warning: 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
    }
    optim_regul_tik0 = toml["optim_regul_tik0"].value_or(ConfigDefaults::OPTIM_REGUL_TIK0);

    datadir = toml["datadir"].value_or(ConfigDefaults::DATADIR);

    output_to_write.resize(num_osc); // Empty vectors by default
    auto write_array = validators::getArrayOfTables(toml, "write");

    // Allowed keys for write table
    const std::string write_type_key = "type";
    const std::set<std::string> write_allowed_keys = {OSC_ID_KEY, write_type_key};

    auto write_str =
        parseOscillatorSettings<std::string>(write_array, num_osc, {}, write_type_key, write_allowed_keys, "write");
    for (size_t i = 0; i < write_str.size(); i++) {
      for (const auto& str : write_str[i]) {
        auto enum_val = parseEnum(str, OUTPUT_TYPE_MAP);
        if (!enum_val.has_value()) {
          logger.exitWithError("Unknown output type: " + str);
        }
        output_to_write[i].push_back(enum_val.value());
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

  control_segments = parseControlSegmentsCfg(settings.indexed_control_segments);

  // Control initialization
  if (settings.indexed_control_init.has_value()) {
    auto init_map = settings.indexed_control_init.value();
    if (init_map.find(0) != init_map.end() && !init_map[0].empty() && init_map[0][0].filename.has_value()) {
      control_initialization_file = init_map[0][0].filename;
      control_initializations.resize(num_osc);
      // Populate with default initialization for each oscillator, extended to match segments
      ControlSegmentInitialization default_init = ControlSegmentInitialization{
          ControlSegmentInitType::CONSTANT, ConfigDefaults::CONTROL_INIT_AMPLITUDE, ConfigDefaults::CONTROL_INIT_PHASE};
      std::vector<ControlSegmentInitialization> default_initialization = {default_init};
      for (size_t i = 0; i < num_osc; i++) {
        control_initializations[i] = default_initialization;
        size_t num_segments = control_segments[i].size();
        copyLast(control_initializations[i], num_segments);
      }
    } else {
      control_initializations = parseControlInitializationsCfg(settings.indexed_control_init);
    }
  } else {
    // Initialize with defaults when no control initialization is provided
    control_initializations.resize(num_osc);
    ControlSegmentInitialization default_init = ControlSegmentInitialization{
        ControlSegmentInitType::CONSTANT, ConfigDefaults::CONTROL_INIT_AMPLITUDE, ConfigDefaults::CONTROL_INIT_PHASE};
    std::vector<ControlSegmentInitialization> default_initialization = {default_init};
    for (size_t i = 0; i < num_osc; i++) {
      control_initializations[i] = default_initialization;
      size_t num_segments = control_segments[i].size();
      copyLast(control_initializations[i], num_segments);
    }
  }

  control_bounds = parseOscillatorSettingsCfg<double>(settings.indexed_control_bounds, control_segments.size(),
                                                      {ConfigDefaults::CONTROL_BOUND});
  // Extend bounds to match number of control segments
  for (size_t i = 0; i < control_bounds.size(); i++) {
    copyLast(control_bounds[i], control_segments[i].size());
  }

  carrier_frequencies.resize(num_osc);
  carrier_frequencies =
      parseOscillatorSettingsCfg<double>(settings.indexed_carrier_frequencies, num_osc, {ConfigDefaults::CARRIER_FREQ});
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

  optim_regul = settings.optim_regul.value_or(ConfigDefaults::OPTIM_REGUL);

  optim_penalty_leakage = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY_LEAKAGE);
  optim_penalty_weightedcost = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST);
  optim_penalty_weightedcost_width = settings.optim_penalty_param.value_or(ConfigDefaults::OPTIM_PENALTY_WEIGHTEDCOST_WIDTH);
  optim_penalty_dpdm = settings.optim_penalty_dpdm.value_or(ConfigDefaults::OPTIM_PENALTY_DPDM);
  optim_penalty_energy = settings.optim_penalty_energy.value_or(ConfigDefaults::OPTIM_PENALTY_ENERGY);
  optim_penalty_variation = settings.optim_penalty_variation.value_or(ConfigDefaults::OPTIM_PENALTY_VARIATION);
  optim_regul_tik0 = settings.optim_regul_tik0.value_or(ConfigDefaults::OPTIM_REGUL_TIK0);
  if (settings.optim_regul_interpolate.has_value()) {
    // Handle deprecated optim_regul_interpolate logic
    optim_regul_tik0 = settings.optim_regul_interpolate.value();
    logger.log("# Warning: 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
  }

  // Output parameters
  datadir = settings.datadir.value_or(ConfigDefaults::DATADIR);
  output_to_write = parseOscillatorSettingsCfg<OutputType>(settings.indexed_output, num_osc);
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

std::string toString(const ControlSegmentInitialization& seg_init) {
  std::string str = "type = \"";
  str += enumToString(seg_init.type, CONTROL_SEGMENT_INIT_TYPE_MAP);
  str += "\"\n";
  str += "amplitude = " + std::to_string(seg_init.amplitude) + "\n";
  str += "phase = " + std::to_string(seg_init.phase);
  return str;
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
  log << "optim_regul = " << optim_regul << "\n";
  log << "optim_penalty_leakage = " << optim_penalty_leakage << "\n";
  log << "optim_penalty_weightedcost = " << optim_penalty_weightedcost << "\n";
  log << "optim_penalty_weightedcost_width = " << optim_penalty_weightedcost_width << "\n";
  log << "optim_penalty_dpdm = " << optim_penalty_dpdm << "\n";
  log << "optim_penalty_energy = " << optim_penalty_energy << "\n";
  log << "optim_penalty_variation = " << optim_penalty_variation << "\n";
  log << "optim_regul_tik0 = " << (optim_regul_tik0 ? "true" : "false") << "\n";

  // Output parameters
  log << "datadir = \"" << datadir << "\"\n";
  log << "output_frequency = " << output_frequency << "\n";
  log << "optim_monitor_frequency = " << optim_monitor_frequency << "\n";
  log << "runtype = \"" << enumToString(runtype, RUN_TYPE_MAP) << "\"\n";
  log << "usematfree = " << (usematfree ? "true" : "false") << "\n";
  log << "linearsolver_type = \"" << enumToString(linearsolver_type, LINEAR_SOLVER_TYPE_MAP) << "\"\n";
  log << "linearsolver_maxiter = " << linearsolver_maxiter << "\n";
  log << "timestepper = \"" << enumToString(timestepper_type, TIME_STEPPER_TYPE_MAP) << "\"\n";
  log << "rand_seed = " << rand_seed << "\n";

  // Section 2: All array-of-tables at the end
  log << "\n";

  // Apply pi-pulse array of tables
  for (size_t i = 0; i < apply_pipulse.size(); ++i) {
    for (const auto& segment : apply_pipulse[i]) {
      log << "[[apply_pipulse]]\n";
      log << "oscID = " << i << "\n";
      log << "tstart = " << segment.tstart << "\n";
      log << "tstop = " << segment.tstop << "\n";
      log << "amp = " << segment.amp << "\n";
      log << "\n";
    }
  }

  // Control segments as array of tables
  for (size_t i = 0; i < control_segments.size(); ++i) {
    if (!control_segments[i].empty()) {
      const auto& seg = control_segments[i][0];
      log << "[[control_segments]]\n";
      log << "oscID = " << i << "\n";
      log << "type = \"" << enumToString(seg.type, CONTROL_TYPE_MAP) << "\"\n";

      // Add segment-specific parameters
      if (seg.type == ControlType::BSPLINE || seg.type == ControlType::BSPLINE0) {
        log << "num = " << seg.nspline.value() << "\n";
        log << "tstart = " << seg.tstart.value() << "\n";
        log << "tstop = " << seg.tstop.value() << "\n";
      } else if (seg.type == ControlType::BSPLINEAMP) {
        log << "num = " << seg.nspline.value() << "\n";
        log << "scaling = " << seg.scaling.value() << "\n";
        log << "tstart = " << seg.tstart.value() << "\n";
        log << "tstop = " << seg.tstop.value() << "\n";
      } else if (seg.type == ControlType::STEP) {
        log << "step_amp1 = " << seg.step_amp1.value() << "\n";
        log << "step_amp2 = " << seg.step_amp2.value() << "\n";
        log << "tramp = " << seg.tramp.value() << "\n";
        log << "tstart = " << seg.tstart.value() << "\n";
        log << "tstop = " << seg.tstop.value() << "\n";
      }
      log << "\n";
    }
  }

  // Control initialization
  if (!control_initialization_file.has_value()) {
    for (size_t i = 0; i < control_initializations.size(); ++i) {
      if (!control_initializations[i].empty()) {
        const auto& init = control_initializations[i][0];
        log << "[[control_initialization]]\n";
        log << "oscID = " << i << "\n";
        log << toString(init) << "\n\n";
      }
    }
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

  // Output write specifications as array of tables
  for (size_t i = 0; i < output_to_write.size(); ++i) {
    if (!output_to_write[i].empty()) {
      log << "[[write]]\n";
      log << "oscID = " << i << "\n";
      log << "type = [";
      for (size_t j = 0; j < output_to_write[i].size(); ++j) {
        log << "\"" << enumToString(output_to_write[i][j], OUTPUT_TYPE_MAP) << "\"";
        if (j < output_to_write[i].size() - 1) log << ", ";
      }
      log << "]\n\n";
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
    PiPulseSegment segment = {tstart, tstop, amp};
    apply_pipulse[oscilID].push_back(segment);

    logger.log("Applying PiPulse to oscillator " + std::to_string(oscilID) + " in [" + std::to_string(tstart) + ", " +
               std::to_string(tstop) + "]: |p+iq|=" + std::to_string(amp) + "\n");

    // Set zero control for all other oscillators during this pipulse
    for (size_t i = 0; i < getNumOsc(); i++) {
      if (i != oscilID) {
        PiPulseSegment zero_segment = {tstart, tstop, 0.0};
        apply_pipulse[i].push_back(zero_segment);
      }
    }
  }
}

ControlSegment Config::parseControlSegment(const toml::table& table) const {
  ControlSegment segment;

  std::string type_str = validators::field<std::string>(table, "type").value();
  std::optional<ControlType> type = parseEnum(type_str, CONTROL_TYPE_MAP);
  if (!type.has_value()) {
    logger.exitWithError("Unrecognized type '" + type_str + "' in control segment.");
  }
  segment.type = *type;

  const std::string type_key = "type";
  const std::string num_key = "num";
  const std::string scaling_key = "scaling";
  const std::string tstart_key = "tstart";
  const std::string tstop_key = "tstop";
  const std::string step_amp1_key = "step_amp1";
  const std::string step_amp2_key = "step_amp2";
  const std::string tramp_key = "tramp";

  std::set<std::string> allowed_keys = {OSC_ID_KEY, type_key,      num_key,       scaling_key, tstart_key,
                                        tstop_key,  step_amp1_key, step_amp2_key, tramp_key};
  validateTableKeys(table, allowed_keys, "control_segments");

  switch (*type) {
    case ControlType::BSPLINE:
    case ControlType::BSPLINE0: {
      segment.nspline = validators::field<size_t>(table, num_key).value();
      segment.tstart = validators::field<double>(table, tstart_key).valueOr(ConfigDefaults::CONTROL_SEG_TSTART);
      segment.tstop = validators::field<double>(table, tstop_key).valueOr(getTotalTime());
      break;
    }
    case ControlType::BSPLINEAMP: {
      segment.nspline = validators::field<size_t>(table, num_key).value();
      segment.scaling = validators::field<double>(table, scaling_key).value();
      segment.tstart = validators::field<double>(table, tstart_key).valueOr(ConfigDefaults::CONTROL_SEG_TSTART);
      segment.tstop = validators::field<double>(table, tstop_key).valueOr(getTotalTime());
      break;
    }
    case ControlType::STEP:
      segment.step_amp1 = validators::field<double>(table, step_amp1_key).value();
      segment.step_amp2 = validators::field<double>(table, step_amp2_key).value();
      segment.tramp = validators::field<double>(table, tramp_key).value();
      segment.tstart = validators::field<double>(table, tstart_key).valueOr(ConfigDefaults::CONTROL_SEG_TSTART);
      segment.tstop = validators::field<double>(table, tstop_key).valueOr(getTotalTime());
      break;
    case ControlType::NONE:
      logger.exitWithError("Unexpected control type " + type_str);
  }

  return segment;
}

std::vector<std::vector<ControlSegment>> Config::parseControlSegments(const toml::array& array_of_tables,
                                                                      size_t num_entries) const {
  ControlSegment default_segment;
  default_segment.type = ControlType::BSPLINE;
  default_segment.nspline = ConfigDefaults::CONTROL_SEG_SPLINE_COUNT;
  default_segment.tstart = ConfigDefaults::CONTROL_SEG_TSTART;
  default_segment.tstop = getTotalTime();

  std::vector<std::vector<ControlSegment>> result(num_entries);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    ControlSegment segment = parseControlSegment(table);

    if (auto osc_id_node = table[OSC_ID_KEY]) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, OSC_ID_KEY).lessThan(num_entries).value();
      result[osc_id].push_back(segment);
    } else {
      // Apply to ALL oscillators
      for (size_t i = 0; i < num_entries; i++) {
        result[i].push_back(segment);
      }
    }
  }

  for (auto& elem : result) {
    if (elem.empty()) {
      elem.push_back(default_segment);
    }
  }

  return result;
}

std::vector<std::vector<ControlSegmentInitialization>> Config::parseControlInitializations(
    const toml::array& array_of_tables, size_t num_entries, std::optional<std::string>& control_init_file) const {
  ControlSegmentInitialization default_init = ControlSegmentInitialization{
      ControlSegmentInitType::CONSTANT, ConfigDefaults::CONTROL_INIT_AMPLITUDE, ConfigDefaults::CONTROL_INIT_PHASE};

  const std::string type_key = "type";
  const std::string filename_key = "filename";
  const std::string amplitude_key = "amplitude";
  const std::string phase_key = "phase";
  const std::set allowed_keys = {OSC_ID_KEY, type_key, filename_key, amplitude_key, phase_key};
  std::vector<std::vector<ControlSegmentInitialization>> result(num_entries);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    validateTableKeys(table, allowed_keys, "control_initialization");
    std::string type = validators::field<std::string>(table, type_key).value();

    auto type_enum = parseEnum(type, CONTROL_SEGMENT_INIT_TYPE_MAP);
    if (!type_enum.has_value()) {
      logger.exitWithError("Unknown control initialization type: " + type);
    }

    ControlSegmentInitialization init;
    init.type = type_enum.value();

    switch (type_enum.value()) {
      case ControlSegmentInitType::FILE: {
        std::string filename = validators::field<std::string>(table, filename_key).value();
        control_init_file = filename;
        break;
      }
      case ControlSegmentInitType::CONSTANT: {
        init.amplitude = validators::field<double>(table, amplitude_key).value();
        init.phase = validators::field<double>(table, phase_key).valueOr(ConfigDefaults::CONTROL_INIT_PHASE);
        break;
      }
      case ControlSegmentInitType::RANDOM: {
        init.amplitude =
            validators::field<double>(table, amplitude_key).valueOr(ConfigDefaults::CONTROL_INIT_RANDOM_AMPLITUDE);
        init.phase = validators::field<double>(table, phase_key).valueOr(ConfigDefaults::CONTROL_INIT_PHASE);
        break;
      }
    }

    if (table.contains(OSC_ID_KEY)) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, OSC_ID_KEY).lessThan(num_entries).value();
      result[osc_id].push_back(init);
    } else {
      // Apply to ALL oscillators
      for (size_t i = 0; i < num_entries; i++) {
        result[i].push_back(init);
      }
    }
  }

  for (auto& elem : result) {
    if (elem.empty()) {
      elem.push_back(default_init);
    }
  }

  // Extend initializations to match number of control segments
  for (size_t i = 0; i < control_initializations.size(); i++) {
    size_t num_segments = control_segments[i].size();
    copyLast(result[i], num_segments);
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
        logger.exitWithError("Optimization target of type FROMFILE must have a filename");
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

std::vector<std::vector<ControlSegment>> Config::parseControlSegmentsCfg(
    const std::optional<std::map<int, std::vector<ControlSegmentData>>>& segments_opt) const {
  ControlSegment default_segment;
  default_segment.type = ControlType::BSPLINE;
  default_segment.nspline = ConfigDefaults::CONTROL_SEG_SPLINE_COUNT;
  default_segment.tstart = ConfigDefaults::CONTROL_SEG_TSTART;
  default_segment.tstop = getTotalTime();
  std::vector<ControlSegment> default_segments = {{default_segment}};

  if (!segments_opt.has_value()) {
    return std::vector<std::vector<ControlSegment>>(nlevels.size(), default_segments);
  }
  const auto segments = segments_opt.value();
  auto parsed_segments = std::vector<std::vector<ControlSegment>>(nlevels.size());
  for (size_t i = 0; i < parsed_segments.size(); i++) {
    if (segments.find(static_cast<int>(i)) != segments.end()) {
      std::vector<ControlSegment> parsed;
      for (const auto& seg_config : segments.at(i)) {
        parsed.push_back(parseControlSegmentCfg(seg_config));
      }
      parsed_segments[i] = parsed;
      default_segments = parsed;
    } else {
      parsed_segments[i] = default_segments;
    }
  }
  return parsed_segments;
}

ControlSegment Config::parseControlSegmentCfg(const ControlSegmentData& seg_config) const {
  const auto& params = seg_config.parameters;

  ControlSegment segment;
  segment.type = seg_config.control_type;

  if (seg_config.control_type == ControlType::BSPLINE || seg_config.control_type == ControlType::BSPLINE0) {
    assert(params.size() >= 1); // nspline is required, should be validated in CfgParser
    segment.nspline = static_cast<size_t>(params[0]);
    segment.tstart = params.size() > 1 ? params[1] : ConfigDefaults::CONTROL_SEG_TSTART;
    segment.tstop = params.size() > 2 ? params[2] : getTotalTime();
  } else if (seg_config.control_type == ControlType::BSPLINEAMP) {
    assert(params.size() >= 2); // nspline and scaling are required, should be validated in CfgParser
    segment.nspline = static_cast<size_t>(params[0]);
    segment.scaling = static_cast<double>(params[1]);
    segment.tstart = params.size() > 2 ? params[2] : ConfigDefaults::CONTROL_SEG_TSTART;
    segment.tstop = params.size() > 3 ? params[3] : getTotalTime();
  } else if (seg_config.control_type == ControlType::STEP) {
    assert(params.size() >= 3); // step_amp1, step_amp2, tramp are required, should be validated in CfgParser
    segment.step_amp1 = static_cast<double>(params[0]);
    segment.step_amp2 = static_cast<double>(params[1]);
    segment.tramp = static_cast<double>(params[2]);
    segment.tstart = params.size() > 3 ? params[3] : ConfigDefaults::CONTROL_SEG_TSTART;
    segment.tstop = params.size() > 4 ? params[4] : getTotalTime();
  }

  return segment;
}

std::vector<std::vector<ControlSegmentInitialization>> Config::parseControlInitializationsCfg(
    const std::optional<std::map<int, std::vector<ControlInitializationData>>>& init_configs) const {
  ControlSegmentInitialization default_init = ControlSegmentInitialization{
      ControlSegmentInitType::CONSTANT, ConfigDefaults::CONTROL_INIT_AMPLITUDE, ConfigDefaults::CONTROL_INIT_PHASE};

  std::vector<std::vector<ControlSegmentInitialization>> control_initializations(nlevels.size());
  for (size_t i = 0; i < nlevels.size(); i++) {
    if (!init_configs.has_value() || init_configs->find(static_cast<int>(i)) == init_configs->end()) {
      control_initializations[i] = {default_init};
      continue;
    }
    for (const auto& init_config : init_configs->at(static_cast<int>(i))) {
      ControlSegmentInitialization init =
          ControlSegmentInitialization{init_config.init_seg_type, init_config.amplitude.value(),
                                       init_config.phase.value_or(ConfigDefaults::CONTROL_INIT_PHASE)};

      default_init = init;
      control_initializations[i].push_back(init);
    }
  }
  return control_initializations;
}
