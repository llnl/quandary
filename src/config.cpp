#include "config.hpp"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <vector>

#include "config_defaults.hpp"
#include "config_validators.hpp"
#include "util.hpp"

Config::Config(const MPILogger& logger, const toml::table& toml) : logger(logger) {
  try {
    // General options
    nlevels = validators::vectorField<size_t>(toml, "nlevels").minLength(1).positive().value();

    size_t num_osc = nlevels.size();
    size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;

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

    crosskerr = validators::vectorField<double>(toml, "crosskerr")
                    .minLength(1)
                    .valueOr(std::vector<double>(num_pairs_osc, ConfigDefaults::CROSSKERR));
    copyLast(crosskerr, num_pairs_osc);

    Jkl = validators::vectorField<double>(toml, "Jkl")
              .minLength(1)
              .valueOr(std::vector<double>(num_pairs_osc, ConfigDefaults::JKL));
    copyLast(Jkl, num_pairs_osc);

    rotfreq = validators::vectorField<double>(toml, "rotfreq").minLength(1).value();
    copyLast(rotfreq, num_osc);

    collapse_type =
        parseEnum(toml["collapse_type"].value<std::string>(), LINDBLAD_TYPE_MAP, ConfigDefaults::COLLAPSE_TYPE_ENUM);

    decay_time = validators::vectorField<double>(toml, "decay_time")
                     .valueOr(std::vector<double>(num_osc, ConfigDefaults::DECAY_TIME));
    copyLast(decay_time, num_osc);

    dephase_time = validators::vectorField<double>(toml, "dephase_time")
                       .valueOr(std::vector<double>(num_osc, ConfigDefaults::DEPHASE_TIME));
    copyLast(dephase_time, num_osc);

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
      size_t oscilID = validators::field<size_t>(pipulse_table, "oscID").value();
      double tstart = validators::field<double>(pipulse_table, "tstart").value();
      double tstop = validators::field<double>(pipulse_table, "tstop").value();
      double amp = validators::field<double>(pipulse_table, "amp").value();

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
    control_bounds =
        parseOscillatorSettings<double>(control_bounds_array, num_osc, {ConfigDefaults::CONTROL_BOUND}, "values");
    // Extend bounds to match number of control segments
    for (size_t i = 0; i < control_bounds.size(); i++) {
      copyLast(control_bounds[i], control_segments[i].size());
    }

    // Parse carrier frequencies
    auto carrier_freq_array = validators::getArrayOfTables(toml, "carrier_frequency");
    carrier_frequencies =
        parseOscillatorSettings<double>(carrier_freq_array, num_osc, {ConfigDefaults::CARRIER_FREQ}, "values");

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
    optim_weights = parseOptimWeights(optim_weights_opt);

    optim_atol = validators::field<double>(toml, "optim_atol").positive().valueOr(ConfigDefaults::OPTIM_ATOL);
    optim_rtol = validators::field<double>(toml, "optim_rtol").positive().valueOr(ConfigDefaults::OPTIM_RTOL);
    optim_ftol = validators::field<double>(toml, "optim_ftol").positive().valueOr(ConfigDefaults::OPTIM_FTOL);
    optim_inftol = validators::field<double>(toml, "optim_inftol").positive().valueOr(ConfigDefaults::OPTIM_INFTOL);
    optim_maxiter = validators::field<size_t>(toml, "optim_maxiter").positive().valueOr(ConfigDefaults::OPTIM_MAXITER);
    optim_regul =
        validators::field<double>(toml, "optim_regul").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_REGUL);

    optim_penalty =
        validators::field<double>(toml, "optim_penalty").greaterThanEqual(0.0).valueOr(ConfigDefaults::OPTIM_PENALTY);
    optim_penalty_param = validators::field<double>(toml, "optim_penalty_param")
                              .greaterThanEqual(0.0)
                              .valueOr(ConfigDefaults::OPTIM_PENALTY_PARAM);
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
    auto write_str = parseOscillatorSettings<std::string>(write_array, num_osc, {}, "type");
    for (size_t i = 0; i < write_str.size(); i++) {
      for (const auto& str : write_str[i]) {
        auto enum_val = parseEnum(str, OUTPUT_TYPE_MAP);
        if (!enum_val) {
          logger.exitWithError("Unknown enum value: " + str);
        }
        output_to_write[i].push_back(*enum_val);
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

  // Finalize and validate
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

  optim_weights = parseOptimWeights(settings.optim_weights);

  optim_atol = settings.optim_atol.value_or(ConfigDefaults::OPTIM_ATOL);
  optim_rtol = settings.optim_rtol.value_or(ConfigDefaults::OPTIM_RTOL);
  optim_ftol = settings.optim_ftol.value_or(ConfigDefaults::OPTIM_FTOL);
  optim_inftol = settings.optim_inftol.value_or(ConfigDefaults::OPTIM_INFTOL);
  optim_maxiter = settings.optim_maxiter.value_or(ConfigDefaults::OPTIM_MAXITER);

  optim_regul = settings.optim_regul.value_or(ConfigDefaults::OPTIM_REGUL);

  optim_penalty = settings.optim_penalty.value_or(ConfigDefaults::OPTIM_PENALTY);
  optim_penalty_param = settings.optim_penalty_param.value_or(ConfigDefaults::OPTIM_PENALTY_PARAM);
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
  setRandSeed(settings.rand_seed);

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
  toml::table config = toml::parse_file(filename);
  return Config(logger, config);
}

Config Config::fromTomlString(const std::string& toml_content, const MPILogger& logger) {
  toml::table config = toml::parse(toml_content);
  return Config(logger, config);
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

template <typename T>
std::string printVector(std::vector<T> vec) {
  std::string out = "[";
  for (size_t i = 0; i < vec.size(); ++i) {
    out += std::to_string(vec[i]);
    if (i < vec.size() - 1) {
      out += ", ";
    }
  }
  out += "]";
  return out;
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
  log << "crosskerr = " << printVector(crosskerr) << "\n";
  log << "Jkl = " << printVector(Jkl) << "\n";
  log << "rotfreq = " << printVector(rotfreq) << "\n";
  log << "collapse_type = \"" << enumToString(collapse_type, LINDBLAD_TYPE_MAP) << "\"\n";
  log << "decay_time = " << printVector(decay_time) << "\n";
  log << "dephase_time = " << printVector(dephase_time) << "\n";
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
  log << "optim_atol = " << optim_atol << "\n";
  log << "optim_rtol = " << optim_rtol << "\n";
  log << "optim_ftol = " << optim_ftol << "\n";
  log << "optim_inftol = " << optim_inftol << "\n";
  log << "optim_maxiter = " << optim_maxiter << "\n";
  log << "optim_regul = " << optim_regul << "\n";
  log << "optim_penalty = " << optim_penalty << "\n";
  log << "optim_penalty_param = " << optim_penalty_param << "\n";
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

void Config::setRandSeed(std::optional<int> rand_seed_) {
  rand_seed = rand_seed_.value_or(ConfigDefaults::RAND_SEED);
  if (rand_seed < 0) {
    std::random_device rd;
    rand_seed = rd(); // random non-reproducable seed
  }
}

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

template <typename T>
std::vector<std::vector<T>> Config::parseOscillatorSettings(const toml::array& array_of_tables, size_t num_entries,
                                                            std::vector<T> default_values,
                                                            const std::string& field_name) const {
  std::vector<std::vector<T>> result(num_entries, default_values);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    std::vector<T> values = validators::vectorField<T>(table, field_name).value();

    if (auto osc_id_node = table["oscID"]) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, "oscID").lessThan(num_entries).value();
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

ControlSegment Config::parseControlSegment(const toml::table& table) const {
  ControlSegment segment;

  std::string type_str = validators::field<std::string>(table, "type").value();
  std::optional<ControlType> type = parseEnum(type_str, CONTROL_TYPE_MAP);
  if (!type.has_value()) {
    logger.exitWithError("Unrecognized type '" + type_str + "' in control segment.");
  }
  segment.type = *type;

  switch (*type) {
    case ControlType::BSPLINE:
    case ControlType::BSPLINE0: {
      segment.nspline = validators::field<size_t>(table, "num").value();
      segment.tstart = validators::field<double>(table, "tstart").valueOr(ConfigDefaults::CONTROL_SEG_TSTART);
      segment.tstop = validators::field<double>(table, "tstop").valueOr(getTotalTime());
      break;
    }
    case ControlType::BSPLINEAMP: {
      segment.nspline = validators::field<size_t>(table, "num").value();
      segment.scaling = validators::field<double>(table, "scaling").value();
      segment.tstart = validators::field<double>(table, "tstart").valueOr(ConfigDefaults::CONTROL_SEG_TSTART);
      segment.tstop = validators::field<double>(table, "tstop").valueOr(getTotalTime());
      break;
    }
    case ControlType::STEP:
      segment.step_amp1 = validators::field<double>(table, "step_amp1").value();
      segment.step_amp2 = validators::field<double>(table, "step_amp2").value();
      segment.tramp = validators::field<double>(table, "tramp").value();
      segment.tstart = validators::field<double>(table, "tstart").valueOr(ConfigDefaults::CONTROL_SEG_TSTART);
      segment.tstop = validators::field<double>(table, "tstop").valueOr(getTotalTime());
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
  std::vector<ControlSegment> default_segments = {default_segment};

  std::vector<std::vector<ControlSegment>> result(num_entries);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    ControlSegment segment = parseControlSegment(table);

    if (auto osc_id_node = table["oscID"]) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, "oscID").lessThan(num_entries).value();
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
      elem = default_segments;
    }
  }

  return result;
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

std::vector<std::vector<ControlSegmentInitialization>> Config::parseControlInitializations(
    const toml::array& array_of_tables, size_t num_entries, std::optional<std::string>& control_init_file) const {
  ControlSegmentInitialization default_init = ControlSegmentInitialization{
      ControlSegmentInitType::CONSTANT, ConfigDefaults::CONTROL_INIT_AMPLITUDE, ConfigDefaults::CONTROL_INIT_PHASE};
  std::vector<ControlSegmentInitialization> default_initialization = {default_init};

  std::vector<std::vector<ControlSegmentInitialization>> result(num_entries);

  for (auto& elem : array_of_tables) {
    auto table = *elem.as_table();
    std::string type = validators::field<std::string>(table, "type").value();

    auto type_enum = parseEnum(type, CONTROL_SEGMENT_INIT_TYPE_MAP);
    if (!type_enum.has_value()) {
      logger.exitWithError("Unknown control initialization type: " + type);
    }

    ControlSegmentInitialization init;
    init.type = type_enum.value();

    switch (type_enum.value()) {
      case ControlSegmentInitType::FILE: {
        std::string filename = validators::field<std::string>(table, "filename").value();
        control_init_file = filename;
        break;
      }
      case ControlSegmentInitType::CONSTANT: {
        init.amplitude = validators::field<double>(table, "amplitude").value();
        init.phase = validators::field<double>(table, "phase").valueOr(ConfigDefaults::CONTROL_INIT_PHASE);
        break;
      }
      case ControlSegmentInitType::RANDOM: {
        init.amplitude =
            validators::field<double>(table, "amplitude").valueOr(ConfigDefaults::CONTROL_INIT_RANDOM_AMPLITUDE);
        init.phase = validators::field<double>(table, "phase").valueOr(ConfigDefaults::CONTROL_INIT_PHASE);
        break;
      }
    }

    if (auto osc_id_node = table["oscID"]) {
      // Apply to specific oscillator
      size_t osc_id = validators::field<size_t>(table, "oscID").lessThan(num_entries).value();
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
      elem = default_initialization;
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

std::vector<double> Config::parseOptimWeights(const std::optional<std::vector<double>>& optim_weights_) const {
  // Set optimization weights, default to uniform weights summing to one
  std::vector<double> optim_weights = optim_weights_.value_or(std::vector<double>{ConfigDefaults::OPTIM_WEIGHT});
  copyLast(optim_weights, n_initial_conditions);

  if (optim_weights.size() != n_initial_conditions) {
    logger.log("Warning: optim_weights size must be less than or equal to number of initial conditions");
    optim_weights.resize(n_initial_conditions);
  }
  // Scale the weights such that they sum up to one: beta_i <- beta_i / (\sum_i beta_i)
  double scaleweights = 0.0;
  for (size_t i = 0; i < n_initial_conditions; i++) scaleweights += optim_weights[i];
  for (size_t i = 0; i < n_initial_conditions; i++) optim_weights[i] = optim_weights[i] / scaleweights;
  return optim_weights;
}
