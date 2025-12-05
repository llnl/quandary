#include "config.hpp"

#include <cassert>
#include <cstddef>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <toml++/impl/forward_declarations.hpp>
#include <toml++/impl/table.hpp>
#include <vector>

#include "cfgparser.hpp"
#include "config_types.hpp"
#include "config_validators.hpp"
#include "defs.hpp"
#include "util.hpp"

namespace {

// Helper function to convert enum back to string using existing enum maps
template <typename EnumType>
std::string enumToString(EnumType value, const std::map<std::string, EnumType>& type_map) {
  for (const auto& [str, enum_val] : type_map) {
    if (enum_val == value) return str;
  }
  return "unknown";
}

} // namespace

Config::Config(const MPILogger& logger, const toml::table& table) : logger(logger) {
  try {
    nlevels = validators::vectorField<size_t>(table, "nlevels").minLength(1).positive().value();

    size_t num_osc = nlevels.size();
    size_t num_pairs_osc = (num_osc - 1) * num_osc / 2;

    nessential = validators::vectorField<size_t>(table, "nessential").minLength(1).positive().valueOr(nlevels);
    copyLast(nessential, num_osc);

    ntime = validators::field<size_t>(table, "ntime").positive().valueOr(ntime);

    dt = validators::field<double>(table, "dt").positive().valueOr(dt);

    transfreq = validators::vectorField<double>(table, "transfreq").minLength(1).value();
    copyLast(transfreq, num_osc);

    selfkerr =
        validators::vectorField<double>(table, "selfkerr").minLength(1).valueOr(std::vector<double>(num_osc, 0.0));
    copyLast(selfkerr, num_osc);

    crosskerr = validators::vectorField<double>(table, "crosskerr")
                    .minLength(1)
                    .valueOr(std::vector<double>(num_pairs_osc, 0.0));
    copyLast(crosskerr, num_pairs_osc);

    Jkl = validators::vectorField<double>(table, "Jkl").minLength(1).valueOr(std::vector<double>(num_pairs_osc, 0.0));
    copyLast(Jkl, num_pairs_osc);

    rotfreq = validators::vectorField<double>(table, "rotfreq").minLength(1).value();
    copyLast(rotfreq, num_osc);

    std::string collapse_type_str = validators::field<std::string>(table, "collapse_type").valueOr("none");
    collapse_type = parseEnum(collapse_type_str, LINDBLAD_TYPE_MAP).value_or(LindbladType::NONE);

    decay_time = validators::vectorField<double>(table, "decay_time").valueOr(std::vector<double>(num_osc, 0.0));
    copyLast(decay_time, num_osc);

    dephase_time = validators::vectorField<double>(table, "dephase_time").valueOr(std::vector<double>(num_osc, 0.0));
    copyLast(dephase_time, num_osc);

    auto init_cond_table = validators::getRequiredTable(table, "initial_condition");
    std::string type_str = validators::field<std::string>(init_cond_table, "type").value();
    std::optional<std::vector<size_t>> levels = validators::getOptionalVector<size_t>(init_cond_table["levels"]);
    std::optional<std::vector<size_t>> osc_IDs = validators::getOptionalVector<size_t>(init_cond_table["oscIDs"]);
    std::optional<std::string> filename = init_cond_table["filename"].value<std::string>();
    initial_condition = parseInitialCondition({type_str, osc_IDs, levels, filename});
    n_initial_conditions = computeNumInitialConditions();

    apply_pipulse = std::vector<std::vector<PiPulseSegment>>(nlevels.size());
    auto apply_pipulse_array_of_tables = validators::getOptionalArrayOfTables(table, "apply_pipulse");
    for (auto& elem : apply_pipulse_array_of_tables) {
      auto pipulse_table = *elem.as_table();
      size_t oscilID = validators::field<size_t>(pipulse_table, "oscID").value();
      double tstart = validators::field<double>(pipulse_table, "tstart").value();
      double tstop = validators::field<double>(pipulse_table, "tstop").value();
      double amp = validators::field<double>(pipulse_table, "amp").value();

      addPiPulseSegment(apply_pipulse, oscilID, tstart, tstop, amp);
    }

    hamiltonian_file_Hsys = table["hamiltonian_file_Hsys"].value<std::string>();
    hamiltonian_file_Hc = table["hamiltonian_file_Hc"].value<std::string>();

    control_segments.resize(num_osc);
    control_initializations.resize(num_osc);
    control_bounds.resize(num_osc);
    carrier_frequencies.resize(num_osc);

    std::map<size_t, std::vector<ControlSegment>> control_segments_parsed;
    auto control_seg_node = table["control_segments"];
    if (control_seg_node.is_array_of_tables()) {
      for (auto& elem : *control_seg_node.as_array()) {
        auto seg_table = *elem.as_table();
        size_t oscilID = validators::field<size_t>(seg_table, "oscID").value();
        ControlSegment control_seg = parseControlSegment(seg_table);
        control_segments_parsed[oscilID].push_back(control_seg);
      }
    }

    // Parse control initialization
    std::map<size_t, std::vector<ControlSegmentInitialization>> osc_inits;

    if (table.contains("control_initialization")) {
      auto init_node = table["control_initialization"];
      if (init_node.is_array_of_tables()) {
        for (auto& elem : *init_node.as_array()) {
          auto init_table = *elem.as_table();
          std::string type = validators::field<std::string>(init_table, "type").value();

          auto type_enum = parseEnum(type, CONTROL_SEGMENT_INIT_TYPE_MAP);
          if (!type_enum.has_value()) {
            logger.exitWithError("Unknown control initialization type: " + type);
          }

          switch (type_enum.value()) {
            case ControlSegmentInitType::FILE: {
              std::string filename = validators::field<std::string>(init_table, "filename").value();
              control_initialization_file = filename;
              break;
            }
            case ControlSegmentInitType::CONSTANT: {
              size_t oscID = validators::field<size_t>(init_table, "oscID").value();
              double amplitude = validators::field<double>(init_table, "amplitude").value();
              double phase = validators::field<double>(init_table, "phase").valueOr(0.0);
              ControlSegmentInitialization init = {ControlSegmentInitType::CONSTANT, amplitude, phase};
              osc_inits[oscID].push_back(init);
              break;
            }
            case ControlSegmentInitType::RANDOM: {
              size_t oscID = validators::field<size_t>(init_table, "oscID").value();
              double amplitude = validators::field<double>(init_table, "amplitude").valueOr(0.1);
              double phase = validators::field<double>(init_table, "phase").valueOr(0.0);
              ControlSegmentInitialization init = {ControlSegmentInitType::RANDOM, amplitude, phase};
              osc_inits[oscID].push_back(init);
              break;
            }
          }
        }
      }
    }

    std::optional<std::map<int, std::vector<double>>> control_bounds_opt = std::nullopt;
    auto control_bounds_node = table["control_bounds"];
    if (control_bounds_node.is_array_of_tables()) {
      control_bounds_opt = std::map<int, std::vector<double>>();
      for (auto& elem : *control_bounds_node.as_array()) {
        auto bounds_table = *elem.as_table();
        size_t oscilID = validators::field<size_t>(bounds_table, "oscID").value();
        (*control_bounds_opt)[oscilID] = validators::vectorField<double>(bounds_table, "values").value();
      }
    }
    std::optional<std::map<int, std::vector<double>>> carrier_freq_opt = std::nullopt;
    auto carrier_freq_node = table["carrier_frequency"];
    if (carrier_freq_node.is_array_of_tables()) {
      carrier_freq_opt = std::map<int, std::vector<double>>();
      for (auto& elem : *carrier_freq_node.as_array()) {
        auto freq_table = *elem.as_table();
        size_t oscilID = validators::field<size_t>(freq_table, "oscID").value();
        (*carrier_freq_opt)[oscilID] = validators::vectorField<double>(freq_table, "values").value();
      }
    }

    std::vector<ControlSegment> default_segments = {
        {ControlType::BSPLINE, SplineParams{DEFAULT_SPLINE_COUNT, 0.0, ntime * dt}}};
    std::vector<ControlSegmentInitialization> default_initialization = {ControlSegmentInitialization{
        ControlSegmentInitType::CONSTANT, DEFAULT_CONTROL_INIT_AMPLITUDE, DEFAULT_CONTROL_INIT_PHASE}};

    for (size_t i = 0; i < control_segments.size(); i++) {
      if (control_segments_parsed.find(i) != control_segments_parsed.end()) {
        control_segments[i] = control_segments_parsed[i];
        default_segments = control_segments_parsed[i];
      } else {
        control_segments[i] = default_segments;
      }
      if (osc_inits.find(i) != osc_inits.end()) {
        control_initializations[i] = osc_inits[i];
        default_initialization = osc_inits[i];
      } else {
        control_initializations[i] = default_initialization;
        size_t num_segments = control_segments[i].size();
        copyLast(control_initializations[i], num_segments);
      }
    }

    control_bounds =
        parseIndexedWithDefaults<double>(control_bounds_opt, control_segments.size(), {DEFAULT_CONTROL_BOUND});
    // Extend bounds to match number of control segments
    for (size_t i = 0; i < control_bounds.size(); i++) {
      copyLast(control_bounds[i], control_segments[i].size());
    }

    carrier_frequencies = parseIndexedWithDefaults<double>(carrier_freq_opt, num_osc, {DEFAULT_CARRIER_FREQ});

    control_enforceBC = validators::field<bool>(table, "control_enforceBC").valueOr(control_enforceBC);

    // optim_target
    std::optional<OptimTargetData> optim_target_config;
    if (table.contains("optim_target")) {
      auto target_table = *table["optim_target"].as_table();
      std::string type_str = validators::field<std::string>(target_table, "target_type").value();
      std::optional<std::string> gate_type_str = target_table["gate_type"].value<std::string>();
      std::optional<std::string> gate_file = target_table["gate_file"].value<std::string>();
      std::optional<std::vector<size_t>> levels = validators::getOptionalVector<size_t>(target_table["levels"]);
      std::optional<std::string> filename = target_table["filename"].value<std::string>();
      optim_target_config = {type_str, gate_type_str, filename, gate_file, levels};
    }
    optim_target = parseOptimTarget(optim_target_config, nlevels);

    gate_rot_freq =
        validators::vectorField<double>(table, "gate_rot_freq").valueOr(std::vector<double>(num_osc, 0.0));
    copyLast(gate_rot_freq, num_osc);

    std::string optim_objective_str = validators::field<std::string>(table, "optim_objective").valueOr("");
    optim_objective = parseEnum(optim_objective_str, OBJECTIVE_TYPE_MAP).value_or(optim_objective);

    std::optional<std::vector<double>> optim_weights_opt =
        validators::getOptionalVector<double>(table["optim_weights"]);
    optim_weights = parseOptimWeights(optim_weights_opt);

    tolerance.atol = validators::field<double>(table, "optim_atol").positive().valueOr(tolerance.atol);
    tolerance.rtol = validators::field<double>(table, "optim_rtol").positive().valueOr(tolerance.rtol);
    tolerance.ftol = validators::field<double>(table, "optim_ftol").positive().valueOr(tolerance.ftol);
    tolerance.inftol = validators::field<double>(table, "optim_inftol").positive().valueOr(tolerance.inftol);
    tolerance.maxiter = validators::field<size_t>(table, "optim_maxiter").positive().valueOr(tolerance.maxiter);
    optim_regul = validators::field<double>(table, "optim_regul").greaterThanEqual(0.0).valueOr(optim_regul);

    penalty.penalty =
        validators::field<double>(table, "optim_penalty").greaterThanEqual(0.0).valueOr(penalty.penalty);
    penalty.penalty_param = validators::field<double>(table, "optim_penalty_param")
                                .greaterThanEqual(0.0)
                                .valueOr(penalty.penalty_param);
    penalty.penalty_dpdm = validators::field<double>(table, "optim_penalty_dpdm")
                               .greaterThanEqual(0.0)
                               .valueOr(penalty.penalty_dpdm);
    penalty.penalty_energy = validators::field<double>(table, "optim_penalty_energy")
                                 .greaterThanEqual(0.0)
                                 .valueOr(penalty.penalty_energy);
    penalty.penalty_variation = validators::field<double>(table, "optim_penalty_variation")
                                    .greaterThanEqual(0.0)
                                    .valueOr(penalty.penalty_variation);

    if (!table.contains("optim_regul_tik0") && table.contains("optim_regul_interpolate")) {
      // Handle deprecated optim_regul_interpolate logic
      optim_regul_tik0 = validators::field<bool>(table, "optim_regul_interpolate").value();
      logger.log("# Warning: 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
    }
    optim_regul_tik0 = validators::field<bool>(table, "optim_regul_tik0").valueOr(optim_regul_tik0);

    datadir = validators::field<std::string>(table, "datadir").valueOr(datadir);

    output_to_write.resize(num_osc); // Empty vectors by default
    auto write_array = validators::getOptionalArrayOfTables(table, "write");
    for (auto& elem : write_array) {
      auto write_table = *elem.as_table();
      size_t oscilID = validators::field<size_t>(write_table, "oscID").lessThan(num_osc).value();

      std::vector<std::string> types_str = validators::vectorField<std::string>(write_table, "type").value();
      output_to_write[oscilID] = convertStringVectorToEnum(types_str, OUTPUT_TYPE_MAP);
    }

    output_frequency = validators::field<size_t>(table, "output_frequency").positive().valueOr(output_frequency);
    optim_monitor_frequency =
        validators::field<size_t>(table, "optim_monitor_frequency").positive().valueOr(optim_monitor_frequency);

    std::string runtype_str = validators::field<std::string>(table, "runtype").valueOr("");
    runtype = parseEnum(runtype_str, RUN_TYPE_MAP).value_or(runtype);

    usematfree = validators::field<bool>(table, "usematfree").valueOr(usematfree);

    std::string linearsolver_type_str = validators::field<std::string>(table, "linearsolver_type").valueOr("");
    linearsolver_type = parseEnum(linearsolver_type_str, LINEAR_SOLVER_TYPE_MAP).value_or(linearsolver_type);

    linearsolver_maxiter =
        validators::field<size_t>(table, "linearsolver_maxiter").positive().valueOr(linearsolver_maxiter);

    std::string timestepper_type_str = validators::field<std::string>(table, "timestepper").valueOr("");
    timestepper_type = parseEnum(timestepper_type_str, TIME_STEPPER_TYPE_MAP).value_or(TimeStepperType::IMR);

    int rand_seed_ = validators::field<int>(table, "rand_seed").valueOr(-1);
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

  if (settings.ntime.has_value()) ntime = settings.ntime.value();

  if (settings.dt.has_value()) dt = settings.dt.value();

  if (!settings.transfreq.has_value()) {
    logger.exitWithError("transfreq cannot be empty");
  }
  transfreq = settings.transfreq.value();
  copyLast(transfreq, num_osc);

  selfkerr = settings.selfkerr.value_or(std::vector<double>(num_osc, 0.0));
  copyLast(selfkerr, num_osc);

  crosskerr = settings.crosskerr.value_or(std::vector<double>(num_pairs_osc, 0.0));
  copyLast(crosskerr, num_pairs_osc);

  Jkl = settings.Jkl.value_or(std::vector<double>(num_pairs_osc, 0.0));
  copyLast(Jkl, num_pairs_osc);

  if (!settings.rotfreq.has_value()) {
    logger.exitWithError("rotfreq cannot be empty");
  }
  rotfreq = settings.rotfreq.value();
  copyLast(rotfreq, num_osc);

  if (settings.collapse_type.has_value()) collapse_type = settings.collapse_type.value();

  decay_time = settings.decay_time.value_or(std::vector<double>(num_osc, 0.0));
  copyLast(decay_time, num_osc);

  dephase_time = settings.dephase_time.value_or(std::vector<double>(num_osc, 0.0));
  copyLast(dephase_time, num_osc);

  if (!settings.initialcondition.has_value()) {
    logger.exitWithError("initialcondition cannot be empty");
  }
  initial_condition = parseInitialCondition(settings.initialcondition.value());
  n_initial_conditions = computeNumInitialConditions();

  apply_pipulse = parsePiPulsesFromCfg(settings.apply_pipulse);

  hamiltonian_file_Hsys = settings.hamiltonian_file_Hsys;
  hamiltonian_file_Hc = settings.hamiltonian_file_Hc;

  // Control and optimization parameters
  carrier_frequencies.resize(num_osc);

  control_segments = parseControlSegments(settings.indexed_control_segments);

  // Control initialization
  if (settings.indexed_control_init.has_value()) {
    auto init_map = settings.indexed_control_init.value();
    if (init_map.find(0) != init_map.end() && !init_map[0].empty() && init_map[0][0].filename.has_value()) {
      control_initialization_file = init_map[0][0].filename;
      control_initializations.resize(num_osc);
      // Populate with default initialization for each oscillator, extended to match segments
      ControlSegmentInitialization default_init = ControlSegmentInitialization{
          ControlSegmentInitType::CONSTANT, DEFAULT_CONTROL_INIT_AMPLITUDE, DEFAULT_CONTROL_INIT_PHASE};
      std::vector<ControlSegmentInitialization> default_initialization = {default_init};
      for (size_t i = 0; i < num_osc; i++) {
        control_initializations[i] = default_initialization;
        size_t num_segments = control_segments[i].size();
        copyLast(control_initializations[i], num_segments);
      }
    } else {
      control_initializations = parseControlInitializations(settings.indexed_control_init);
    }
  }

  if (settings.control_enforceBC.has_value()) control_enforceBC = settings.control_enforceBC.value();
  control_bounds = parseIndexedWithDefaults<double>(settings.indexed_control_bounds, control_segments.size(),
                                                    {DEFAULT_CONTROL_BOUND});
  // Extend bounds to match number of control segments
  for (size_t i = 0; i < control_bounds.size(); i++) {
    copyLast(control_bounds[i], control_segments[i].size());
  }

  carrier_frequencies =
      parseIndexedWithDefaults<double>(settings.indexed_carrier_frequencies, num_osc, {DEFAULT_CARRIER_FREQ});
  optim_target = parseOptimTarget(settings.optim_target, nlevels);

  if (settings.gate_rot_freq.has_value()) gate_rot_freq = settings.gate_rot_freq.value();
  copyLast(gate_rot_freq, num_osc);

  if (settings.optim_objective.has_value()) optim_objective = settings.optim_objective.value();

  optim_weights = parseOptimWeights(settings.optim_weights);

  tolerance = OptimTolerance{};
  if (settings.optim_atol.has_value()) tolerance.atol = settings.optim_atol.value();
  if (settings.optim_rtol.has_value()) tolerance.rtol = settings.optim_rtol.value();
  if (settings.optim_ftol.has_value()) tolerance.ftol = settings.optim_ftol.value();
  if (settings.optim_inftol.has_value()) tolerance.inftol = settings.optim_inftol.value();
  if (settings.optim_maxiter.has_value()) tolerance.maxiter = settings.optim_maxiter.value();

  if (settings.optim_regul.has_value()) optim_regul = settings.optim_regul.value();

  penalty = OptimPenalty{};
  if (settings.optim_penalty.has_value()) penalty.penalty = settings.optim_penalty.value();
  if (settings.optim_penalty_param.has_value()) penalty.penalty_param = settings.optim_penalty_param.value();
  if (settings.optim_penalty_dpdm.has_value()) penalty.penalty_dpdm = settings.optim_penalty_dpdm.value();
  if (settings.optim_penalty_energy.has_value()) penalty.penalty_energy = settings.optim_penalty_energy.value();
  if (settings.optim_penalty_variation.has_value())
    penalty.penalty_variation = settings.optim_penalty_variation.value();

  if (settings.optim_regul_tik0.has_value()) {
    optim_regul_tik0 = settings.optim_regul_tik0.value();
  } else if (settings.optim_regul_interpolate.has_value()) {
    // Handle deprecated optim_regul_interpolate logic
    optim_regul_tik0 = settings.optim_regul_interpolate.value();
    logger.log("# Warning: 'optim_regul_interpolate' is deprecated. Please use 'optim_regul_tik0' instead.\n");
  }

  // Output parameters
  if (settings.datadir.has_value()) datadir = settings.datadir.value();
  output_to_write = parseIndexedWithDefaults<OutputType>(settings.indexed_output, num_osc);
  if (settings.output_frequency.has_value()) output_frequency = settings.output_frequency.value();
  if (settings.optim_monitor_frequency.has_value()) optim_monitor_frequency = settings.optim_monitor_frequency.value();
  if (settings.runtype.has_value()) runtype = settings.runtype.value();
  if (settings.usematfree.has_value()) usematfree = settings.usematfree.value();
  if (settings.linearsolver_type.has_value()) linearsolver_type = settings.linearsolver_type.value();
  if (settings.linearsolver_maxiter.has_value()) linearsolver_maxiter = settings.linearsolver_maxiter.value();
  if (settings.timestepper_type.has_value()) timestepper_type = settings.timestepper_type.value();
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

std::string print(const InitialCondition& initial_condition) {
  return std::visit([](const auto& opt) { return opt.toString(); }, initial_condition);
}

} //namespace

std::string ControlSegmentInitialization::toString() const {
  std::string str = "type = \"";
  str += enumToString(type, CONTROL_SEGMENT_INIT_TYPE_MAP);
  str += "\"\n";
  str += "amplitude = " + std::to_string(amplitude) + "\n";
  str += "phase = " + std::to_string(phase);
  return str;
}

std::string FromFileInitialCondition::toString() const { return "{type = \"file\", filename = \"" + filename + "}"; }

std::string PureInitialCondition::toString() const {
  std::string out = "{type = \"pure\", levels = ";
  out += printVector(levels);
  out += "}";
  return out;
}

std::string OscillatorIDsInitialCondition::toString(std::string name) const {
  std::string out = "{type = \"" + name + "\", oscIDs = ";
  out += printVector(osc_IDs);
  out += "}";
  return out;
}

void Config::printConfig(std::stringstream& log) const {
  log << "# Configuration settings\n";
  log << "# =============================================\n\n";

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
  log << "initial_condition = " << print(initial_condition) << "\n";

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

  if (hamiltonian_file_Hsys.has_value()) {
    log << "hamiltonian_file_Hsys = \"" << hamiltonian_file_Hsys.value() << "\"\n";
  }
  if (hamiltonian_file_Hc.has_value()) {
    log << "hamiltonian_file_Hc = \"" << hamiltonian_file_Hc.value() << "\"\n";
  }

  // Optimization parameters
  log << "\n";

  log << "control_enforce_BC = " << (control_enforceBC ? "true" : "false") << "\n";

  // Control initialization file
  if (control_initialization_file.has_value()) {
    log << "[control_initialization_file]\n";
    log << "filename = \"" << control_initialization_file.value() << "\"\n\n";
  }

  log << "optim_target = " << toString(optim_target) << "\n";
  log << "gate_rot_freq = " << printVector(gate_rot_freq) << "\n";
  log << "optim_objective = \"" << enumToString(optim_objective, OBJECTIVE_TYPE_MAP) << "\"\n";
  log << "optim_weights = " << printVector(optim_weights) << "\n";
  log << "optim_atol = " << tolerance.atol << "\n";
  log << "optim_rtol = " << tolerance.rtol << "\n";
  log << "optim_ftol = " << tolerance.ftol << "\n";
  log << "optim_inftol = " << tolerance.inftol << "\n";
  log << "optim_maxiter = " << tolerance.maxiter << "\n";
  log << "optim_regul = " << optim_regul << "\n";
  log << "optim_penalty = " << penalty.penalty << "\n";
  log << "optim_penalty_param = " << penalty.penalty_param << "\n";
  log << "optim_penalty_dpdm = " << penalty.penalty_dpdm << "\n";
  log << "optim_penalty_energy = " << penalty.penalty_energy << "\n";
  log << "optim_penalty_variation = " << penalty.penalty_variation << "\n";
  log << "optim_regul_tik0 = " << (optim_regul_tik0 ? "true" : "false") << "\n";

  // Control segments as array of tables
  for (size_t i = 0; i < control_segments.size(); ++i) {
    if (!control_segments[i].empty()) {
      const auto& seg = control_segments[i][0];
      log << "[[control_segments]]\n";
      log << "oscID = " << i << "\n";
      log << "type = \"" << enumToString(seg.type, CONTROL_TYPE_MAP) << "\"\n";

      // Add segment-specific parameters
      if (std::holds_alternative<SplineParams>(seg.params)) {
        auto params = std::get<SplineParams>(seg.params);
        log << "num = " << params.nspline << "\n";
        if (params.tstart != 0.0) log << "tstart = " << params.tstart << "\n";
        if (params.tstop != dt * ntime) log << "tstop = " << params.tstop << "\n";
      } else if (std::holds_alternative<SplineAmpParams>(seg.params)) {
        auto params = std::get<SplineAmpParams>(seg.params);
        log << "num = " << params.nspline << "\n";
        log << "scaling = " << params.scaling << "\n";
        if (params.tstart != 0.0) log << "tstart = " << params.tstart << "\n";
        if (params.tstop != dt * ntime) log << "tstop = " << params.tstop << "\n";
      } else if (std::holds_alternative<StepParams>(seg.params)) {
        auto params = std::get<StepParams>(seg.params);
        log << "step_amp1 = " << params.step_amp1 << "\n";
        log << "step_amp2 = " << params.step_amp2 << "\n";
        log << "tramp = " << params.tramp << "\n";
        log << "tstart = " << params.tstart << "\n";
        log << "tstop = " << params.tstop << "\n";
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
        log << init.toString() << "\n\n";
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

  // Output parameters
  log << "\n";
  log << "datadir = \"" << datadir << "\"\n";
  log << "output_frequency = " << output_frequency << "\n";
  log << "optim_monitor_frequency = " << optim_monitor_frequency << "\n";
  log << "runtype = \"" << enumToString(runtype, RUN_TYPE_MAP) << "\"\n";
  log << "usematfree = " << (usematfree ? "true" : "false") << "\n";
  log << "linearsolver_type = \"" << enumToString(linearsolver_type, LINEAR_SOLVER_TYPE_MAP) << "\"\n";
  log << "linearsolver_maxiter = " << linearsolver_maxiter << "\n";
  log << "timestepper = \"" << enumToString(timestepper_type, TIME_STEPPER_TYPE_MAP) << "\"\n";
  log << "rand_seed = " << rand_seed << "\n";

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

  log << "# =============================================\n\n";
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
}

size_t Config::computeNumInitialConditions() const {
  size_t n_initial_conditions = 0;
  if (std::holds_alternative<FromFileInitialCondition>(initial_condition))
    n_initial_conditions = 1;
  else if (std::holds_alternative<PureInitialCondition>(initial_condition))
    n_initial_conditions = 1;
  else if (std::holds_alternative<PerformanceInitialCondition>(initial_condition))
    n_initial_conditions = 1;
  else if (std::holds_alternative<EnsembleInitialCondition>(initial_condition))
    n_initial_conditions = 1;
  else if (std::holds_alternative<ThreeStatesInitialCondition>(initial_condition))
    n_initial_conditions = 3;
  else if (std::holds_alternative<NPlusOneInitialCondition>(initial_condition)) {
    // compute system dimension N
    n_initial_conditions = 1;
    for (size_t i = 0; i < nlevels.size(); i++) {
      n_initial_conditions *= nlevels[i];
    }
    n_initial_conditions += 1;
  } else if (std::holds_alternative<DiagonalInitialCondition>(initial_condition)) {
    /* Compute ninit = dim(subsystem defined by list of oscil IDs) */
    const auto& diag_init = std::get<DiagonalInitialCondition>(initial_condition);
    const auto& osc_IDs = diag_init.osc_IDs;

    n_initial_conditions = 1;
    for (size_t oscilID : osc_IDs) {
      if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
    }
  } else if (std::holds_alternative<BasisInitialCondition>(initial_condition)) {
    /* Compute ninit = dim(subsystem defined by list of oscil IDs) */
    const auto& basis_init = std::get<BasisInitialCondition>(initial_condition);
    const auto& osc_IDs = basis_init.osc_IDs;

    n_initial_conditions = 1;
    for (size_t oscilID : osc_IDs) {
      if (oscilID < nessential.size()) n_initial_conditions *= nessential[oscilID];
    }
    // if Schroedinger solver: ninit = N, do nothing.
    // else Lindblad solver: ninit = N^2
    if (collapse_type != LindbladType::NONE) {
      n_initial_conditions = (size_t)pow(n_initial_conditions, 2.0);
    }
  }
  logger.log("Number of initial conditions: " + std::to_string(n_initial_conditions) + "\n");
  return n_initial_conditions;
}

void Config::setRandSeed(std::optional<int> rand_seed_) {
  rand_seed = rand_seed_.value_or(-1);
  if (rand_seed < 0) {
    std::random_device rd;
    rand_seed = rd(); // random non-reproducable seed
  }
}

template <typename T>
std::vector<std::vector<T>> Config::parseIndexedWithDefaults(
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

template <typename EnumType>
std::vector<EnumType> Config::convertStringVectorToEnum(const std::vector<std::string>& strings,
                                                        const std::map<std::string, EnumType>& type_map) const {
  std::vector<EnumType> result;
  result.reserve(strings.size());
  for (const auto& str : strings) {
    auto enum_val = parseEnum(str, type_map);
    if (!enum_val) {
      logger.exitWithError("Unknown enum value: " + str);
    }
    result.push_back(*enum_val);
  }
  return result;
}

InitialCondition Config::parseInitialCondition(const InitialConditionData& config) const {
  auto opt_type = parseEnum(config.type, INITCOND_TYPE_MAP);

  if (!opt_type.has_value()) {
    logger.exitWithError("initial condition type not found.");
  }
  InitialConditionType type = opt_type.value();

  /* Sanity check for Schrodinger solver initial conditions */
  if (collapse_type == LindbladType::NONE) {
    if (type == InitialConditionType::ENSEMBLE || type == InitialConditionType::THREESTATES ||
        type == InitialConditionType::NPLUSONE) {
      logger.exitWithError(
          "\n\n ERROR for initial condition setting: \n When running Schroedingers solver"
          " (collapse_type == NONE), the initial condition needs to be either 'pure' or 'from file' or 'diagonal' or "
          "'basis'."
          " Note that 'diagonal' and 'basis' in the Schroedinger case are the same (all unit vectors).\n\n");
    }
  }

  // If no params are given for BASIS, ENSEMBLE, or DIAGONAL, default to all oscillators
  auto init_cond_IDs = config.osc_IDs.value_or(std::vector<size_t>{});
  if (!config.osc_IDs.has_value() &&
      (type == InitialConditionType::BASIS || type == InitialConditionType::ENSEMBLE ||
       type == InitialConditionType::DIAGONAL)) {
    for (size_t i = 0; i < nlevels.size(); i++) {
      init_cond_IDs.push_back(i);
    }
  }

  switch (type) {
    case InitialConditionType::FROMFILE:
      if (!config.filename.has_value()) {
        logger.exitWithError("initialcondition of type FROMFILE must have a filename");
      }
      return FromFileInitialCondition{config.filename.value()};
    case InitialConditionType::PURE:
      if (!config.levels.has_value()) {
        logger.exitWithError("initialcondition of type PURE must have 'levels'");
      }
      if (config.levels.value().size() != nlevels.size()) {
        logger.exitWithError("initialcondition of type PURE must have exactly " + std::to_string(nlevels.size()) +
                             " parameters, got " + std::to_string(config.levels.value().size()));
      }
      for (size_t k = 0; k < config.levels.value().size(); k++) {
        if (config.levels.value()[k] >= nlevels[k]) {
          logger.exitWithError("ERROR in config setting. The requested pure state initialization " +
                               std::to_string(config.levels.value()[k]) +
                               " exceeds the number of allowed levels for that oscillator (" +
                               std::to_string(nlevels[k]) + ").\n");
        }
      }
      return PureInitialCondition{config.levels.value()};

    case InitialConditionType::BASIS:
      if (collapse_type == LindbladType::NONE) {
        // DIAGONAL and BASIS initial conditions in the Schroedinger case are the same. Overwrite it to DIAGONAL
        return DiagonalInitialCondition{init_cond_IDs};
      }
      return BasisInitialCondition{init_cond_IDs};

    case InitialConditionType::ENSEMBLE:
      if (init_cond_IDs.back() >= nlevels.size()) {
        logger.exitWithError("Last element in initialcondition params exceeds number of oscillators");
      }

      for (size_t i = 1; i < init_cond_IDs.size() - 1; i++) {
        if (init_cond_IDs[i] + 1 != init_cond_IDs[i + 1]) {
          logger.exitWithError("List of oscillators for ensemble initialization should be consecutive!\n");
        }
      }
      return EnsembleInitialCondition{init_cond_IDs};

    case InitialConditionType::DIAGONAL:
      return DiagonalInitialCondition{init_cond_IDs};

    case InitialConditionType::THREESTATES:
      return ThreeStatesInitialCondition{};

    case InitialConditionType::NPLUSONE:
      return NPlusOneInitialCondition{};

    case InitialConditionType::PERFORMANCE:
      return PerformanceInitialCondition{};
  }

  // Should not happen, but make compiler happy
  logger.exitWithError("Internal error: parseInitialCondition reached end of function.");
  return BasisInitialCondition({});
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

std::vector<std::vector<PiPulseSegment>> Config::parsePiPulsesFromCfg(
    const std::optional<std::vector<PiPulseData>>& pulses) const {
  auto apply_pipulse = std::vector<std::vector<PiPulseSegment>>(nlevels.size());

  if (!pulses.has_value()) {
    return apply_pipulse;
  }
  for (const auto& pulse_config : *pulses) {
    addPiPulseSegment(apply_pipulse, pulse_config.oscil_id, pulse_config.tstart, pulse_config.tstop, pulse_config.amp);
  }
  return apply_pipulse;
}

std::vector<std::vector<ControlSegment>> Config::parseControlSegments(
    const std::optional<std::map<int, std::vector<ControlSegmentData>>>& segments_opt) const {
  std::vector<ControlSegment> default_segments = {{ControlType::BSPLINE, SplineParams{10, 0.0, ntime * dt}}};

  if (!segments_opt.has_value()) {
    return std::vector<std::vector<ControlSegment>>(nlevels.size(), default_segments);
  }
  const auto segments = segments_opt.value();
  auto parsed_segments = std::vector<std::vector<ControlSegment>>(nlevels.size());
  for (size_t i = 0; i < parsed_segments.size(); i++) {
    if (segments.find(static_cast<int>(i)) != segments.end()) {
      auto parsed = parseOscControlSegments(segments.at(i));
      parsed_segments[i] = parsed;
      default_segments = parsed;
    } else {
      parsed_segments[i] = default_segments;
    }
  }
  return parsed_segments;
}

std::vector<ControlSegment> Config::parseOscControlSegments(const std::vector<ControlSegmentData>& segments) const {
  std::vector<ControlSegment> control_segs = std::vector<ControlSegment>();

  for (const auto& seg_config : segments) {
    control_segs.push_back(parseControlSegment(seg_config));
  }
  return control_segs;
}

ControlSegment Config::parseControlSegment(const ControlSegmentData& seg_config) const {
  const auto& params = seg_config.parameters;

  // Create appropriate params variant based on type
  ControlSegment segment;
  segment.type = seg_config.control_type;

  if (seg_config.control_type == ControlType::BSPLINE || seg_config.control_type == ControlType::BSPLINE0) {
    SplineParams spline_params;
    assert(params.size() >= 1); // nspline is required, should be validated in CfgParser
    spline_params.nspline = static_cast<size_t>(params[0]);
    spline_params.tstart = params.size() > 1 ? params[1] : 0.0;
    spline_params.tstop = params.size() > 2 ? params[2] : ntime * dt;
    segment.params = spline_params;
  } else if (seg_config.control_type == ControlType::BSPLINEAMP) {
    SplineAmpParams spline_amp_params;
    assert(params.size() >= 2); // nspline and scaling are required, should be validated in CfgParser
    spline_amp_params.nspline = static_cast<size_t>(params[0]);
    spline_amp_params.scaling = static_cast<double>(params[1]);
    spline_amp_params.tstart = params.size() > 2 ? params[2] : 0.0;
    spline_amp_params.tstop = params.size() > 3 ? params[3] : ntime * dt;
    segment.params = spline_amp_params;
  } else if (seg_config.control_type == ControlType::STEP) {
    StepParams step_params;
    assert(params.size() >= 3); // step_amp1, step_amp2, tramp are required, should be validated in CfgParser
    step_params.step_amp1 = static_cast<double>(params[0]);
    step_params.step_amp2 = static_cast<double>(params[1]);
    step_params.tramp = static_cast<double>(params[2]);
    step_params.tstart = params.size() > 3 ? params[3] : 0.0;
    step_params.tstop = params.size() > 4 ? params[4] : ntime * dt;
    segment.params = step_params;
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
      SplineParams spline_params;
      spline_params.nspline = validators::field<size_t>(table, "num").value();
      spline_params.tstart = validators::field<double>(table, "tstart").valueOr(0.0);
      spline_params.tstop = validators::field<double>(table, "tstop").valueOr(ntime * dt);
      segment.params = spline_params;
      break;
    }
    case ControlType::BSPLINEAMP: {
      SplineAmpParams spline_amp_params;
      spline_amp_params.nspline = validators::field<size_t>(table, "num").value();
      spline_amp_params.scaling = validators::field<double>(table, "scaling").value();
      spline_amp_params.tstart = validators::field<double>(table, "tstart").valueOr(0.0);
      spline_amp_params.tstop = validators::field<double>(table, "tstop").valueOr(ntime * dt);
      segment.params = spline_amp_params;
      break;
    }
    case ControlType::STEP:
      StepParams step_params;
      step_params.step_amp1 = validators::field<double>(table, "step_amp1").value();
      step_params.step_amp2 = validators::field<double>(table, "step_amp2").value();
      step_params.tramp = validators::field<double>(table, "tramp").value();
      step_params.tstart = validators::field<double>(table, "tstart").valueOr(0.0);
      step_params.tstop = validators::field<double>(table, "tstop").valueOr(ntime * dt);
      segment.params = step_params;
      break;
    case ControlType::NONE:
      logger.exitWithError("Unexpected control type " + type_str);
  }

  return segment;
}

std::vector<std::vector<ControlSegmentInitialization>> Config::parseControlInitializations(
    const std::optional<std::map<int, std::vector<ControlInitializationData>>>& init_configs) const {
  ControlSegmentInitialization default_init = ControlSegmentInitialization{ControlSegmentInitType::CONSTANT, 0.0, 0.0};

  std::vector<std::vector<ControlSegmentInitialization>> control_initializations(nlevels.size());
  for (size_t i = 0; i < nlevels.size(); i++) {
    if (!init_configs.has_value() || init_configs->find(static_cast<int>(i)) == init_configs->end()) {
      control_initializations[i] = {default_init};
      continue;
    }
    for (const auto& init_config : init_configs->at(static_cast<int>(i))) {
      ControlSegmentInitialization init = ControlSegmentInitialization{
          init_config.init_seg_type, init_config.amplitude.value(), init_config.phase.value_or(0.0)};

      default_init = init;
      control_initializations[i].push_back(init);
    }
  }
  return control_initializations;
}

ControlSegmentInitialization Config::parseControlInitialization(const toml::table& table) const {
  std::string type_str = validators::field<std::string>(table, "type").value();

  std::optional<ControlSegmentInitType> type = parseEnum(type_str, CONTROL_SEGMENT_INIT_TYPE_MAP);
  if (!type.has_value()) {
    logger.exitWithError("Unrecognized type '" + type_str + "' in control initialization.");
  }
  return ControlSegmentInitialization{type.value(), validators::field<double>(table, "amplitude").value(),
                                      validators::field<double>(table, "phase").valueOr(0.0)};
}

OptimTargetSettings Config::parseOptimTarget(const std::optional<OptimTargetData>& opt_config,
                                             const std::vector<size_t>& nlevels) const {
  if (!opt_config.has_value()) {
    return PureOptimTarget{};
  }

  const OptimTargetData& config = opt_config.value();

  // Convert target type string to enum
  auto type = parseEnum(config.target_type, TARGET_TYPE_MAP);
  if (!type.has_value()) {
    logger.exitWithError("Unknown optimization target type: " + config.target_type);
  }

  switch (*type) {
    case TargetType::GATE: {
      GateOptimTarget gate_target;
      gate_target.gate_type = config.gate_type.has_value()
          ? parseEnum(config.gate_type.value(), GATE_TYPE_MAP).value_or(GateType::NONE)
          : GateType::NONE;
      gate_target.gate_file = config.gate_file.value_or("");
      return gate_target;
    }

    case TargetType::PURE: {
      PureOptimTarget pure_target;

      if (!config.levels.has_value() || config.levels->empty()) {
        logger.log(
            "# Warning: You want to prepare a pure state, but didn't specify which one."
            " Taking default: ground-state |0...0> \n");
        pure_target.purestate_levels = std::vector<size_t>(nlevels.size(), 0);
        return pure_target;
      }

      // Copy levels and validate
      for (auto level : config.levels.value()) {
        pure_target.purestate_levels.push_back(static_cast<size_t>(level));
      }
      pure_target.purestate_levels.resize(nlevels.size(), nlevels.back());

      for (size_t i = 0; i < nlevels.size(); i++) {
        if (pure_target.purestate_levels[i] >= nlevels[i]) {
          logger.exitWithError("ERROR in config setting. The requested pure state target |" +
                               std::to_string(pure_target.purestate_levels[i]) +
                               "> exceeds the number of modeled levels for that oscillator (" +
                               std::to_string(nlevels[i]) + ").\n");
        }
      }

      return pure_target;
    }

    case TargetType::FROMFILE: {
      FileOptimTarget file_target;
      file_target.file = config.filename.value_or("");
      return file_target;
    }
  }

  // Should never reach here, but satisfy compiler
  return PureOptimTarget{};
}

std::vector<double> Config::parseOptimWeights(const std::optional<std::vector<double>>& optim_weights_) const {
  // Set optimization weights, default to uniform weights summing to one
  std::vector<double> optim_weights = optim_weights_.value_or(std::vector<double>{1.0});
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
