#include "cfgparser.hpp"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <vector>

#include "defs.hpp"
#include "mpi_logger.hpp"
#include "util.hpp"

CfgParser::CfgParser(const MPILogger& logger) : logger(logger) {
  // Register config parameter setters
  // General options
  registerConfig("nlevels", settings.nlevels);
  registerConfig("nessential", settings.nessential);
  registerConfig("ntime", settings.ntime);
  registerConfig("dt", settings.dt);
  registerConfig("transfreq", settings.transfreq);
  registerConfig("selfkerr", settings.selfkerr);
  registerConfig("crosskerr", settings.crosskerr);
  registerConfig("Jkl", settings.Jkl);
  registerConfig("rotfreq", settings.rotfreq);
  registerConfig("collapse_type", settings.collapse_type);
  registerConfig("decay_time", settings.decay_time);
  registerConfig("dephase_time", settings.dephase_time);
  registerConfig("initialcondition", settings.initialcondition);
  registerConfig("apply_pipulse", settings.apply_pipulse);
  registerConfig("hamiltonian_file_Hsys", settings.hamiltonian_file_Hsys);
  registerConfig("hamiltonian_file_Hc", settings.hamiltonian_file_Hc);

  // Optimization options
  registerConfig("control_enforceBC", settings.control_enforceBC);
  registerConfig("optim_target", settings.optim_target);
  registerConfig("gate_rot_freq", settings.gate_rot_freq);
  registerConfig("optim_objective", settings.optim_objective);
  registerConfig("optim_weights", settings.optim_weights);

  // Indexed settings (per-oscillator)
  registerIndexedConfig("control_segments", settings.indexed_control_segments);
  registerIndexedConfig("control_initialization", settings.indexed_control_init);
  registerIndexedConfig("control_bounds", settings.indexed_control_bounds);
  registerIndexedConfig("carrier_frequency", settings.indexed_carrier_frequencies);
  registerIndexedConfig("output", settings.indexed_output);
  registerConfig("optim_atol", settings.optim_atol);
  registerConfig("optim_rtol", settings.optim_rtol);
  registerConfig("optim_ftol", settings.optim_ftol);
  registerConfig("optim_inftol", settings.optim_inftol);
  registerConfig("optim_maxiter", settings.optim_maxiter);
  registerConfig("optim_regul", settings.optim_regul);
  registerConfig("optim_penalty", settings.optim_penalty);
  registerConfig("optim_penalty_param", settings.optim_penalty_param);
  registerConfig("optim_penalty_dpdm", settings.optim_penalty_dpdm);
  registerConfig("optim_penalty_energy", settings.optim_penalty_energy);
  registerConfig("optim_penalty_variation", settings.optim_penalty_variation);
  registerConfig("optim_regul_tik0", settings.optim_regul_tik0);
  registerConfig("optim_regul_interpolate", optim_regul_interpolate);

  // Output and runtypes
  registerConfig("datadir", settings.datadir);
  // Note: "output" is handled as indexed config in the indexed settings section above
  registerConfig("output_frequency", settings.output_frequency);
  registerConfig("optim_monitor_frequency", settings.optim_monitor_frequency);
  registerConfig("runtype", settings.runtype);
  registerConfig("usematfree", settings.usematfree);
  registerConfig("linearsolver_type", settings.linearsolver_type);
  registerConfig("linearsolver_maxiter", settings.linearsolver_maxiter);
  registerConfig("timestepper", settings.timestepper_type);
  registerConfig("rand_seed", settings.rand_seed);
}

std::vector<std::vector<double>> CfgParser::convertIndexedToVectorVector(
    const std::map<int, std::vector<double>>& indexed_map, size_t num_oscillators) {
  std::vector<std::vector<double>> result(num_oscillators);
  for (const auto& [osc_idx, values] : indexed_map) {
    if (static_cast<size_t>(osc_idx) < result.size()) {
      result[osc_idx] = values;
    }
  }
  return result;
}

std::vector<std::vector<OutputType>> CfgParser::convertIndexedToOutputVector(
    const std::map<int, std::vector<OutputType>>& indexed_map, size_t num_oscillators) {
  std::vector<std::vector<OutputType>> result(num_oscillators);
  for (const auto& [osc_idx, values] : indexed_map) {
    if (static_cast<size_t>(osc_idx) < result.size()) {
      result[osc_idx] = values;
    }
  }
  return result;
}

namespace {

std::string trimWhitespace(std::string s) {
  s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); }), s.end());
  return s;
}

bool isComment(const std::string& line) { return line.size() > 0 && (line[0] == '#' || line[0] == '/'); }

bool isValidControlType(const std::string& str) {
  return CONTROL_TYPE_MAP.find(toLower(str)) != CONTROL_TYPE_MAP.end();
}

bool isValidControlSegmentInitType(const std::string& str) {
  return CONTROL_SEGMENT_INIT_TYPE_MAP.find(toLower(str)) != CONTROL_SEGMENT_INIT_TYPE_MAP.end();
}

} // namespace

std::vector<std::string> CfgParser::split(const std::string& str, char delimiter) {
  std::vector<std::string> result;
  std::stringstream ss(str);
  std::string item;
  while (std::getline(ss, item, delimiter)) {
    if (!item.empty()) {
      result.push_back(item);
    }
  }
  return result;
}

void CfgParser::applyConfigLine(const std::string& line) {
  std::string trimmedLine = trimWhitespace(line);
  if (!trimmedLine.empty() && !isComment(trimmedLine)) {
    int pos = trimmedLine.find('=');
    std::string key = trimmedLine.substr(0, pos);
    std::string value = trimmedLine.substr(pos + 1);

    // First try exact key match
    if (setters.count(key)) {
      try {
        setters[key](value);
      } catch (const std::exception& e) {
        logger.exitWithError("Error parsing '" + key + "': " + e.what());
      }
    } else {
      // Try to handle indexed settings (e.g., control_segments0, output1)
      bool handled = handleIndexedSetting(key, value);
      if (!handled) {
        logger.exitWithError("Unknown option '" + key + "'");
      }
    }
  }
}

bool CfgParser::handleIndexedSetting(const std::string& key, const std::string& value) {
  // Check if key ends with a digit (indexed setting)
  if (key.empty() || !std::isdigit(key.back())) {
    return false;
  }

  // Find where the index starts
  size_t index_pos = key.find_last_not_of("0123456789") + 1;
  if (index_pos == std::string::npos || index_pos >= key.length()) {
    return false;
  }

  std::string base_key = key.substr(0, index_pos);
  int index = std::stoi(key.substr(index_pos));

  if (indexed_setters.count(base_key)) {
    try {
      indexed_setters[base_key](index, value);
      return true;
    } catch (const std::exception& e) {
      logger.exitWithError("Error parsing indexed setting '" + key + "': " + e.what());
    }
  }

  return false;
}

ParsedConfigData CfgParser::parseFile(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    logger.exitWithError("Unable to read the file " + filename);
  }

  loadFromStream(file);
  file.close();
  return settings;
}

ParsedConfigData CfgParser::parseString(const std::string& config_content) {
  std::istringstream stream(config_content);
  loadFromStream(stream);
  return settings;
}

// Enum converter implementations
template <>
RunType CfgParser::convertFromString<RunType>(const std::string& str) {
  auto it = RUN_TYPE_MAP.find(toLower(str));
  if (it == RUN_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown run type: " + str + ".\n");
  }
  return it->second;
}

template <>
LindbladType CfgParser::convertFromString<LindbladType>(const std::string& str) {
  auto it = LINDBLAD_TYPE_MAP.find(toLower(str));
  if (it == LINDBLAD_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown Lindblad type: " + str + ".\n");
  }
  return it->second;
}

template <>
LinearSolverType CfgParser::convertFromString<LinearSolverType>(const std::string& str) {
  auto it = LINEAR_SOLVER_TYPE_MAP.find(toLower(str));
  if (it == LINEAR_SOLVER_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown linear solver type: " + str + ".\n");
  }
  return it->second;
}

template <>
TimeStepperType CfgParser::convertFromString<TimeStepperType>(const std::string& str) {
  auto it = TIME_STEPPER_TYPE_MAP.find(toLower(str));
  if (it == TIME_STEPPER_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown time stepper type: " + str + ".\n");
  }
  return it->second;
}

template <>
TargetType CfgParser::convertFromString<TargetType>(const std::string& str) {
  auto it = TARGET_TYPE_MAP.find(toLower(str));
  if (it == TARGET_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown target type: " + str + ".\n");
  }
  return it->second;
}

template <>
InitialConditionType CfgParser::convertFromString<InitialConditionType>(const std::string& str) {
  auto it = INITCOND_TYPE_MAP.find(toLower(str));
  if (it == INITCOND_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown initial condition type: " + str + ".\n");
  }
  return it->second;
}

template <>
GateType CfgParser::convertFromString<GateType>(const std::string& str) {
  auto it = GATE_TYPE_MAP.find(toLower(str));
  if (it == GATE_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown gate type: " + str + ".\n");
  }
  return it->second;
}

template <>
OutputType CfgParser::convertFromString<OutputType>(const std::string& str) {
  auto it = OUTPUT_TYPE_MAP.find(toLower(str));
  if (it == OUTPUT_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown output type: " + str + ".\n");
  }
  return it->second;
}

template <>
ObjectiveType CfgParser::convertFromString<ObjectiveType>(const std::string& str) {
  auto it = OBJECTIVE_TYPE_MAP.find(toLower(str));
  if (it == OBJECTIVE_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown objective type: " + str + ".\n");
  }
  return it->second;
}

template <>
ControlType CfgParser::convertFromString<ControlType>(const std::string& str) {
  auto it = CONTROL_TYPE_MAP.find(toLower(str));
  if (it == CONTROL_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown control type: " + str + ".\n");
  }
  return it->second;
}

template <>
ControlSegmentInitType CfgParser::convertFromString<ControlSegmentInitType>(const std::string& str) {
  auto it = CONTROL_SEGMENT_INIT_TYPE_MAP.find(toLower(str));
  if (it == CONTROL_SEGMENT_INIT_TYPE_MAP.end()) {
    logger.exitWithError("\n\n ERROR: Unknown control segment initialization type: " + str + ".\n");
  }
  return it->second;
}

// Struct converter implementations
template <>
InitialCondition CfgParser::convertFromString<InitialCondition>(const std::string& str) {
  auto parts = split(str);
  if (parts.empty()) {
    logger.exitWithError("Empty initialcondition specification");
  }

  InitialCondition init_cond;
  auto type = convertFromString<InitialConditionType>(parts[0]);
  init_cond.type = type;

  if (type == InitialConditionType::FROMFILE) {
    if (parts.size() < 2) {
      logger.exitWithError("initialcondition of type FROMFILE must have a filename");
    }
    init_cond.filename = parts[1];
  } else if (type == InitialConditionType::PURE) {
    init_cond.levels = std::vector<size_t>();
    for (size_t i = 1; i < parts.size(); ++i) {
      init_cond.levels.value().push_back(convertFromString<int>(parts[i]));
    }
  } else if (type == InitialConditionType::ENSEMBLE || type == InitialConditionType::DIAGONAL ||
             type == InitialConditionType::BASIS) {
    if (parts.size() > 1) {
      init_cond.osc_IDs = std::vector<size_t>();
      for (size_t i = 1; i < parts.size(); ++i) {
        init_cond.osc_IDs.value().push_back(convertFromString<int>(parts[i]));
      }
    }
  }

  return init_cond;
}

template <>
OptimTargetSettings CfgParser::convertFromString<OptimTargetSettings>(const std::string& str) {
  auto parts = split(str);
  if (parts.empty()) {
    logger.exitWithError("optim_target must have at least a target type specified.");
  }

  OptimTargetSettings target_settings;
  auto target_type = convertFromString<TargetType>(parts[0]);
  target_settings.type = target_type;

  switch (target_type) {
    case TargetType::GATE: {
      if (parts.size() < 2) {
        logger.exitWithError("Target type 'gate' requires a gate name.");
      }
      auto gate_type = convertFromString<GateType>(parts[1]);
      target_settings.gate_type = gate_type;

      if (gate_type == GateType::FILE) {
        if (parts.size() < 3) {
          logger.exitWithError("Gate type 'file' requires a filename.");
        }
        target_settings.gate_file = parts[2];
      }
      break;
    }
    case TargetType::PURE:
      target_settings.levels = std::vector<size_t>{};
      for (size_t i = 1; i < parts.size(); ++i) {
        target_settings.levels->push_back(convertFromString<int>(parts[i]));
      }
      break;
    case TargetType::FROMFILE:
      if (parts.size() < 2) {
        logger.exitWithError("Target type 'file' requires a filename.");
      }
      target_settings.file = parts[1];
      break;
  }

  return target_settings;
}

template <>
std::vector<PiPulseData> CfgParser::convertFromString<std::vector<PiPulseData>>(const std::string& str) {
  auto parts = split(str);
  if (parts.size() % 4 != 0) {
    logger.exitWithError("PiPulse vector requires multiples of 4 parameters: oscil_id, tstart, tstop, amp");
  }

  std::vector<PiPulseData> configs;
  for (size_t i = 0; i < parts.size(); i += 4) {
    PiPulseData config;
    config.oscil_id = convertFromString<size_t>(parts[i]);
    config.tstart = convertFromString<double>(parts[i + 1]);
    config.tstop = convertFromString<double>(parts[i + 2]);
    config.amp = convertFromString<double>(parts[i + 3]);
    configs.push_back(config);
  }

  return configs;
}

template <>
std::vector<ControlSegmentData> CfgParser::convertFromString<std::vector<ControlSegmentData>>(const std::string& str) {
  const auto parts = split(str);

  std::vector<ControlSegmentData> segments;
  size_t i = 0;

  while (i < parts.size()) {
    if (!isValidControlType(parts[i])) {
      logger.exitWithError("Expected control type, got: " + parts[i]);
    }

    ControlSegmentData segment;
    segment.control_type = convertFromString<ControlType>(parts[i++]);

    // Parse parameters until next ControlType or end
    while (i < parts.size() && !isValidControlType(parts[i])) {
      segment.parameters.push_back(convertFromString<double>(parts[i++]));
    }

    // Validate minimum parameter count
    size_t min_params = 1;
    switch (segment.control_type) {
      case ControlType::STEP:
        min_params = 3; // step_amp1, step_amp2, tramp
        break;
      case ControlType::BSPLINE:
      case ControlType::BSPLINE0:
        min_params = 1; // num_basis_functions
        break;
      case ControlType::BSPLINEAMP:
        min_params = 2; // nspline, scaling
        break;
      case ControlType::NONE:
        logger.exitWithError("Control segment type NONE is not valid for configuration.");
        break;
    }

    if (segment.parameters.size() < min_params) {
      logger.exitWithError("Control type requires at least " + std::to_string(min_params) + " parameters, got " +
                           std::to_string(segment.parameters.size()));
    }

    segments.push_back(segment);
  }

  return segments;
}

template <>
std::vector<ControlInitializationData> CfgParser::convertFromString<std::vector<ControlInitializationData>>(
    const std::string& str) {
  const auto parts = split(str);

  std::vector<ControlInitializationData> initializations;
  size_t i = 0;

  while (i < parts.size()) {
    std::string type_str = parts[i++];

    // Validate minimum parameter count
    if (parts.size() <= 1) {
      logger.exitWithError("Expected control_initialization to have a type and at least one parameter.");
    }

    ControlInitializationData initialization;

    auto type_enum = parseEnum(type_str, CONTROL_SEGMENT_INIT_TYPE_MAP);
    if (!type_enum.has_value()) {
      logger.exitWithError("Expected control initialization type (file, constant, random), got: " + type_str);
    }

    switch (type_enum.value()) {
      case ControlSegmentInitType::FILE: {
        initialization.filename = parts[i];
        initializations.push_back(initialization);
        return initializations;
      }
      case ControlSegmentInitType::CONSTANT:
      case ControlSegmentInitType::RANDOM: {
        initialization.init_seg_type = type_enum.value();
        initialization.amplitude = convertFromString<double>(parts[i++]);
        if (i < parts.size() && !isValidControlSegmentInitType(parts[i])) {
          initialization.phase = convertFromString<double>(parts[i++]);
        }
        initializations.push_back(initialization);
        break;
      }
    }
  }

  return initializations;
}
