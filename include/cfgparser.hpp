#pragma once

#include <petsc.h>

#include <cstddef>
#include <functional>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "defs.hpp"
#include "mpi_logger.hpp"
#include "util.hpp"

struct PiPulseData {
  size_t oscil_id; ///< Oscillator ID
  double tstart; ///< Start time
  double tstop; ///< Stop time
  double amp; ///< Amplitude
};

struct ControlSegmentData {
  ControlType control_type; ///< Type of control segment
  std::vector<double> parameters; ///< Parameters for control segment
};

struct ControlInitializationData {
  ControlInitializationType init_seg_type; ///< Type of initialization per segment
  std::optional<double> amplitude; ///< Initial amplitude
  std::optional<double> phase; ///< Initial phase (optional)
  std::optional<std::string> filename; ///< Filename (for file init - mutually exclusive with above)
};

/**
 * @brief Configuration settings passed to Config constructor.
 *
 * Contains all optional configuration parameters that can be provided
 * to configure a Config object. Used by CfgParser to pass settings.
 */
struct ParsedConfigData {
  // General parameters
  std::optional<std::vector<size_t>> nlevels;
  std::optional<std::vector<size_t>> nessential;
  std::optional<size_t> ntime;
  std::optional<double> dt;
  std::optional<std::vector<double>> transfreq;
  std::optional<std::vector<double>> selfkerr;
  std::optional<std::vector<double>> crosskerr;
  std::optional<std::vector<double>> Jkl;
  std::optional<std::vector<double>> rotfreq;
  std::optional<LindbladType> collapse_type;
  std::optional<std::vector<double>> decay_time;
  std::optional<std::vector<double>> dephase_time;
  std::optional<InitialCondition> initialcondition;
  std::optional<std::vector<PiPulseData>> apply_pipulse;
  std::optional<std::string> hamiltonian_file_Hsys;
  std::optional<std::string> hamiltonian_file_Hc;

  // Control and optimization parameters
  std::optional<std::map<int, std::vector<ControlSegmentData>>> indexed_control_segments;
  std::optional<bool> control_enforceBC;
  std::optional<std::map<int, std::vector<ControlInitializationData>>> indexed_control_init;
  std::optional<std::map<int, std::vector<double>>> indexed_control_bounds;
  std::optional<std::map<int, std::vector<double>>> indexed_carrier_frequencies;
  std::optional<OptimTargetSettings> optim_target;
  std::optional<std::vector<double>> gate_rot_freq;
  std::optional<ObjectiveType> optim_objective;
  std::optional<std::vector<double>> optim_weights;
  std::optional<double> optim_tol_grad_abs;
  std::optional<double> optim_tol_grad_rel;
  std::optional<double> optim_tol_finalcost;
  std::optional<double> optim_tol_infidelity;
  std::optional<size_t> optim_maxiter;
  std::optional<double> optim_regul;
  std::optional<double> optim_penalty;
  std::optional<double> optim_penalty_param;
  std::optional<double> optim_penalty_dpdm;
  std::optional<double> optim_penalty_energy;
  std::optional<double> optim_penalty_variation;
  std::optional<bool> optim_regul_tik0;
  std::optional<bool> optim_regul_interpolate; // deprecated

  // Output parameters
  std::optional<std::string> datadir;
  std::optional<std::map<int, std::vector<OutputType>>> indexed_output;
  std::optional<size_t> output_frequency;
  std::optional<size_t> optim_monitor_frequency;
  std::optional<RunType> runtype;
  std::optional<bool> usematfree;
  std::optional<LinearSolverType> linearsolver_type;
  std::optional<size_t> linearsolver_maxiter;
  std::optional<TimeStepperType> timestepper_type;
  std::optional<int> rand_seed;
};

namespace {
// Vector detection
template <typename T>
struct is_vector : std::false_type {};
template <typename T>
struct is_vector<std::vector<T>> : std::true_type {};
template <typename T>
inline constexpr bool is_vector_v = is_vector<T>::value;

// Helper for static_assert(false) in dependent contexts
template <typename>
inline constexpr bool always_false_v = false;
} // namespace

/**
 * @brief Configuration parser that converts raw config to validated Config
 *
 * Handles the logic for applying defaults, parsing mixed-type strings,
 * and validating configuration parameters.
 * TODO cfg: delete this class when .cfg format is removed.
 */
class CfgParser {
 private:
  std::unordered_map<std::string, std::function<void(const std::string&)>> setters; ///< Setters from config string
  std::unordered_map<std::string, std::function<void(int, const std::string&)>>
      indexed_setters; ///< Setters for indexed config strings

  const MPILogger& logger;

  // Configuration settings storage
  ParsedConfigData settings; ///< All configuration settings in one place
  std::optional<bool> optim_regul_interpolate; ///< Deprecated version of optim_regul_tik0

 public:
  CfgParser(const MPILogger& logger);
  ParsedConfigData parseFile(const std::string& filename);
  ParsedConfigData parseString(const std::string& config_content);

 private:
  std::vector<std::string> split(const std::string& str, char delimiter = ',');
  void applyConfigLine(const std::string& line);
  bool handleIndexedSetting(const std::string& key, const std::string& value);

  template <typename StreamType>
  void loadFromStream(StreamType& stream) {
    std::string line;
    while (getline(stream, line)) {
      applyConfigLine(line);
    }
  }

  std::vector<std::vector<double>> convertIndexedToVectorVector(const std::map<int, std::vector<double>>& indexed_map, size_t num_oscillators);

  template <typename T>
  void registerConfig(const std::string& key, std::optional<T>& member) {
    setters[key] = [this, &member](const std::string& value) { member = convertFromString<T>(value); };
  }

  template <typename T>
  void registerIndexedConfig(const std::string& base_key, std::optional<std::map<int, T>>& storage) {
    indexed_setters[base_key] = [this, &storage](int index, const std::string& value) {
      if (!storage.has_value()) {
        storage = std::map<int, T>{};
      }
      (*storage)[index] = convertFromString<T>(value);
    };
  }

  template <typename T>
  T convertFromString(const std::string& str) {
    if constexpr (std::is_same_v<T, std::string>) {
      return str;
    } else if constexpr (std::is_same_v<T, bool>) {
      const std::set<std::string> trueValues = {"true", "yes", "1"};
      std::string lowerStr = toLower(str);
      return trueValues.find(lowerStr) != trueValues.end();
    } else if constexpr (std::is_same_v<T, int>) {
      return std::stoi(str);
    } else if constexpr (std::is_same_v<T, size_t>) {
      return static_cast<size_t>(std::stoul(str));
    } else if constexpr (std::is_same_v<T, double>) {
      return std::stod(str);
    } else if constexpr (is_vector_v<T>) {
      return parseVector<T>(str);
    } else {
      static_assert(always_false_v<T>, "Unsupported type for convertFromString");
    }
  }

  template <typename VectorType>
  VectorType parseVector(const std::string& str) {
    using ElementType = typename VectorType::value_type;
    VectorType vec;
    auto parts = split(str);
    vec.reserve(parts.size());
    for (const auto& part : parts) {
      vec.push_back(convertFromString<ElementType>(part));
    }
    return vec;
  }
};
