#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <toml++/toml.hpp>
#include <vector>

/**
 * @brief Validation utilities for TOML configuration parsing.
 *
 * Provides a chainable API for type-safe TOML field extraction and validation. Chain validation
 * methods together, then extract the final value in one operation with clear error messages.
 * 
 * ## Basic Usage 
 *
 * 1. Create a validator for option "key" within toml_table using
 *      - `validators::field<T>(toml_table, "key")`: for scalar fields
 *      - `validators::vectorField<T>(toml_table, "key")`: for vector fields
 * 2. Append chainable validation methods from below as needed:
 *    For scalar fields:
 *      - `greaterThan(value)`: Field must be > value
 *      - `greaterThanEqual(value)`: Field must be >= value
 *      - `lessThan(value)`: Field must be < value
 *      - `positive()`: Field must be > 0 (shorthand for greaterThan(0))
 *    For vector fields:
 *      - `minLength(size)`: Vector must have at least size elements
 *      - `hasLength(size)`: Vector must have exactly size elements
 *      - `positive()`: All elements must be > 0
 * 3. Extract the value by appending 
 *         .value()           - for required fields
 *         .valueOr(default)  - for optional fields with default
 * 
 * ## Examples
 * @code  
 * // Parse a required positive double with key "step_size".
 * double step_size = validators::field<double>(toml, "step_size").positive().value();
 *
 * // Parse an optional string with key "output_dir" with default value "./data_out"
 * std::string output_dir = validators::field<std::string>(toml, "output_dir").valueOr("./data_out");
 *
 * // Parse a required vector of positive doubles with key "frequencies", defaulting to {1.0, 2.0}
 * std::vector<double> frequencies = validators::vectorField<double>(toml, "frequencies").minLength(1).positive().valueOr(std::vector<double>{1.0, 2.0});
 * 
 * @endcode
 *
 * ## Supported Types
 *
 * Scalar fields: int, size_t, double, std::string, bool, etc.
 * Vector fields: std::vector<T> where T is any supported scalar type
 */
namespace validators {

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
 * @brief Exception thrown when configuration validation fails.
 */
class ValidationError : public std::runtime_error {
 public:
  ValidationError(const std::string& field, const std::string& message)
      : std::runtime_error("Validation error for field '" + field + "': " + message) {}
};

/**
 * @brief Chainable validator for scalar TOML fields.
 *
 * This class provides type-safe extraction and validation of single-value fields from
 * TOML configuration using the method chaining pattern described above.
 *
 * ## Available Validation Methods
 *
 * If any of the following methods are called, the corresponding validation is applied
 * when extracting the value, and an error is thrown if the validation fails.
 *
 * - `greaterThan(value)`: Field must be > value
 * - `greaterThanEqual(value)`: Field must be >= value
 * - `lessThan(value)`: Field must be < value
 * - `positive()`: Field must be > 0 (shorthand for greaterThan(0))
 * 
 *  ## Extraction Methods
 *
 * - `value()`: Extract scalar or throws an error if field is missing or invalid
 * - `valueOr(default)`: Extract scalar or return default if field missing
 *
 * @tparam T Type of field to validate (int, double, string, bool, etc.)
 */
template <typename T>
class Validator {
 private:
  const toml::table& config;
  std::string key;
  std::optional<T> greater_than;
  std::optional<T> greater_than_equal;
  std::optional<T> less_than;

 public:
  Validator(const toml::table& config_, const std::string& key_) : config(config_), key(key_) {}

  /**
   * @brief Requires value to be strictly greater than threshold.
   *
   * @param greater_than_ Threshold value (exclusive)
   * @return Reference to this validator for chaining
   */
  Validator& greaterThan(T greater_than_) {
    greater_than = greater_than_;
    return *this;
  }

  /**
   * @brief Requires value to be greater than or equal to threshold.
   *
   * @param greater_than_equal_ Threshold value (inclusive)
   * @return Reference to this validator for chaining
   */
  Validator& greaterThanEqual(T greater_than_equal_) {
    greater_than_equal = greater_than_equal_;
    return *this;
  }

  /**
   * @brief Requires value to be strictly less than threshold.
   *
   * @param less_than_ Threshold value (exclusive)
   * @return Reference to this validator for chaining
   */
  Validator& lessThan(T less_than_) {
    less_than = less_than_;
    return *this;
  }

  /**
   * @brief Requires value to be strictly positive (> 0).
   *
   * @return Reference to this validator for chaining
   */
  Validator& positive() {
    greaterThan(T{0});
    return *this;
  }

 private:
  std::optional<T> extractValue() {
    // If key doesn't exist, return nullopt
    if (!config.contains(key)) {
      return std::nullopt;
    }

    // Key exists, try to extract value with type checking
    auto val = config[key].template value<T>();
    if (!val) {
      // Key exists but wrong type - always an error
      throw ValidationError(key, "wrong type (expected " + getTypeName<T>() + ")");
    }

    return val;
  }

  T validateValue(T result) {
    if (greater_than && result <= *greater_than) {
      std::ostringstream oss;
      oss << "must be > " << *greater_than << ", got " << result;
      throw ValidationError(key, oss.str());
    }

    if (greater_than_equal && result < *greater_than_equal) {
      std::ostringstream oss;
      oss << "must be >= " << *greater_than_equal << ", got " << result;
      throw ValidationError(key, oss.str());
    }

    if (less_than && result >= *less_than) {
      std::ostringstream oss;
      oss << "must be < " << *less_than << ", got " << result;
      throw ValidationError(key, oss.str());
    }

    return result;
  }

 public:
  /**
   * @brief Extracts and validates the field value.
   *
   * Throws ValidationError if the field is missing
   * or if any validation rules fail.
   *
   * @return The validated field value
   * @throws ValidationError If validation fails
   */
  T value() {
    auto val = extractValue();

    if (!val) {
      throw ValidationError(key, "field not found");
    }

    return validateValue(*val);
  }

  /**
   * @brief Extracts field value or returns default if missing.
   *
   * If the field is present, validates it (may throw). If absent,
   * returns the provided default without validation.
   *
   * @param default_value_ Default value to use if field is missing
   * @return The field value or default
   * @throws ValidationError If field exists but validation fails
   */
  T valueOr(T default_value_) {
    auto val = extractValue();
    if (!val) return default_value_; // Key doesn't exist - use default

    return validateValue(*val); // Key exists - validate it (will throw on wrong type)
  }
};

/**
 * @brief Chainable validator for vector TOML fields.
 *
 * This class provides type-safe extraction and validation of array fields from
 * TOML configuration using the method chaining pattern described above.
 *
 * ## Available Validation Methods
 *
 * - `minLength(size)`: Vector must have at least `size` elements
 * - `hasLength(size)`: Vector must have exactly `size` elements
 * - `positive()`: All elements must be > 0 (only for numeric types)
 *
 * ## Extraction Methods
 *
 * - `value()`: Extract validated array (throws if field missing or invalid)
 * - `valueOr(default)`: Extract array or return default if field missing
 *
 * ## Type Requirements
 *
 * The element type T must be a type supported by the TOML library (int, double, string, bool).
 * For numeric types, you can use the positive() validation. For other types, only
 * array-level validations (length) are available.
 *
 * @tparam T Element type of the vector (int, double, string, bool, etc.)
 */
template <typename T>
class VectorValidator {
 private:
  const toml::table& config;
  std::string key;
  std::optional<size_t> min_length;
  std::optional<size_t> exact_length;
  bool is_positive = false;

 public:
  VectorValidator(const toml::table& config_, const std::string& key_) : config(config_), key(key_) {}

  /**
   * @brief Requires minimum vector length.
   *
   * @param min_len_ Minimum number of elements (inclusive)
   * @return Reference to this validator for chaining
   */
  VectorValidator& minLength(size_t min_len_) {
    min_length = min_len_;
    return *this;
  }

  /**
   * @brief Requires all vector elements to be strictly positive (> 0).
   *
   * @return Reference to this validator for chaining
   */
  VectorValidator& positive() {
    is_positive = true;
    return *this;
  }

  /**
   * @brief Requires vector to have exactly the specified length.
   *
   * @param exact_len_ Exact number of elements required
   * @return Reference to this validator for chaining
   */
  VectorValidator& hasLength(size_t exact_len_) {
    exact_length = exact_len_;
    return *this;
  }

 private:
  std::optional<std::vector<T>> extractVector() {
    // If key doesn't exist, return nullopt
    if (!config.contains(key)) {
      return std::nullopt;
    }

    // Key exists, check if it's an array
    auto* arr = config[key].as_array();
    if (!arr) {
      // Key exists but wrong type - always an error
      throw ValidationError(key, "wrong type (expected array)");
    }

    // Extract and validate array elements
    std::vector<T> result;
    for (size_t i = 0; i < arr->size(); ++i) {
      auto val = arr->at(i).template value<T>();
      if (!val) {
        std::ostringstream oss;
        oss << "element [" << i << "] wrong type (expected " << getTypeName<T>() << ")";
        throw ValidationError(key, oss.str());
      }
      result.push_back(*val);
    }

    return result;
  }

  std::vector<T> validateVector(std::vector<T> result) {
    if (min_length && result.size() < *min_length) {
      std::ostringstream oss;
      oss << "must have at least " << *min_length << " elements, got " << result.size();
      throw ValidationError(key, oss.str());
    }

    for (size_t i = 0; i < result.size(); ++i) {
      T& element = result[i];

      if (is_positive && element <= T{0}) {
        std::ostringstream oss;
        oss << "element [" << i << "] must be positive, got " << element;
        throw ValidationError(key, oss.str());
      }
    }

    return result;
  }

 public:
  /**
   * @brief Extracts and validates the vector field.
   *
   * @return The validated vector
   * @throws ValidationError If validation fails
   */
  std::vector<T> value() {
    auto val = extractVector();

    if (!val) {
      throw ValidationError(key, "field not found");
    }

    return validateVector(*val);
  }

  /**
   * @brief Extracts vector or returns default if missing.
   *
   * @param default_value_ Default vector to use if field is missing
   * @return The field vector or default
   * @throws ValidationError If field exists but validation fails
   */
  std::vector<T> valueOr(const std::vector<T>& default_value_) {
    auto val = extractVector();
    if (!val) return default_value_; // Key doesn't exist - use default

    return validateVector(*val); // Key exists - validate it (will throw on wrong type)
  }
};

/**
 * @brief Creates a scalar field validator.
 *
 * Helper function to start a validation chain for scalar fields.
 *
 * @tparam T Type of field to validate
 * @param config_ TOML table containing the field
 * @param key_ Name of the field to validate
 * @return A Validator for chaining validation rules
 */
template <typename T>
Validator<T> field(const toml::table& config_, const std::string& key_) {
  return Validator<T>(config_, key_);
}

/**
 * @brief Creates a vector field validator.
 *
 * Helper function to start a validation chain for array/vector fields.
 *
 * @tparam T Element type of the vector
 * @param config_ TOML table containing the field
 * @param key_ Name of the field to validate
 * @return A VectorValidator for chaining validation rules
 */
template <typename T>
VectorValidator<T> vectorField(const toml::table& config_, const std::string& key_) {
  return VectorValidator<T>(config_, key_);
}

/**
 * @brief Extracts an optional scalar value from a TOML node.
 *
 * Helper for extracting optional scalar values when the validator API doesn't fit well
 * (e.g., nested structures or conditional parsing). If the node doesn't exist or
 * contains a type mismatch, returns nullopt.
 *
 * @tparam T Type of the scalar value
 * @tparam NodeType Type of the TOML node
 * @param node TOML node that may contain a value
 * @return Value if node exists and has matching type, nullopt otherwise
 */
template <typename T, typename NodeType>
std::optional<T> getOptional(const toml::node_view<NodeType>& node) {
  return node.template value<T>();
}

/**
 * @brief Extracts an optional vector from a TOML node.
 *
 * Helper for extracting vectors when the validator API doesn't fit well
 * (e.g., nested structures or conditional parsing). If the node is not
 * an array or contains type mismatches, returns nullopt.
 *
 * @tparam T Element type of the vector
 * @tparam NodeType Type of the TOML node
 * @param node TOML node that may contain an array
 * @return Vector if node is a valid array with matching types, nullopt otherwise
 */
template <typename T, typename NodeType>
std::optional<std::vector<T>> getOptionalVector(const toml::node_view<NodeType>& node) {
  auto* arr = node.as_array();
  if (!arr) return std::nullopt;

  std::vector<T> result;
  for (size_t i = 0; i < arr->size(); ++i) {
    auto val = arr->at(i).template value<T>();
    if (!val) return std::nullopt; // Type mismatch in array element
    result.push_back(*val);
  }

  return result;
}

/**
 * @brief Extracts a required table from a TOML configuration.
 *
 * Validates that the specified key exists and contains a table.
 *
 * @param config Parent TOML table
 * @param key Name of the table field
 * @return Reference to the table
 * @throws ValidationError if table is missing or wrong type
 */
inline const toml::table& getRequiredTable(const toml::table& config, const std::string& key) {
  if (!config.contains(key)) {
    throw ValidationError(key, "table is required");
  }

  auto* table = config[key].as_table();
  if (!table) {
    throw ValidationError(key, "must be a table");
  }

  return *table;
}

/**
 * @brief Extracts an array of tables from a TOML configuration.
 *
 * Returns an empty array if the key doesn't exist or isn't an array of tables.
 *
 * @param config Parent TOML table
 * @param key Name of the array of tables field
 * @return The array of tables if it exists, otherwise an empty array
 */
inline toml::array getArrayOfTables(const toml::table& config, const std::string& key) {
  if (!config.contains(key) || !config[key].is_array_of_tables()) {
    return toml::array{};
  }

  return *config[key].as_array();
}

/**
 * @brief Parses a field that can be either a scalar (applied to all) or an exact-size array.
 *
 * For per-oscillator settings, this allows users to specify either:
 * - A single scalar value: `transfreq = 4.1` (applied to all oscillators)
 * - An array of exact size: `transfreq = [4.1, 4.2, 4.3]` (one per oscillator)
 *
 * @tparam T Element type (int, double, etc.)
 * @param config TOML table containing the field
 * @param key Name of the field
 * @param expected_size Expected array size (e.g., num_oscillators)
 * @return Vector of expected_size elements
 * @throws ValidationError if field is missing, wrong type, or array has wrong size
 */
template <typename T>
std::vector<T> scalarOrVector(const toml::table& config, const std::string& key, size_t expected_size) {
  if (!config.contains(key)) {
    throw ValidationError(key, "field not found");
  }

  if (auto* arr = config[key].as_array()) {
    // Array: must have exact expected size
    if (arr->size() != expected_size) {
      std::ostringstream oss;
      oss << "array must have exactly " << expected_size << " elements, got " << arr->size();
      throw ValidationError(key, oss.str());
    }
    std::vector<T> result;
    for (size_t i = 0; i < arr->size(); ++i) {
      auto val = arr->at(i).template value<T>();
      if (!val) {
        std::ostringstream oss;
        oss << "element [" << i << "] wrong type (expected " << getTypeName<T>() << ")";
        throw ValidationError(key, oss.str());
      }
      result.push_back(*val);
    }
    return result;
  }

  if (auto val = config[key].template value<T>()) {
    // Scalar: fill vector with this value
    return std::vector<T>(expected_size, *val);
  }

  throw ValidationError(key, "must be either a scalar value or an array of " + getTypeName<T>());
}

/**
 * @brief Parses an optional field that can be either a scalar or an exact-size array.
 *
 * Same as scalarOrVector but returns default_value if field is missing.
 *
 * @tparam T Element type (int, double, etc.)
 * @param config TOML table containing the field
 * @param key Name of the field
 * @param expected_size Expected array size (e.g., num_oscillators)
 * @param default_value Default vector to return if field is missing
 * @return Vector of expected_size elements
 * @throws ValidationError if field has wrong type or array has wrong size
 */
template <typename T>
std::vector<T> scalarOrVectorOr(const toml::table& config, const std::string& key, size_t expected_size, const std::vector<T>& default_value) {
  if (!config.contains(key)) {
    return default_value;
  }
  return scalarOrVector<T>(config, key, expected_size);
}

} // namespace validators
