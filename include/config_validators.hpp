#pragma once

#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * @brief Validation utilities for configuration parsing.
 *
 * Provides a chainable API for type-safe value validation. Chain validation
 * methods together, then extract the final value in one operation with clear error messages.
 *
 * ## Basic Usage
 *
 * 1. Create a validator using:
 *      - `validators::field<T>(value, "key")`: for scalar fields
 *      - `validators::vectorField<T>(value, "key")`: for vector fields
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
 * // Validate a required positive double with key "step_size".
 * double step_size = validators::field<double>(step_size_opt, "step_size").positive().value();
 *
 * // Validate an optional string with key "output_dir" with default value "./data_out"
 * std::string output_dir = validators::field<std::string>(output_dir_opt, "output_dir").valueOr("./data_out");
 *
 * // Validate a vector of positive doubles with key "frequencies", defaulting to {1.0, 2.0}
 * std::vector<double> frequencies = validators::vectorField<double>(frequencies_opt, "frequencies").minLength(1).positive().valueOr(std::vector<double>{1.0, 2.0});
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
 * @brief Exception thrown when configuration validation fails.
 */
class ValidationError : public std::runtime_error {
 public:
  ValidationError(const std::string& field, const std::string& message)
      : std::runtime_error("Validation error for field '" + field + "': " + message) {}
};

/**
 * @brief Chainable validator for scalar fields.
 *
 * This class provides type-safe validation of single values using the method chaining
 * pattern described above.
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
  std::optional<T> val;
  std::string key;
  std::optional<T> greater_than;
  std::optional<T> greater_than_equal;
  std::optional<T> less_than;

 public:
  Validator(std::optional<T> val_, const std::string& key_) : val(std::move(val_)), key(key_) {}

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
    if (!val) {
      return default_value_;
    }
    return validateValue(*val);
  }
};

/**
 * @brief Chainable validator for vector fields.
 *
 * This class provides type-safe validation of vector/array fields using the method
 * chaining pattern described above.
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
 * @tparam T Element type of the vector (int, double, string, bool, etc.)
 */
template <typename T>
class VectorValidator {
 private:
  std::optional<std::vector<T>> val;
  std::string key;
  std::optional<size_t> min_length;
  std::optional<size_t> exact_length;
  bool is_positive = false;

 public:
  VectorValidator(std::optional<std::vector<T>> val_, const std::string& key_) : val(std::move(val_)), key(key_) {}

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
  std::vector<T> validateVector(std::vector<T> result) {
    if (exact_length && result.size() != *exact_length) {
      std::ostringstream oss;
      oss << "must have exactly " << *exact_length << " elements, got " << result.size();
      throw ValidationError(key, oss.str());
    }

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
    if (!val) {
      return default_value_;
    }
    return validateVector(*val);
  }
};

/**
 * @brief Creates a scalar field validator.
 *
 * Helper function to start a validation chain for scalar fields.
 *
 * @tparam T Type of field to validate
 * @param val_ Optional value to validate
 * @param key_ Name of the field (for error messages)
 * @return A Validator for chaining validation rules
 */
template <typename T>
Validator<T> field(std::optional<T> val_, const std::string& key_ = "value") {
  return Validator<T>(std::move(val_), key_);
}

/**
 * @brief Creates a vector field validator.
 *
 * Helper function to start a validation chain for array/vector fields.
 *
 * @tparam T Element type of the vector
 * @param val_ Optional vector to validate
 * @param key_ Name of the field (for error messages)
 * @return A VectorValidator for chaining validation rules
 */
template <typename T>
VectorValidator<T> vectorField(std::optional<std::vector<T>> val_, const std::string& key_ = "value") {
  return VectorValidator<T>(std::move(val_), key_);
}

} // namespace validators
