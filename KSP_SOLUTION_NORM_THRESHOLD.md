# KSP Solution Norm Threshold Feature

## Overview

A new termination criterion has been added to the Gauss-Newton KSP solver that monitors the solution norm and terminates the solver early when it exceeds a specified threshold.

## Motivation

This feature helps prevent the KSP solver from:
- Diverging to unreasonably large solutions
- Wasting computational effort on solutions that are clearly unsuitable
- Producing numerical instabilities in the outer Gauss-Newton optimization loop

## Implementation Details

### Configuration Parameter

**Parameter name:** `gn_ksp_solution_norm_threshold`

**Location in TOML:** `[optimization]` section

**Default value:** `1e12` (effectively disabled)

**Typical usage value:** `0.1` (as requested)

### How It Works

1. At each KSP iteration, the solver checks:
   - Standard convergence criteria (residual tolerance, max iterations)
   - **NEW:** Solution norm threshold

2. If `||solution||_2 > threshold`, the solver terminates with:
   - Convergence reason: `KSP_DIVERGED_DTOL` (divergence tolerance exceeded)
   - Diagnostic message printed to screen and convergence history file

3. The convergence test is implemented using PETSc's `KSPSetConvergenceTest` API

### Files Modified

1. **[include/config_defaults.hpp](include/config_defaults.hpp)**
   - Added `OPTIM_GN_KSP_SOLUTION_NORM_THRESHOLD = 1e12`

2. **[include/config.hpp](include/config.hpp)**
   - Added member variable `optim_gn_ksp_solution_norm_threshold`
   - Added getter `getOptimGnKspSolutionNormThreshold()`

3. **[src/config.cpp](src/config.cpp)**
   - Added parsing logic with validation (`greaterThan(0.0)`)
   - Added printing in config log output

4. **[include/optimproblem.hpp](include/optimproblem.hpp)**
   - Added member variable `ksp_solution_norm_threshold`
   - Added getter `getKspSolutionNormThreshold()`

5. **[src/optimproblem.cpp](src/optimproblem.cpp)**
   - Added `KSPConvergenceTestSolutionNorm()` callback function
   - Initialize threshold from config in constructor
   - Set custom convergence test in `solveGaussNewtonKSP()`

6. **[config_template.toml](config_template.toml)**
   - Added documentation and example usage

## Usage Example

```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_type = "cg"
gn_ksp_rtol = 1e-2
gn_ksp_maxiter = 100
gn_ksp_damping = 1e-3
gn_ksp_solution_norm_threshold = 0.1  # Terminate if solution norm exceeds 0.1
```

## Output

When the threshold is triggered, you'll see:

**Console output:**
```
KSP terminated: solution norm 1.234567e-01 exceeds threshold 1.000000e-01
```

**Convergence history file (`ksp_convergence_*.dat`):**
```
# KSP terminated: solution norm 1.234567e-01 exceeds threshold 1.000000e-01
```

## Testing

See example configuration files:
- [example_gn_ksp_damping.toml](example_gn_ksp_damping.toml) - Updated with new parameter
- [example_gn_solution_norm_threshold.toml](example_gn_solution_norm_threshold.toml) - Dedicated example

## Technical Notes

### Interaction with Other Convergence Criteria

The solution norm check is performed **after** the default convergence tests. This means:

1. If KSP has already converged or diverged for another reason, that reason takes precedence
2. The solution norm check only applies when KSP would otherwise continue iterating
3. This prevents false positives from initial iterations with large solutions

### Custom Convergence Test Function

```cpp
PetscErrorCode KSPConvergenceTestSolutionNorm(KSP ksp, PetscInt it, 
                                               PetscReal rnorm, 
                                               KSPConvergedReason *reason, 
                                               void *ctx)
```

This function:
- First calls `KSPConvergedDefault()` to apply standard tests
- Then checks solution norm only if still iterating
- Uses `KSP_DIVERGED_DTOL` reason code for exceeded threshold

### Why Large Default Value?

The default of `1e12` effectively disables this feature by default because:
- Typical solution norms in quantum control problems are O(1e-3) to O(1e1)
- A value of `1e12` won't be reached in practice
- This preserves backward compatibility
- Users must opt-in by setting a reasonable threshold

## Choosing the Threshold

Guidelines for selecting an appropriate threshold value:

1. **Too small (< 0.01):** May terminate prematurely during normal convergence
2. **Too large (> 100):** Won't catch divergence until it's severe
3. **Recommended range:** 0.1 to 10.0
4. **Problem-dependent:** Inspect typical solution norms in converged cases

## Related Parameters

This feature works in conjunction with:
- `gn_ksp_damping`: Damping parameter that helps stabilize the system
- `gn_ksp_rtol`: Relative tolerance for residual norm
- `gn_ksp_maxiter`: Maximum iterations (another safety mechanism)

For best results, tune these parameters together based on your specific problem.
