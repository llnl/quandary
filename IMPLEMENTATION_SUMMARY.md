# Implementation Summary: Gauss-Newton KSP Enhancements

This document summarizes all the enhancements made to the Gauss-Newton KSP solver in Quandary.

## Features Implemented

### 1. Tunable KSP Damping Parameter
**Status:** ✅ Complete

The `ksp_damping` parameter was previously hardcoded to `1e-3`. It is now configurable via TOML.

**Configuration:**
```toml
[optimization]
gn_ksp_damping = 1e-3  # Default: 1e-3
```

**Purpose:** Controls the regularization `(J^T J + damping * I)` applied to stabilize the Gauss-Newton Hessian approximation.

---

### 2. Solution Norm Threshold Termination
**Status:** ✅ Complete

A new termination criterion monitors the KSP solution norm and stops if it exceeds a threshold.

**Configuration:**
```toml
[optimization]
gn_ksp_solution_norm_threshold = 0.1  # Default: 1e12 (disabled)
```

**Purpose:** Prevents divergence by detecting when the solution is growing unreasonably large. When `||solution|| > threshold`, the KSP solver terminates with reason `KSP_DIVERGED_DTOL`.

**Output:**
```
KSP terminated: solution norm 1.234e-01 exceeds threshold 1.000e-01
```

---

### 3. MINRES-QLP Variant
**Status:** ✅ Complete

Enables PETSc's MINRES-QLP variant, which is more robust than standard MINRES for ill-conditioned systems.

**Configuration:**
```toml
[optimization]
gn_ksp_type = "minres"
gn_minres_qlp = true  # Default: false
```

**Alternatives:**
- Command line: `./quandary config.toml -gn_ksp_minres_qlp true`
- PETSc option: `PETSC_OPTIONS="-gn_ksp_minres_qlp true"`

**Purpose:** Uses QLP factorization instead of LQ factorization in the Lanczos process, providing better numerical stability for ill-conditioned Gauss-Newton systems.

---

## Configuration Summary

All new parameters are in the `[optimization]` section:

```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_type = "minres"

# New parameters:
gn_ksp_damping = 1e-3                     # Regularization parameter
gn_ksp_solution_norm_threshold = 0.1      # Divergence detection
gn_minres_qlp = true                      # Enable QLP variant

# Existing parameters (unchanged):
gn_ksp_rtol = 1e-3
gn_ksp_maxiter = 100
gn_ksp_warmstart = false
gn_pc_type = "none"
```

---

## Files Modified

### Configuration System
- `include/config_defaults.hpp` - Added 3 new defaults
- `include/config.hpp` - Added 3 member variables and getters
- `src/config.cpp` - Added parsing, validation, and printing logic

### Optimization Problem
- `include/optimproblem.hpp` - Added member variable and getter for threshold
- `src/optimproblem.cpp` - 
  - Initialize all parameters from config
  - Custom convergence test for solution norm monitoring
  - Enable MINRES-QLP via PETSc options

### Documentation
- `config_template.toml` - Documented all new parameters
- `example_gn_ksp_damping.toml` - Example showing damping and threshold
- `example_gn_solution_norm_threshold.toml` - Focused threshold example
- `example_minresqlp.toml` - MINRES-QLP configuration example
- `KSP_SOLUTION_NORM_THRESHOLD.md` - Detailed threshold documentation
- `MINRES_QLP.md` - Comprehensive MINRES-QLP guide
- `IMPLEMENTATION_SUMMARY.md` - This file

---

## Usage Examples

### Example 1: Basic Damping Configuration
```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_damping = 5e-4  # Reduce from default 1e-3 for less regularization
```

### Example 2: Divergence Detection
```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_solution_norm_threshold = 0.1  # Stop if solution norm > 0.1
```

### Example 3: Robust Ill-Conditioned Solving
```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_type = "minres"
gn_minres_qlp = true                    # Use QLP variant
gn_ksp_damping = 1e-2                   # Increased regularization
gn_ksp_solution_norm_threshold = 1.0    # Safety threshold
```

### Example 4: Command-Line Override
```bash
# Enable MINRES-QLP without editing config
./quandary config.toml -gn_ksp_minres_qlp true

# Try different damping values
./quandary config.toml -gn_ksp_damping 1e-4
./quandary config.toml -gn_ksp_damping 1e-2

# Set solution norm threshold
./quandary config.toml -gn_ksp_solution_norm_threshold 0.5
```

---

## Backward Compatibility

All features are **backward compatible**:

- Default values preserve existing behavior
- Old config files work without modification
- New parameters are optional
- Command-line options don't break existing workflows

**Migration:** No changes needed to existing config files. New features are opt-in.

---

## Testing Recommendations

### Test 1: Verify Damping Works
```bash
# Run with different damping values
./quandary config.toml -gn_ksp_damping 1e-4  # Less regularization
./quandary config.toml -gn_ksp_damping 1e-2  # More regularization

# Compare convergence histories
diff ksp_convergence_*.dat
```

### Test 2: Verify Solution Norm Threshold
```bash
# Set a low threshold that should trigger
./quandary config.toml -gn_ksp_solution_norm_threshold 0.01

# Check for termination message:
# "KSP terminated: solution norm X.XXe-XX exceeds threshold X.XXe-XX"
```

### Test 3: Compare MINRES vs MINRES-QLP
```bash
# Standard MINRES
./quandary config.toml -gn_ksp_type minres

# MINRES-QLP
./quandary config.toml -gn_ksp_type minres -gn_ksp_minres_qlp true

# Compare:
# - Total KSP iterations
# - Convergence behavior (residual decay)
# - Final optimization result
```

---

## Performance Impact

### Damping Parameter
- **Overhead:** None (same computation, different value)
- **Effect:** Can improve or worsen convergence depending on problem

### Solution Norm Threshold
- **Overhead:** Negligible (~1 extra norm computation per KSP iteration)
- **Effect:** Prevents wasted computation on diverging solutions

### MINRES-QLP
- **Overhead:** ~10-20% per KSP iteration
- **Effect:** Often takes fewer total iterations for ill-conditioned problems
- **Net impact:** Problem-dependent; can be faster for difficult problems

---

## Future Enhancements

Potential future work:

1. **Adaptive damping:** Automatically adjust `ksp_damping` during optimization
2. **Multiple thresholds:** Separate thresholds for warning vs termination
3. **Preconditioner support:** Add preconditioner options for MINRES-QLP
4. **Performance profiling:** Built-in timing for KSP vs MatVec operations
5. **Auto-selection:** Automatically choose between MINRES and MINRES-QLP based on condition number

---

## References

### MINRES-QLP
- Choi, S.-C., Paige, C. C., & Saunders, M. A. (2011). "MINRES-QLP: A Krylov subspace method for indefinite or singular symmetric systems." SIAM Journal on Scientific Computing, 33(4), 1810-1836.
- PETSc MINRES documentation: https://petsc.org/release/manualpages/KSP/KSPMINRES/

### Gauss-Newton Methods
- Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (2nd ed.). Springer.
- Chapter on Trust Region Methods and Gauss-Newton approximations

---

## Build Status

✅ **Compiles successfully** with no errors (verified 2026-09-23)

**Build command:**
```bash
cd build && make -j4
```

**Test availability:**
- Unit tests: Not yet implemented for new features
- Integration tests: Manual testing recommended
- Example configs: Provided in `example_*.toml` files
