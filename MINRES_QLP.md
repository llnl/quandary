# MINRES-QLP for Gauss-Newton System

## Overview

MINRES-QLP (Minimum Residual with QLP factorization) is a variant of the MINRES iterative solver that is more robust for ill-conditioned or nearly singular symmetric systems. This feature leverages PETSc's built-in MINRES-QLP implementation.

## What is MINRES-QLP?

MINRES-QLP improves upon standard MINRES by using a QLP factorization of the tridiagonal matrix generated during the Lanczos process, rather than the standard LQ factorization. This provides:

1. **Better numerical stability** for ill-conditioned problems
2. **More accurate residual norms** throughout the iteration
3. **Earlier detection of convergence** or stagnation
4. **Better handling** of nearly singular systems

**Reference:** Choi, Paige, and Saunders (2011), "MINRES-QLP: A Krylov subspace method for indefinite or singular symmetric systems", SIAM Journal on Scientific Computing, 33(4), 1810-1836.

## When to Use MINRES-QLP

Use MINRES-QLP instead of standard MINRES when:

- The Gauss-Newton Hessian approximation becomes ill-conditioned during optimization
- Standard MINRES shows slow or erratic convergence
- The optimization is near a saddle point or singular region
- You want more robust behavior at the cost of slightly more computation per iteration

## Configuration

### Option 1: TOML Configuration File

Add to your `[optimization]` section:

```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_type = "minres"       # Must be "minres"
gn_minres_qlp = true         # Enable QLP variant
gn_ksp_rtol = 1e-2
gn_ksp_maxiter = 100
```

### Option 2: Command-Line Override

You can enable MINRES-QLP without modifying the config file:

```bash
./quandary config.toml -gn_ksp_minres_qlp true
```

Or combine with other KSP options:

```bash
./quandary config.toml -gn_ksp_type minres -gn_ksp_minres_qlp true -gn_ksp_rtol 1e-3
```

### Option 3: PETSc Direct Option

PETSc's command-line option also works:

```bash
./quandary config.toml -gn_ksp_minres_qlp true
```

This is the same as the TOML `gn_minres_qlp = true` setting.

## Implementation Details

### How It Works

1. When `gn_minres_qlp = true` and `gn_ksp_type = "minres"`, the code sets:
   ```
   -gn_ksp_minres_qlp true
   ```
   in PETSc's options database before calling `KSPSetFromOptions()`.

2. PETSc's MINRES solver then internally switches to the QLP variant.

3. The KSP convergence history and monitoring remain the same - you still get:
   - Iteration-by-iteration residual norms
   - Solution norm tracking
   - Convergence files (`ksp_convergence_minres_iterXXXX.dat`)

### Configuration Hierarchy

Settings are applied in this order (later ones override earlier):

1. **TOML config file defaults** (e.g., `gn_minres_qlp = true`)
2. **Command-line options** (e.g., `-gn_ksp_minres_qlp false`)
3. **PETSc environment** (e.g., `PETSC_OPTIONS="-gn_ksp_minres_qlp true"`)

This allows you to:
- Set defaults in the TOML file for reproducibility
- Override on the command line for experimentation
- Use environment variables for cluster/batch job settings

## Performance Comparison

### Standard MINRES vs MINRES-QLP

| Aspect | Standard MINRES | MINRES-QLP |
|--------|----------------|------------|
| **Cost per iteration** | Lower | ~10-20% higher |
| **Ill-conditioning tolerance** | Moderate | Excellent |
| **Convergence reliability** | Good | Better |
| **Memory** | Similar | Similar |
| **Total iterations** | May need more | Often fewer |

For well-conditioned problems, standard MINRES is slightly faster. For ill-conditioned problems, MINRES-QLP often wins by taking fewer iterations.

## Example: Comparing Standard vs QLP

### Run with Standard MINRES

```bash
./quandary config.toml
```

or explicitly:

```toml
gn_ksp_type = "minres"
gn_minres_qlp = false  # Standard MINRES
```

### Run with MINRES-QLP

```bash
./quandary config.toml -gn_ksp_minres_qlp true
```

or in TOML:

```toml
gn_ksp_type = "minres"
gn_minres_qlp = true  # QLP variant
```

### Compare Convergence

Check the convergence history files:
```bash
# Standard MINRES
cat ksp_convergence_minres_iter0000.dat

# With QLP
cat ksp_convergence_minres_iter0000.dat  # Same filename, different behavior
```

Look for:
- **Fewer iterations** to converge
- **Smoother residual decay** (less erratic behavior)
- **Better final accuracy** for the same tolerance

## Diagnostic Output

When MINRES-QLP is enabled, you'll see:

```
Enabling MINRES-QLP variant for Gauss-Newton KSP solver
```

during initialization. Then during solving:

```
KSP it 0: residual norm = 1.234567e-02, solution norm = 0.000000e+00
KSP it 1: residual norm = 5.678901e-03, solution norm = 1.234567e-03
...
Gauss-Newton minres stats: iterations = 25, MatVec counter = 25, residual norm = 9.876543e-04
```

The output format is identical to standard MINRES - only the internal algorithm differs.

## Troubleshooting

### MINRES-QLP Not Available

If you see an error like:
```
Unknown option: -gn_ksp_minres_qlp
```

Your PETSc version may not support MINRES-QLP. The option was added in PETSc 3.4. Update PETSc or use standard MINRES.

### No Improvement vs Standard MINRES

If MINRES-QLP doesn't help:

1. **Try a preconditioner** - ill-conditioning may need preconditioning:
   ```toml
   gn_pc_type = "jacobi"  # or "ilu"
   ```

2. **Increase damping** - the Gauss-Newton matrix may need more regularization:
   ```toml
   gn_ksp_damping = 1e-2  # Larger than default 1e-3
   ```

3. **Check solution norm threshold** - solutions may be diverging:
   ```toml
   gn_ksp_solution_norm_threshold = 0.1
   ```

### QLP Slower Than Expected

MINRES-QLP has ~10-20% overhead per iteration. If you're seeing much worse performance:

1. Check that you're comparing total solve time, not per-iteration time
2. MINRES-QLP may take fewer iterations, compensating for the per-iteration cost
3. For well-conditioned problems, standard MINRES is faster - use QLP only when needed

## Related Parameters

MINRES-QLP works well with these other settings:

```toml
[optimization]
solver_type = "gauss_newton"
gn_ksp_type = "minres"
gn_minres_qlp = true                     # Enable QLP

# These parameters work with both standard MINRES and MINRES-QLP:
gn_ksp_rtol = 1e-2                       # Convergence tolerance
gn_ksp_maxiter = 100                     # Max iterations
gn_ksp_damping = 1e-3                    # Regularization
gn_ksp_solution_norm_threshold = 0.1     # Divergence detection
gn_ksp_warmstart = true                  # Reuse previous solution
```

## Files Modified

- **[include/config_defaults.hpp](include/config_defaults.hpp)** - Added `OPTIM_GN_MINRES_QLP = false`
- **[include/config.hpp](include/config.hpp)** - Added `optim_gn_minres_qlp` member and getter
- **[src/config.cpp](src/config.cpp)** - Added parsing and printing logic
- **[src/optimproblem.cpp](src/optimproblem.cpp)** - Added PETSc option setting
- **[config_template.toml](config_template.toml)** - Added documentation
- **[example_minresqlp.toml](example_minresqlp.toml)** - Example configuration

## See Also

- PETSc MINRES documentation: https://petsc.org/release/manualpages/KSP/KSPMINRES/
- Original MINRES-QLP paper: https://doi.org/10.1137/100787921
- [KSP_SOLUTION_NORM_THRESHOLD.md](KSP_SOLUTION_NORM_THRESHOLD.md) - Solution norm monitoring
- [config_template.toml](config_template.toml) - Full configuration reference
