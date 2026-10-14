# Gauss-Newton Solver Calling Sequences

## Overview

Quandary implements three Gauss-Newton solver modes, selected via `solver_type` in the configuration:

- **Primal mode**: `solver_type = "gauss_newton"` with `gn_solver_mode = "primal"`
- **Dual mode**: `solver_type = "gauss_newton"` with `gn_solver_mode = "dual"` 
- **Least-squares mode**: `solver_type = "least_squares"`

All three modes solve the Gauss-Newton subproblem to compute a search direction at each optimization iteration, but use different formulations.

---

## 1. Primal Mode (Normal Equations)

**Formulation:** Solves the primal normal equations in design space:
```
(L^T L + λI) v = -G
```
where L is the Jacobian (tangent_dim × ndesign), λ is damping, and G is the gradient.

### Calling Sequence

```
solve() [optimproblem.cpp:1688]
├─ For each optimization iteration:
│  ├─ evalGradF(x, G)                    // Compute gradient G at current x
│  │
│  ├─ solveGaussNewtonNormalEqKSP(x, v_zero, -G, v_out)
│  │  ├─ setGaussNewtonEvalPoint(x)     // Cache evaluation point
│  │  ├─ MatCreateVecs(GN_NormalEq_MatShell, ...) // Create vectors (size ndesign)
│  │  ├─ KSPSolve(ksp_GN_NormalEq, -G, v_out)    // Solve with PETSc KSP
│  │  │  └─ [KSP internally calls matrix-vector products:]
│  │  │     └─ GN_NormalEq_MatMult(A, v, Av)      // Compute Av = (L^T L + λI) v
│  │  │        ├─ evalLinearizedForward(x, v)     // Compute L*v (forward linearization)
│  │  │        │  └─ timestepper->solveLinODE()   // Linearized Schrödinger equation
│  │  │        ├─ Project L*v to tangent space    // Get coordinates in Lie algebra
│  │  │        ├─ Reconstruct from tangent space  // Back to full state space
│  │  │        ├─ solveAdjointODE()               // Compute L^T (L*v)
│  │  │        └─ Add damping: Av += λ*v
│  │  └─ v_out now contains search direction
│  │
│  ├─ armijoLineSearch(x, f, G, v_out, x_new)    // Backtracking line search
│  └─ x ← x_new                                   // Update iterate
```

**Key Functions:**
- **`GN_NormalEq_MatMult`**: Matrix-vector product for `(L^T L + λI) v`
- **`evalLinearizedForward`**: Forward Jacobian application `L*v`
- **`solveAdjointODE`**: Adjoint Jacobian application `L^T*w`

**Matrix dimensions:**
- System matrix: `ndesign × ndesign` (e.g., 288 × 288)
- Vectors: size `ndesign`

---

## 2. Dual Mode (Tangent Space)

**Formulation:** Solves the dual formulation in tangent space, then recovers design space direction:
```
(L L^T + λI) y = b    where b = -∇_U J (projected to tangent space)
v = L^T y              (recover design space direction)
```
where the system is size `tangent_dim = dim² - 1` (e.g., 15 for 2×2 system).

### Calling Sequence

```
solve() [optimproblem.cpp:1688]
├─ For each optimization iteration:
│  ├─ evalGradF(x, G)                    // Compute gradient G
│  │
│  ├─ solveGaussNewtonDual(x, y_zero, v_out)  OR
│  │  solveGaussNewtonDualProjected(x, v_out)
│  │  │
│  │  ├─ setGaussNewtonEvalPoint(x)     // Cache evaluation point
│  │  ├─ evalNonlinearForward(x)        // Evaluate U(T) for current x
│  │  │
│  │  ├─ Compute RHS b = -∇_U J (projected to tangent space):
│  │  │  ├─ evalJ(U_final, ...)         // Evaluate objective J(U(T))
│  │  │  ├─ evalJ_diff(...)             // Compute ∇_U J via reverse-mode AD
│  │  │  └─ projectVecToTangentSpace()  // Project to Lie algebra (size tangent_dim)
│  │  │
│  │  ├─ GN_Dual_MatCreateVecs(GN_Dual_MatShell, ...) // Create vectors (size tangent_dim)
│  │  │
│  │  ├─ [For projected solver only:]
│  │  │  ├─ FormGaussNewtonDualMatrix()          // Build dual GN matrix explicitly
│  │  │  ├─ EPSSolve()                           // Compute top r eigenvalues/vectors
│  │  │  ├─ Direct solve in rank-r subspace      // Solve V_r^T A V_r * α = V_r^T b
│  │  │  └─ KSP refinement iterations            // Iteratively refine solution
│  │  │
│  │  ├─ [For standard dual solver:]
│  │  │  └─ KSPSolve(ksp_GN_Dual, b, y)          // Solve with PETSc KSP
│  │  │     └─ [KSP internally calls matrix-vector products:]
│  │  │        └─ GN_Dual_MatMult(A, y, Ay)       // Compute Ay = (L L^T + λI) y
│  │  │           ├─ GN_Jacobian_MatMultTranspose(A, y, v_temp)  // v_temp = L^T y
│  │  │           │  ├─ reconstructVecFromTangentSpace()  // Tangent space → full state
│  │  │           │  └─ solveAdjointODE()                // Apply adjoint Jacobian
│  │  │           ├─ GN_Jacobian_MatMult(A, v_temp, Ay)  // Ay = L * v_temp = L L^T y
│  │  │           │  ├─ evalLinearizedForward(x, v_temp) // Apply forward Jacobian
│  │  │           │  └─ projectVecToTangentSpace()       // Full state → tangent space
│  │  │           └─ Add damping: Ay += λ*y
│  │  │
│  │  └─ Recover design space direction:
│  │     └─ GN_Jacobian_MatMultTranspose(A, y, v_out)   // v_out = L^T y
│  │
│  ├─ Scale: v_out ← -v_out
│  ├─ armijoLineSearch(x, f, G, v_out, x_new)    // Backtracking line search
│  └─ x ← x_new                                   // Update iterate
```

**Key Functions:**
- **`GN_Dual_MatMult`**: Matrix-vector product for `(L L^T + λI) y`
- **`GN_Jacobian_MatMult`**: Pure Jacobian forward `y = L*v` (no augmentation)
- **`GN_Jacobian_MatMultTranspose`**: Pure Jacobian transpose `v = L^T*y` (no augmentation)
- **`projectVecToTangentSpace`**: Project state to Lie algebra coordinates
- **`reconstructVecFromTangentSpace`**: Reconstruct state from Lie algebra

**Matrix dimensions:**
- Dual system matrix: `tangent_dim × tangent_dim` (e.g., 15 × 15)
- Jacobian L: `tangent_dim × ndesign` (e.g., 15 × 288)
- RHS b: size `tangent_dim`
- Solution y: size `tangent_dim`
- Recovered v: size `ndesign`

---

## 3. Least-Squares Mode

**Formulation:** Solves the linear least-squares problem directly:
```
min_v ||L*v - b||²
```
Optionally augmented with Tikhonov regularization when `leastsquares_damping > 0`:
```
min_v ||[L; λI]*v - [b; 0]||²
```

### Calling Sequence

```
solve() [optimproblem.cpp:1820]
├─ For each optimization iteration:
│  ├─ evalGradF(x, G)                    // Compute gradient G
│  │
│  ├─ solveGaussNewtonLeastSquares(x, v_zero, v_out)
│  │  │
│  │  ├─ setGaussNewtonEvalPoint(x)     // Cache evaluation point
│  │  ├─ evalNonlinearForward(x)        // Evaluate U(T) for current x
│  │  │
│  │  ├─ Compute RHS b = -∇_U J (projected to tangent space):
│  │  │  ├─ GN_LeastSquares_MatCreateVecs(GN_LeastSquares_MatShell, ..., &b)
│  │  │  │  // Size = tangent_dim (no damping) or tangent_dim + ndesign (with damping)
│  │  │  ├─ VecZeroEntries(b)           // Initialize (bottom part stays zero if augmented)
│  │  │  ├─ evalJ(U_final, ...)         // Evaluate objective J(U(T))
│  │  │  ├─ evalJ_diff(...)             // Compute ∇_U J via reverse-mode AD
│  │  │  └─ projectVecToTangentSpace()  // Project to Lie algebra (fills top tangent_dim entries)
│  │  │
│  │  ├─ Solve least-squares system:
│  │  │  │
│  │  │  ├─ [If leastsquares_solver == "lsqr":]
│  │  │  │  └─ KSPSolve(ksp_LeastSquares, b, v_out)  // LSQR algorithm
│  │  │  │     └─ [LSQR internally calls matrix operations:]
│  │  │  │        ├─ GN_LeastSquares_MatMult(A, v, y)    // y = [L; λI]*v
│  │  │  │        │  ├─ evalLinearizedForward(x, v)      // Compute L*v (forward Jacobian)
│  │  │  │        │  ├─ projectVecToTangentSpace()       // Full state → tangent (top part)
│  │  │  │        │  └─ [If damping > 0:] y_bottom = λ*v // Append damping rows
│  │  │  │        │
│  │  │  │        └─ GN_LeastSquares_MatMultTranspose(A, w, v)  // v = [L; λI]^T*w
│  │  │  │           ├─ reconstructVecFromTangentSpace()    // Tangent → full state (top part)
│  │  │  │           ├─ solveAdjointODE()                  // Apply adjoint L^T
│  │  │  │           └─ [If damping > 0:] v += λ*w_bottom  // Add damping contribution
│  │  │  │
│  │  │  └─ [If leastsquares_solver == "brgn":]
│  │  │     └─ TaoSolve(tao_brgn, v_out)  // TAO Bounded Regularized Gauss-Newton
│  │  │        ├─ TaoBRGN_EvalResidual()  // F(v) = L*v - b
│  │  │        │  └─ GN_LeastSquares_MatMult(A, v, F)
│  │  │        │     └─ VecAXPY(F, -1.0, b)  // F -= b
│  │  │        │
│  │  │        └─ TaoBRGN_EvalJacobianResidual()  // Jacobian = L (constant)
│  │  │           └─ Returns GN_LeastSquares_MatShell (matrix-free)
│  │  │
│  │  └─ v_out now contains search direction
│  │
│  ├─ Scale: v_out ← -v_out
│  ├─ armijoLineSearch(x, f, G, v_out, x_new)    // Backtracking line search
│  └─ x ← x_new                                   // Update iterate
```

**Key Functions:**
- **`GN_LeastSquares_MatMult`**: Forward operator `y = [L; λI]*v` (augmented if damping > 0)
- **`GN_LeastSquares_MatMultTranspose`**: Adjoint operator `v = [L^T, λI]*w` (augmented if damping > 0)
- **`evalLinearizedForward`**: Forward Jacobian application
- **`solveAdjointODE`**: Adjoint Jacobian application

**Matrix dimensions:**
- **Without damping** (`leastsquares_damping = 0`):
  - Matrix L: `tangent_dim × ndesign` (e.g., 15 × 288)
  - RHS b: size `tangent_dim`
  
- **With damping** (`leastsquares_damping > 0`):
  - Augmented matrix [L; λI]: `(tangent_dim + ndesign) × ndesign` (e.g., 303 × 288)
  - Augmented RHS [b; 0]: size `tangent_dim + ndesign`

---

## Key Differences Summary

| Aspect | Primal Mode | Dual Mode | Least-Squares Mode |
|--------|-------------|-----------|-------------------|
| **System solved** | `(L^T L + λI) v = -G` | `(L L^T + λI) y = b`<br/>then `v = L^T y` | `min` ‖`L*v - b`‖² <br/>(optionally augmented) |
| **System size** | `ndesign × ndesign` | `tangent_dim × tangent_dim` | `tangent_dim × ndesign` (or augmented) |
| **System matrix** | GN_NormalEq_MatShell | GN_Dual_MatShell | GN_LeastSquares_MatShell |
| **Solver** | KSP (CG, MINRES, etc.) | KSP or Projected eigenspace | LSQR or TAO BRGN |
| **MatMult function** | `GN_NormalEq_MatMult` | `GN_Dual_MatMult` | `GN_LeastSquares_MatMult` |
| **Advantage** | Direct in design space | Smaller system for small tangent space | Natural for ill-conditioned problems |
| **Best for** | Well-conditioned problems | `tangent_dim << ndesign` | General purpose, augmented regularization |

---

## Common Operations

All three modes share these core operations:

- **`evalNonlinearForward(x)`**: Solve Schrödinger equation forward in time to get U(T)
- **`evalLinearizedForward(x, v)`**: Solve linearized Schrödinger equation to compute L*v
- **`solveAdjointODE(...)`**: Solve adjoint Schrödinger equation to compute L^T*w
- **`projectVecToTangentSpace(...)`**: Project state U to Lie algebra coordinates (SU(N) tangent space)
- **`reconstructVecFromTangentSpace(...)`**: Reconstruct state from Lie algebra coordinates
- **`armijoLineSearch(...)`**: Backtracking line search with Armijo condition

The key distinction is **how the Gauss-Newton subproblem is formulated and solved** - primal vs dual formulation, or direct least-squares.
