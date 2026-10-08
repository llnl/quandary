#pragma once
#include <complex>
#include <cmath>
#include <petscsys.h>
#include <petscmat.h>

/**
 * @brief Utilities for generalized Gell-Mann matrices and tangent space projection
 *
 * Tangent space T_{U(T)}SU(N) has real dimension N²-1, spanned by {τ_a = i*U(T)*T_a}
 * where {T_a} are generalized Gell-Mann matrices (Hermitian, traceless, orthonormal).
 *
 * For X in tangent space: X = sum_a h_a τ_a with real coefficients h_a = Im(Tr(T_a * U^† * X))
 *
 * Gell-Mann basis structure:
 * - Symmetric: S_jk = (1/√2)(e_j e_k^T + e_k e_j^T),  j < k
 * - Antisym:   A_jk = (i/√2)(e_k e_j^T - e_j e_k^T), j < k
 * - Diagonal:  D_l = (1/√(l(l+1))) diag(1,...,1,-l,0,...,0), l=1,...,N-1
 *
 * Total: C(N,2) + C(N,2) + (N-1) = N²-1
 */

namespace GellMann {

/**
 * @brief Project complex vector to tangent space at U(T)
 *
 * Computes tangent coefficients h_a from X_vec for column columnID
 *
 * @param X_vec Input vector ([real; imag] blocks)
 * @param U_final_re Real part of U(T)
 * @param U_final_im Imaginary part of U(T)
 * @param N Hilbert dimension
 * @param columnID Column index (0-based)
 * @param h_out Output coefficients (length N²-1)
 */
void projectVecToTangentSpace(Vec X_vec, Mat U_final_re, Mat U_final_im, PetscInt N, PetscInt columnID, PetscScalar *h_out);

/**
 * @brief Reconstruct vector from tangent coefficients (adjoint operation)
 *
 * Computes X_vec = sum_a h_a * τ_a * e_columnID where τ_a = i*U*T_a
 *
 * @param h_in Input coefficients (length N²-1)
 * @param U_final_re Real part of U(T)
 * @param U_final_im Imaginary part of U(T)
 * @param N Hilbert dimension
 * @param columnID Column index (0-based)
 * @param X_vec Output vector ([real; imag] blocks)
 */
void reconstructVecFromTangentSpace(const PetscScalar *h_in, Mat U_final_re, Mat U_final_im, PetscInt N, PetscInt columnID, Vec X_vec);



/* NEW VERSION THAT DOES NOT RELY ON X BEING HERMITIAN */
void projectMatToTangentSpace(Mat X_re, Mat X_im, PetscInt N, Vec h_out);

} // namespace GellMann
