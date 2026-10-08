#include "gellmann.hpp"
#include <cmath>
#include <cassert>

namespace GellMann {

void projectVecToTangentSpace(Vec X_vec, Mat U_final_re, Mat U_final_im, PetscInt dim, PetscInt columnID, PetscScalar *h_out) {

  // Compute omega = U^† * X_vec
  Vec omega;
  VecDuplicate(X_vec, &omega);

  // Set up vector strides for real and imaginary parts
  IS isu, isv;
  PetscInt rank_petsc, mpisize_petsc;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  PetscInt localsize_u = dim / mpisize_petsc;
  PetscInt ilow = rank_petsc * localsize_u;
  ISCreateStride(PETSC_COMM_WORLD, localsize_u, ilow*2, 1, &isu);
  ISCreateStride(PETSC_COMM_WORLD, localsize_u, ilow*2+localsize_u, 1, &isv);

  Vec omega_re, omega_im;
  Vec x_re, x_im;
  VecGetSubVector(X_vec, isu, &x_re);
  VecGetSubVector(X_vec, isv, &x_im);
  VecGetSubVector(omega, isu, &omega_re);
  VecGetSubVector(omega, isv, &omega_im);

  // omega = U^† * X: omega_re = U_re^T * x_re + U_im^T * x_im
  //                  omega_im = U_re^T * x_im - U_im^T * x_re
  Vec temp;
  VecDuplicate(omega_re, &temp);
  MatMultTranspose(U_final_re, x_re, omega_re);
  MatMultTranspose(U_final_im, x_im, temp);
  VecAXPY(omega_re, 1.0, temp);
  MatMultTranspose(U_final_re, x_im, omega_im);
  MatMultTranspose(U_final_im, x_re, temp);
  VecAXPY(omega_im, -1.0, temp);
  VecDestroy(&temp);

  VecRestoreSubVector(omega, isu, &omega_re);
  VecRestoreSubVector(omega, isv, &omega_im);
  VecRestoreSubVector(X_vec, isu, &x_re);
  VecRestoreSubVector(X_vec, isv, &x_im);

  // Extract tangent coefficients h_a from omega
  PetscInt nSA_pairs = dim*(dim-1)/2;
  PetscInt offsetS = 0;
  PetscInt offsetA = nSA_pairs;
  PetscInt offsetD = 2*nSA_pairs;

  auto idxS = [&](PetscInt j, PetscInt k) -> PetscInt {
    PetscInt jj = j-1, kk = k-1;
    return offsetS + jj*dim - jj*(jj+1)/2 + (kk - jj - 1);
  };
  auto idxA = [&](PetscInt j, PetscInt k) -> PetscInt {
    PetscInt jj = j-1, kk = k-1;
    return offsetA + jj*dim - jj*(jj+1)/2 + (kk - jj - 1);
  };
  auto idxD = [&](PetscInt j) -> PetscInt {
    return offsetD + (j-1);
  };

  PetscInt j = columnID + 1;  // Convert to 1-based indexing

  // Off-diagonal: h_S = √2 * Im(omega_k), h_A = -√2 * Re(omega_k) for k > j
  for (PetscInt k = j+1; k <= dim; ++k) {
    PetscScalar Omega_kj_re = 0.0;
    PetscScalar Omega_kj_im = 0.0;

    PetscInt elem_idx = k - 1;
    if (ilow <= elem_idx && elem_idx < ilow + localsize_u) {
      PetscInt local_offset = elem_idx - ilow;
      PetscInt id_global_re = ilow * 2 + local_offset;
      PetscInt id_global_im = ilow * 2 + localsize_u + local_offset;
      VecGetValues(omega, 1, &id_global_re, &Omega_kj_re);
      VecGetValues(omega, 1, &id_global_im, &Omega_kj_im);
    }

    h_out[idxS(j,k)] += PetscSqrtReal(2.0) * Omega_kj_im;
    h_out[idxA(j,k)] += -PetscSqrtReal(2.0) * Omega_kj_re;
  }

  // Diagonal: extract Im(omega_j) and scatter to D_l coefficients
  PetscScalar Omega_jj_im = 0.0;
  PetscInt diag_idx = j - 1;

  if (ilow <= diag_idx && diag_idx < ilow + localsize_u) {
    PetscInt local_offset = diag_idx - ilow;
    PetscInt id_global_im = ilow * 2 + localsize_u + local_offset;
    VecGetValues(omega, 1, &id_global_im, &Omega_jj_im);
  }

  for (PetscInt l = j; l <= dim-1; ++l) {
    h_out[idxD(l)] += Omega_jj_im / PetscSqrtReal((PetscReal)l*(PetscReal)(l+1));
  }

  if (j >= 2) {
    PetscInt l = j - 1;
    h_out[idxD(l)] += -(PetscReal)l * Omega_jj_im / PetscSqrtReal((PetscReal)l*(PetscReal)(l+1));
  }

  ISDestroy(&isu);
  ISDestroy(&isv);
  VecDestroy(&omega);
}

void reconstructVecFromTangentSpace(const PetscScalar *h_in, Mat U_final_re, Mat U_final_im, PetscInt dim, PetscInt columnID, Vec X_vec) {

  // Adjoint of projectVecToTangentSpace
  // Reconstruct omega from h_in, then compute X = U * omega

  VecZeroEntries(X_vec);

  Vec omega;
  VecDuplicate(X_vec, &omega);
  VecZeroEntries(omega);

  IS isu, isv;
  PetscInt rank_petsc, mpisize_petsc;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank_petsc);
  MPI_Comm_size(PETSC_COMM_WORLD, &mpisize_petsc);
  PetscInt localsize_u = dim / mpisize_petsc;
  PetscInt ilow = rank_petsc * localsize_u;
  ISCreateStride(PETSC_COMM_WORLD, localsize_u, ilow*2, 1, &isu);
  ISCreateStride(PETSC_COMM_WORLD, localsize_u, ilow*2+localsize_u, 1, &isv);

  PetscInt nSA_pairs = dim*(dim-1)/2;
  PetscInt offsetS = 0;
  PetscInt offsetA = nSA_pairs;
  PetscInt offsetD = 2*nSA_pairs;

  auto idxS = [&](PetscInt j, PetscInt k) -> PetscInt {
    PetscInt jj = j-1, kk = k-1;
    return offsetS + jj*dim - jj*(jj+1)/2 + (kk - jj - 1);
  };
  auto idxA = [&](PetscInt j, PetscInt k) -> PetscInt {
    PetscInt jj = j-1, kk = k-1;
    return offsetA + jj*dim - jj*(jj+1)/2 + (kk - jj - 1);
  };
  auto idxD = [&](PetscInt j) -> PetscInt {
    return offsetD + (j-1);
  };

  PetscInt j = columnID + 1;  // Convert to 1-based indexing

  // Off-diagonal k > j: Im(omega_k) = h_S / √2, Re(omega_k) = -h_A / √2
  for (PetscInt k = j+1; k <= dim; ++k) {
    PetscInt elem_idx = k - 1;
    PetscScalar h_S = h_in[idxS(j,k)];
    PetscScalar h_A = h_in[idxA(j,k)];
    if (ilow <= elem_idx && elem_idx < ilow + localsize_u) {
      PetscInt local_offset = elem_idx - ilow;
      PetscInt id_global_im = ilow * 2 + localsize_u + local_offset;
      PetscInt id_global_re = ilow * 2 + local_offset;
      VecSetValue(omega, id_global_im, h_S / PetscSqrtReal(2.0), ADD_VALUES);
      VecSetValue(omega, id_global_re, -h_A / PetscSqrtReal(2.0), ADD_VALUES);
    }
  }

  // Off-diagonal k < j: contributions from pairs (k,j) via skew-Hermiticity
  for (PetscInt k = 1; k < j; ++k) {
    PetscScalar h_S = h_in[idxS(k, j)];
    PetscScalar h_A = h_in[idxA(k, j)];
    PetscInt elem_idx = k - 1;
    if (ilow <= elem_idx && elem_idx < ilow + localsize_u) {
      PetscInt local_offset = elem_idx - ilow;
      PetscInt id_global_re = ilow*2 + local_offset;
      PetscInt id_global_im = ilow*2 + localsize_u + local_offset;
      VecSetValue(omega, id_global_im, h_S / PetscSqrtReal(2.0), ADD_VALUES);
      VecSetValue(omega, id_global_re, h_A / PetscSqrtReal(2.0), ADD_VALUES);
    }
  }

  // Diagonal: gather contributions from D_l coefficients to Im(omega_j)
  PetscInt diag_idx = j - 1;
  if (ilow <= diag_idx && diag_idx < ilow + localsize_u) {
    PetscScalar omega_jj_im = 0.0;
    for (PetscInt l = j; l <= dim-1; ++l) {
      omega_jj_im += h_in[idxD(l)] / PetscSqrtReal((PetscReal)l*(PetscReal)(l+1));
    }
    if (j >= 2) {
      PetscInt l = j - 1;
      omega_jj_im += -((PetscReal)l) * h_in[idxD(l)] / PetscSqrtReal((PetscReal)l*(PetscReal)(l+1));
    }
    PetscInt local_offset = diag_idx - ilow;
    PetscInt id_global_im = ilow * 2 + localsize_u + local_offset;
    VecSetValue(omega, id_global_im, omega_jj_im, ADD_VALUES);
  }

  VecAssemblyBegin(omega);
  VecAssemblyEnd(omega);

  // X = U * omega
  Vec x_re, x_im;
  VecGetSubVector(X_vec, isu, &x_re);
  VecGetSubVector(X_vec, isv, &x_im);
  Vec omega_re, omega_im;
  VecGetSubVector(omega, isu, &omega_re);
  VecGetSubVector(omega, isv, &omega_im);
  Vec temp;
  VecDuplicate(x_re, &temp);

  MatMult(U_final_re, omega_re, x_re);
  MatMult(U_final_im, omega_im, temp);
  VecAXPY(x_re, -1.0, temp);
  MatMult(U_final_re, omega_im, x_im);
  MatMult(U_final_im, omega_re, temp);
  VecAXPY(x_im, 1.0, temp);
  VecDestroy(&temp);

  VecRestoreSubVector(X_vec, isu, &x_re);
  VecRestoreSubVector(X_vec, isv, &x_im);
  VecRestoreSubVector(omega, isu, &omega_re);
  VecRestoreSubVector(omega, isv, &omega_im);

  ISDestroy(&isu);
  ISDestroy(&isv);
  VecDestroy(&omega);
}



void projectMatToTangentSpace(Mat X_re, Mat X_im, PetscInt N, Vec h_out){

  // Set up h[a] = Im(Tr(G_a^† * X)) where G_a are the generalized Gell-Mann matrices, a=1,...,N²-1

  /*  1) symmetric part:     1/sqrt(2) * Im( (X_kj + X_jk) )       for j < k  */
  /*  2) antisymmetric part: 1/sqrt(2) * Im( i*(X_jk - X_kj) ) = 1/sqrt(2) * Re(X_jk - X_kj) for j < k  */
  /*  3) diagonal part:      1/sqrt(l(l+1)) * Im( \sum_{j=1}^{l} X_jj - l * X_(l+1)(l+1) ) for l = 1,...,N-1 */

  PetscInt nSA_pairs = N*(N-1)/2; // number of symmetric/antisymmetric pairs

  // Reset h_out
  VecZeroEntries(h_out);

  // Entry (row,col) of X (0-based);
  auto getEntry = [&](PetscInt row, PetscInt col, PetscScalar &re, PetscScalar &im) {
    re = 0.0; im = 0.0;
    MatGetValues(X_re, 1, &row, 1, &col, &re);
    MatGetValues(X_im, 1, &row, 1, &col, &im);
  };

  const PetscReal invsqrt2 = 1.0 / PetscSqrtReal(2.0);

  // Symmetric and antisymmetric parts, ordered as in idxS/idxA of projectVecToTangentSpace
  for (PetscInt j = 1; j <= N; ++j) {
    for (PetscInt k = j+1; k <= N; ++k) {
      PetscInt jj = j-1, kk = k-1;
      PetscInt pair = jj*N - jj*(jj+1)/2 + (kk - jj - 1);
      PetscScalar jk_re, jk_im, kj_re, kj_im;
      getEntry(jj, kk, jk_re, jk_im);
      getEntry(kk, jj, kj_re, kj_im);
      VecSetValue(h_out, pair, invsqrt2 * (kj_im + jk_im), INSERT_VALUES);
      VecSetValue(h_out, nSA_pairs + pair, invsqrt2 * (jk_re - kj_re), INSERT_VALUES);
    }
  }

  // Diagonal part: running sum of Im(X_jj)
  PetscScalar cumsum = 0.0;
  for (PetscInt l = 1; l <= N-1; ++l) {
    PetscScalar re, im_l, im_next;
    getEntry(l-1, l-1, re, im_l);
    getEntry(l, l, re, im_next);
    cumsum += im_l;
    VecSetValue(h_out, 2*nSA_pairs + (l-1), (cumsum - l * im_next) / PetscSqrtReal(l*(l+1)), INSERT_VALUES);
  }
  VecAssemblyBegin(h_out);
  VecAssemblyEnd(h_out);
}

} // namespace GellMann
