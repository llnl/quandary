/**
 * @file test_gellmann_orthonormality.cpp
 * @brief Standalone test for orthonormality of Gell-Mann basis matrices
 *
 * Tests that the generalized Gell-Mann matrices form an orthonormal basis:
 *   Tr(T_a^† * T_b) = δ_ab
 *
 * Gell-Mann basis structure for dimension N:
 * - Symmetric:     S_jk = (1/√2)(|j⟩⟨k| + |k⟩⟨j|),  j < k  [N(N-1)/2 matrices]
 * - Antisymmetric: A_jk = (i/√2)(|k⟩⟨j| - |j⟩⟨k|),  j < k  [N(N-1)/2 matrices]
 * - Diagonal:      D_l = (1/√(l(l+1))) diag(1,...,1,-l,0,...,0), l=1,...,N-1  [N-1 matrices]
 * Total: N²-1 Hermitian, traceless, orthonormal matrices
 *
 * Compile:
 *   g++ -std=c++14 -o test_gellmann test_gellmann_orthonormality.cpp
 *
 * Run:
 *   ./test_gellmann
 */

#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cassert>

using Complex = std::complex<double>;
using Matrix = std::vector<std::vector<Complex>>;

// Helper: Create zero matrix
Matrix zeros(int N) {
    return Matrix(N, std::vector<Complex>(N, 0.0));
}

// Helper: Compute trace of a matrix
Complex trace(const Matrix& M) {
    Complex tr = 0.0;
    for (size_t i = 0; i < M.size(); ++i) {
        tr += M[i][i];
    }
    return tr;
}

// Helper: Compute Hermitian conjugate transpose
Matrix hermitian_conjugate(const Matrix& M) {
    int N = M.size();
    Matrix Mdag = zeros(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Mdag[i][j] = std::conj(M[j][i]);
        }
    }
    return Mdag;
}

// Helper: Matrix multiplication
Matrix matmul(const Matrix& A, const Matrix& B) {
    int N = A.size();
    Matrix C = zeros(N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Helper: Compute inner product Tr(A^† * B)
Complex inner_product(const Matrix& A, const Matrix& B) {
    Matrix Adag = hermitian_conjugate(A);
    Matrix prod = matmul(Adag, B);
    return trace(prod);
}

// Generate Symmetric Gell-Mann matrix S_jk (1-based indexing: j < k)
Matrix create_S(int N, int j, int k) {
    assert(1 <= j && j < k && k <= N);
    Matrix S = zeros(N);
    double norm = 1.0 / std::sqrt(2.0);
    S[j-1][k-1] = norm;  // |j⟩⟨k|
    S[k-1][j-1] = norm;  // |k⟩⟨j|
    return S;
}

// Generate Antisymmetric Gell-Mann matrix A_jk (1-based indexing: j < k)
Matrix create_A(int N, int j, int k) {
    assert(1 <= j && j < k && k <= N);
    Matrix A = zeros(N);
    double norm = 1.0 / std::sqrt(2.0);
    Complex i_unit(0.0, 1.0);
    A[k-1][j-1] = i_unit * norm;   // i|k⟩⟨j|
    A[j-1][k-1] = -i_unit * norm;  // -i|j⟩⟨k|
    return A;
}

// Generate Diagonal Gell-Mann matrix D_l (1-based indexing: l = 1,...,N-1)
Matrix create_D(int N, int l) {
    assert(1 <= l && l <= N-1);
    Matrix D = zeros(N);
    double norm = 1.0 / std::sqrt((double)l * (double)(l+1));

    // First l diagonal elements = 1
    for (int i = 0; i < l; ++i) {
        D[i][i] = norm;
    }
    // (l+1)-th diagonal element = -l
    D[l][l] = -l * norm;
    // Remaining elements = 0 (already initialized)

    return D;
}

// Generate all Gell-Mann matrices for dimension N
std::vector<Matrix> generate_gellmann_basis(int N) {
    std::vector<Matrix> basis;

    // Symmetric matrices: S_jk for j < k
    for (int j = 1; j <= N; ++j) {
        for (int k = j+1; k <= N; ++k) {
            basis.push_back(create_S(N, j, k));
        }
    }

    // Antisymmetric matrices: A_jk for j < k
    for (int j = 1; j <= N; ++j) {
        for (int k = j+1; k <= N; ++k) {
            basis.push_back(create_A(N, j, k));
        }
    }

    // Diagonal matrices: D_l for l = 1,...,N-1
    for (int l = 1; l <= N-1; ++l) {
        basis.push_back(create_D(N, l));
    }

    return basis;
}

// Print a matrix (for debugging)
void print_matrix(const Matrix& M, const std::string& name = "") {
    if (!name.empty()) {
        std::cout << name << ":" << std::endl;
    }
    for (const auto& row : M) {
        for (const auto& elem : row) {
            std::cout << std::setw(12) << elem << " ";
        }
        std::cout << std::endl;
    }
}

// Test orthonormality: Tr(T_a^† * T_b) = δ_ab
bool test_orthonormality(int N, double tol = 1e-12) {
    std::vector<Matrix> basis = generate_gellmann_basis(N);
    int num_basis = basis.size();

    std::cout << "Testing Gell-Mann basis for N = " << N << std::endl;
    std::cout << "Number of basis elements: " << num_basis
              << " (expected: " << N*N - 1 << ")" << std::endl;

    if (num_basis != N*N - 1) {
        std::cerr << "ERROR: Wrong number of basis elements!" << std::endl;
        return false;
    }

    // Print basis matrices for N=2, N=3, and N=4
    if (N == 2 || N == 3 || N == 4) {
        std::cout << "\nBasis matrices:" << std::endl;
        for (int a = 0; a < num_basis; ++a) {
            std::cout << "\nT_" << a << ":" << std::endl;
            print_matrix(basis[a]);
        }
        std::cout << std::endl;
    }

    bool all_passed = true;
    double max_diag_error = 0.0;
    double max_offdiag_error = 0.0;

    for (int a = 0; a < num_basis; ++a) {
        for (int b = 0; b < num_basis; ++b) {
            Complex ip = inner_product(basis[a], basis[b]);
            double expected = (a == b) ? 1.0 : 0.0;
            double error = std::abs(ip - expected);

            if (a == b) {
                max_diag_error = std::max(max_diag_error, error);
                if (error > tol) {
                    std::cout << "FAIL: Tr(T_" << a << "^† * T_" << a << ") = "
                              << ip << " (expected 1.0, error = " << error << ")" << std::endl;
                    all_passed = false;
                }
            } else {
                max_offdiag_error = std::max(max_offdiag_error, error);
                if (error > tol) {
                    std::cout << "FAIL: Tr(T_" << a << "^† * T_" << b << ") = "
                              << ip << " (expected 0.0, error = " << error << ")" << std::endl;
                    all_passed = false;
                }
            }
        }
    }

    std::cout << "\nOrthonormality test" << std::endl;
    std::cout << std::scientific << std::setprecision(3);
    std::cout << "Max diagonal error:     " << max_diag_error << std::endl;
    std::cout << "Max off-diagonal error: " << max_offdiag_error << std::endl;

    return all_passed;
}

// Test that all matrices are Hermitian
bool test_hermiticity(int N, double tol = 1e-12) {
    std::vector<Matrix> basis = generate_gellmann_basis(N);

    std::cout << "\nTesting Hermiticity (T_a^† = T_a)..." << std::endl;

    bool all_passed = true;
    double max_error = 0.0;

    for (size_t a = 0; a < basis.size(); ++a) {
        Matrix T = basis[a];
        Matrix Tdag = hermitian_conjugate(T);

        // Check T^† = T
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                double error = std::abs(T[i][j] - Tdag[i][j]);
                max_error = std::max(max_error, error);
                if (error > tol) {
                    std::cout << "FAIL: T_" << a << " is not Hermitian at ("
                              << i << "," << j << "), error = " << error << std::endl;
                    all_passed = false;
                }
            }
        }
    }

    std::cout << "Max Hermiticity error: " << max_error << std::endl;
    return all_passed;
}

// Test that all matrices are traceless
bool test_traceless(int N, double tol = 1e-12) {
    std::vector<Matrix> basis = generate_gellmann_basis(N);

    std::cout << "\nTesting Traceless property (Tr(T_a) = 0)..." << std::endl;

    bool all_passed = true;
    double max_error = 0.0;

    for (size_t a = 0; a < basis.size(); ++a) {
        Complex tr = trace(basis[a]);
        double error = std::abs(tr);
        max_error = std::max(max_error, error);

        if (error > tol) {
            std::cout << "FAIL: Tr(T_" << a << ") = " << tr
                      << " (expected 0.0, error = " << error << ")" << std::endl;
            all_passed = false;
        }
    }

    std::cout << "Max trace error: " << max_error << std::endl;
    return all_passed;
}

int main() {
    std::cout << "=== Gell-Mann Basis Orthonormality Test ===" << std::endl;
    std::cout << std::endl;

    bool all_tests_passed = true;

    // Test for different dimensions
    std::vector<int> test_dims = {2, 3, 4, 5};

    for (int N : test_dims) {
        std::cout << "========================================" << std::endl;

        bool ortho_pass = test_orthonormality(N);
        bool herm_pass = test_hermiticity(N);
        bool trace_pass = test_traceless(N);

        if (ortho_pass && herm_pass && trace_pass) {
            std::cout << "✓ ALL TESTS PASSED for N = " << N << std::endl;
        } else {
            std::cout << "✗ SOME TESTS FAILED for N = " << N << std::endl;
            all_tests_passed = false;
        }

        std::cout << std::endl;
    }

    std::cout << "========================================" << std::endl;
    if (all_tests_passed) {
        std::cout << "✓✓✓ ALL TESTS PASSED ✓✓✓" << std::endl;
        return 0;
    } else {
        std::cout << "✗✗✗ SOME TESTS FAILED ✗✗✗" << std::endl;
        return 1;
    }
}
