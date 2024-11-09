#include "Solver.h"
#include <iostream>
#include <vector>
#include <stdexcept>

// Gaussian elimination solver for the system K_global * u = F
std::vector<double> Solver::solveSystem(const std::vector<std::vector<double>>& K_global, const std::vector<double>& F) {
    int n = K_global.size();
    std::vector<std::vector<double>> A = K_global;  // Create a copy of K_global for manipulation
    std::vector<double> b = F;  // Create a copy of F (the right-hand side vector)

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Pivoting
        double maxEl = std::abs(A[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > maxEl) {
                maxEl = std::abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap the maximum row with the current row (pivoting)
        for (int k = i; k < n; ++k) {
            std::swap(A[maxRow][k], A[i][k]);
        }
        std::swap(b[maxRow], b[i]);

        // Make all rows below this one 0 in the current column
        for (int k = i + 1; k < n; ++k) {
            double c = -A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                if (i == j) {
                    A[k][j] = 0;
                }
                else {
                    A[k][j] += c * A[i][j];
                }
            }
            b[k] += c * b[i];
        }
    }

    // Solve equation A*u = b for an upper triangular matrix A
    std::vector<double> u(n);
    for (int i = n - 1; i >= 0; --i) {
        u[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= A[k][i] * u[i];
        }
    }

    return u;  // Return the solution vector (displacement/flow field)
}
