#pragma once
#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

class Solver {
public:
    // Function to solve the linear system K_global * u = F
    std::vector<double> solveSystem(const std::vector<std::vector<double>>& K_global, const std::vector<double>& F);
};

#endif // SOLVER_H
