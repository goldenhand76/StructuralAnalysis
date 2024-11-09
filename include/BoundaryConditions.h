#pragma once
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>

// Structure to represent a load applied at a specific node
struct Load {
    int node;      // The node index where the load is applied
    double value;  // The value of the load
};

class BoundaryConditions {
public:
    // Function to apply boundary conditions
    void applyBoundaryConditions(std::vector<std::vector<double>>& K_global, std::vector<double>& F, const std::vector<int>& fixedNodes, const std::vector<Load>& loads);
};

#endif // BOUNDARYCONDITIONS_H
