#include "PostProcessor.h"
#include <iostream>
#include <vector>
#include <cmath>

// Function to compute the stress field for all elements
std::vector<std::vector<double>> PostProcessor::computeStressField(const std::vector<double>& u, const std::vector<Triangle>& elements, const std::vector<Point>& nodes, double E, double nu) {
    std::vector<std::vector<double>> stress_field;

    // Loop over each element in the mesh
    for (const auto& element : elements) {
        // Get the displacements for the nodes of the element
        std::vector<double> u_element = {
            u[2 * element.a.x], u[2 * element.a.y],
            u[2 * element.b.x], u[2 * element.b.y],
            u[2 * element.c.x], u[2 * element.c.y]
        };

        // Compute the stress for this element
        std::vector<double> stress = computeElementStress(u_element, element, nodes, E, nu);

        // Store the stress in the stress field
        stress_field.push_back(stress);
    }

    return stress_field;  // Return the stress field for all elements
}

// Helper function to compute the stress for a single element
std::vector<double> PostProcessor::computeElementStress(const std::vector<double>& u_element, const Triangle& element, const std::vector<Point>& nodes, double E, double nu) {
    // Elasticity matrix D for plane stress (can be different for plane strain)
    double factor = E / (1 - nu * nu);
    std::vector<std::vector<double>> D = {
        { factor, factor * nu, 0 },
        { factor * nu, factor, 0 },
        { 0, 0, factor * (1 - nu) / 2 }
    };

    // Compute strain-displacement matrix B (simplified, actual FEM would compute based on geometry)
    // B is derived from the geometry of the triangle and relates the element displacements to strain
    std::vector<std::vector<double>> B = {
        { 1.0, 0.0, -1.0, 0.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0, 0.0, 1.0, -1.0 },
        { 0.0, 0.0, 1.0, -1.0, 0.0, 0.0 }
    };

    // Compute the strain for the element (epsilon = B * u_element)
    std::vector<double> strain(3, 0.0);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 6; ++j) {
            strain[i] += B[i][j] * u_element[j];
        }
    }

    // Compute the stress for the element (sigma = D * strain)
    std::vector<double> stress(3, 0.0);
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            stress[i] += D[i][j] * strain[j];
        }
    }

    return stress;  // Return the stress vector for this element
}
