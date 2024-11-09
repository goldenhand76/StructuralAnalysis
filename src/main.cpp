#include "DelaunayTriangulation.h"
#include "StiffnessMatrixAssembler.h"
#include "BoundaryConditions.h"
#include "Solver.h"
#include "PostProcessor.h"

int main() {
    Graph graph;

    for (int i = 0; i < 20; i++) {
        for (int j = 0; j < 20; j++) {
            graph.addPoint(Point(i, j));
        }
    }

    // Print triangles
    graph.printTriangles();

    StiffnessMatrixAssembler assembler;
    BoundaryConditions bc;
    Solver solver;
    PostProcessor postProcessor;


    // Define material properties
    double E = 210e9; // Young's modulus in Pascals
    double nu = 0.3;  // Poisson's ratio

    // Create a global stiffness matrix (initialize as an empty 2D vector)
    std::vector<std::vector<double>> K_global;
    std::vector<double> F(graph.points.size(), 0.0); // Initialize force vector with zeros

    // Assemble the global stiffness matrix
    assembler.assembleGlobalStiffness(graph, E, nu, K_global);

    // Print the global stiffness matrix
    std::cout << "\nGlobal Stiffness Matrix:" << std::endl;
    for (const auto& row : K_global) {
        for (const auto& value : row) {
            std::cout << value << "\t";
        }
        std::cout << std::endl;
    }

    // Define fixed nodes (for example: fix node 0 and node 1)
    std::vector<int> fixedNodes = { 0, 1 };

    // Define loads (for example: apply a load of 1000 N at node 2)
    std::vector<Load> loads = {
        { 2, 1000.0 } // Load of 1000 N applied at node 2
    };

    // Apply boundary conditions
    bc.applyBoundaryConditions(K_global, F, fixedNodes, loads);

    // Print the modified global stiffness matrix and force vector
    std::cout << "\nGlobal Stiffness Matrix after Boundary Conditions:" << std::endl;
    for (const auto& row : K_global) {
        for (const auto& value : row) {
            std::cout << value << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << "\nForce Vector:" << std::endl;
    for (const auto& value : F) {
        std::cout << value << "\t";
    }
    std::cout << std::endl;

    // Solve the system for displacements
    std::vector<double> displacements = solver.solveSystem(K_global, F);

    // Print the results
    std::cout << "\nDisplacement Vector:" << std::endl;
    for (const auto& value : displacements) {
        std::cout << value << "\t";
    }
    std::cout << std::endl;

    // Post-process to compute the stress field
    std::vector<std::vector<double>> stress_field = postProcessor.computeStressField(displacements, graph.triangles, graph.points, E, nu);

    // Print the stress results
    std::cout << "\nStress Field:" << std::endl;
    for (const auto& stress : stress_field) {
        std::cout << "Stress: ";
        for (const auto& s : stress) {
            std::cout << s << "\t";
        }
        std::cout << std::endl;
    }

    return 0; 
}

// 100 Points : 14 Sec
// 400 Points : more than 126 Min
// CPU Usage  : 8%