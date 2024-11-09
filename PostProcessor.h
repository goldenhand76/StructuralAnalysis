#pragma once
#ifndef POSTPROCESSOR_H
#define POSTPROCESSOR_H

#include <vector>
#include "DelaunayTriangulation.h"

class PostProcessor {
public:
    // Function to compute stress field
    std::vector<std::vector<double>> computeStressField(const std::vector<double>& u, const std::vector<Triangle>& elements, const std::vector<Point>& nodes, double E, double nu);

private:
    // Helper function to compute stress for an individual element
    std::vector<double> computeElementStress(const std::vector<double>& u_element, const Triangle& element, const std::vector<Point>& nodes, double E, double nu);
};

#endif // POSTPROCESSOR_H
