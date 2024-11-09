#ifndef DELAUNAYTRIANGULATION_H
#define DELAUNAYTRIANGULATION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>

class Point {
public:
    float x, y;

    Point(float x, float y);
    std::vector<float> pos() const;
    bool isEqual(const Point& other) const;
    std::string pointToStr() const;
};

class Edge {
public:
    Point a, b;

    Edge(const Point& a, const Point& b);
    bool isEqual(const Edge& other) const;
    float length() const;
    bool edgeIntersection(const Edge& other) const;
};

class Triangle {
public:
    Point a, b, c;

    Triangle(const Point& a, const Point& b, const Point& c);
    bool isEqual(const Triangle& other) const;
};

class Graph {
private:
    std::vector<float> circumcircle(const Triangle& tri) const;
    bool pointInCircle(const Point& point, const std::vector<float>& circle) const;

public:
    std::vector<Point> points;
    std::vector<Triangle> triangles;
    std::vector<Edge> edges;

    int getNodeIndex(const Point& point) const;
    void addPoint(const Point& point);
    bool triangleIsDelaunay(const Triangle& triangle) const;
    void printTriangles() const;
    bool areCollinear(const Point& a, const Point& b, const Point& c);
};

#endif // DELAUNAYTRIANGULATION_H
