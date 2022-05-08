#include "graph.hpp"
#include <random>
#include <map>
#include <utility>
#include <cassert>
#include <algorithm>
#include <set>
#include <iostream>

int SEED = 0;

// ----- WeightedGraph -----

WeightedGraph::WeightedEdge::WeightedEdge(int to, long long weight):to(to), weight(weight) {}

WeightedGraph::WeightedTriangleEdge::WeightedTriangleEdge(int vertex1, int vertex2, long long weight):vertex1(vertex1), vertex2(vertex2), weight(weight) {}

WeightedGraph::WeightedGraph(int n_vertices):n_vertices(n_vertices) {
    edges.resize(n_vertices);
    triangles.resize(n_vertices);
    vertices_connected_by_triangle.resize(n_vertices);
    return;
}

void WeightedGraph::add_edge(int u, int v, long long weight) {
    assert(0 <= u && u < n_vertices && 0 <= v && v < n_vertices);
    edges[u].emplace_back(v, weight);
    edges[v].emplace_back(u, weight);
    return;
}

void WeightedGraph::add_triangle(int a, int b, int c, long long weight) {
    assert(0 <= a && a < n_vertices && 0 <= b && b < n_vertices && 0 <= c && c < n_vertices);

    triangles[a].emplace_back(b, c, weight);
    triangles[b].emplace_back(a, c, weight);
    triangles[c].emplace_back(a, b, weight);

    return;
}

// O(\Delta \ln(n) / \lambda^* \eps^2)
double WeightedGraph::greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps) {
    double delta = 0;
    for (int i = 0;i < n_vertices; i++) {
        double delta_i = 0;
        for (const auto &edge : edges[i]) {
            delta_i += edge.weight;
        }
        delta = std::max(delta, delta_i);
    }

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "ln(n) : " << log(n_vertices) << std::endl;
    std::cout << "lamda_star : " << lamda_star_lower_bound << std::endl;

    return (delta * log(n_vertices)) / (lamda_star_lower_bound * eps * eps);
}

double WeightedGraph::triangle_greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps) {
    double delta = 0;
    for (int i = 0;i < n_vertices; i++) {
        double delta_i = 0;
        for (const auto &triangle : triangles[i]) {
            delta_i += triangle.weight;
        }
        delta = std::max(delta, delta_i);
    }

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "ln(n) : " << log(n_vertices) << std::endl;
    std::cout << "lamda_star : " << lamda_star_lower_bound << std::endl;

    return (delta * log(n_vertices)) / (lamda_star_lower_bound * eps * eps);
}

// extract triangle O(n^3) or O(m^{3/2})
void WeightedGraph::extract_triangle(void) {

    std::map<std::pair<int, int>, long long> edge_to_weight;
    for (int i = 0;i < n_vertices; i++) {
        for (auto edge : edges[i]) {
            auto edge_normalized = std::minmax(i, edge.to);
            edge_to_weight[edge_normalized] = edge.weight;
        }
    }

    std::vector<std::set<int>> vertexset_connected_by_triangle(n_vertices);
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &edge1 : edges[i]) {
            for (const auto &edge2 : edges[i]) {

                if (edge1.to == edge2.to) break;
                if (i < edge1.to && i < edge2.to) {
                    auto other_edge = std::minmax(edge1.to, edge2.to);
                    if (edge_to_weight.count(other_edge) > 0) {
                        long long weight = edge1.weight + edge2.weight + edge_to_weight[other_edge];
                        add_triangle(i, edge1.to, edge2.to, weight);

                        vertexset_connected_by_triangle[i].insert(edge1.to);
                        vertexset_connected_by_triangle[i].insert(edge2.to);

                        vertexset_connected_by_triangle[edge1.to].insert(i);
                        vertexset_connected_by_triangle[edge1.to].insert(edge2.to);

                        vertexset_connected_by_triangle[edge2.to].insert(i);
                        vertexset_connected_by_triangle[edge2.to].insert(edge1.to);
                    }
                }

            }
        }
    }

    for (int i = 0;i < n_vertices; i++) {
        vertices_connected_by_triangle.reserve(vertexset_connected_by_triangle[i].size());
        for (auto v : vertexset_connected_by_triangle[i]) {
            vertices_connected_by_triangle[i].push_back(v);
        }
    }


    return;
}


// ----- Graph -----

Graph::Edge::Edge(int to):to(to) {}
Graph::TriangleEdge::TriangleEdge(int vertex1, int vertex2):vertex1(vertex1), vertex2(vertex2) {}

Graph::Graph(int n_vertices):n_vertices(n_vertices) {
    edges.resize(n_vertices);
    triangles.resize(n_vertices);
    vertices_connected_by_triangle.resize(n_vertices);
    return;
}

void Graph::add_edge(int u, int v) {
    assert(0 <= u && u < n_vertices && 0 <= v && v < n_vertices);
    edges[u].emplace_back(v);
    edges[v].emplace_back(u);
    return;
}

void Graph::add_triangle(int a, int b, int c) {
    assert(0 <= a && a < n_vertices && 0 <= b && b < n_vertices && 0 <= c && c < n_vertices);

    triangles[a].emplace_back(b, c);
    triangles[b].emplace_back(a, c);
    triangles[c].emplace_back(a, b);

    return;
}

// O(\Delta \ln(n) / \lambda^* \eps^2)
double Graph::greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps) {
    int delta = 0;
    for (int i = 0;i < n_vertices; i++) {
        delta = std::max(delta, (int) edges[i].size());
    }

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "ln(n) : " << log(n_vertices) << std::endl;
    std::cout << "lamda_star : " << lamda_star_lower_bound << std::endl;

    return (delta * log(n_vertices)) / (lamda_star_lower_bound * eps * eps);
}

double Graph::triangle_greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps) {
    int delta = 0;
    for (int i = 0;i < n_vertices; i++) {
        delta = std::max(delta, (int) triangles[i].size());
    }

    std::cout << "Delta : " << delta << std::endl;
    std::cout << "ln(n) : " << log(n_vertices) << std::endl;
    std::cout << "lamda_star : " << lamda_star_lower_bound << std::endl;

    return (delta * log(n_vertices)) / (lamda_star_lower_bound * eps * eps);
}

WeightedGraph Graph::put_uniform_weight(void) {
    // temporary weighted distribtion
    // think more about better or good weighting distribution
    std::mt19937 engine(SEED);

    std::uniform_int_distribution<> dist(0.0, 10000.0);

    WeightedGraph weighted_graph(n_vertices);

    for (int i = 0;i < n_vertices; i++) {
        for (const auto& edge : edges[i]) {

            if (i < edge.to) {
                long long weight = dist(engine);
                weighted_graph.add_edge(i, edge.to, weight);
            }

        }
    }

    return weighted_graph;
}

// extract triangle O(n^3) or O(m^{3/2})
void Graph::extract_triangle(void) {

    std::set<std::pair<int, int>> edges_set;
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &edge : edges[i]) {
            if (i < edge.to) {
                edges_set.emplace(i, edge.to);
            }
        }
    }

    std::vector<std::set<int>> vertexset_connected_by_triangle(n_vertices);
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &edge1 : edges[i]) {
            for (const auto &edge2 : edges[i]) {

                if (edge1.to == edge2.to) break;
                if (i < edge1.to && i < edge2.to) {
                    auto other_edge = std::minmax(edge1.to, edge2.to);
                    if (edges_set.count(other_edge) > 0) {
                        add_triangle(i, edge1.to, edge2.to);

                        vertexset_connected_by_triangle[i].insert(edge1.to);
                        vertexset_connected_by_triangle[i].insert(edge2.to);

                        vertexset_connected_by_triangle[edge1.to].insert(i);
                        vertexset_connected_by_triangle[edge1.to].insert(edge2.to);

                        vertexset_connected_by_triangle[edge2.to].insert(i);
                        vertexset_connected_by_triangle[edge2.to].insert(edge1.to);
                    }
                }

            }
        }
    }

    for (int i = 0;i < n_vertices; i++) {
        vertices_connected_by_triangle.reserve(vertexset_connected_by_triangle[i].size());
        for (auto v : vertexset_connected_by_triangle[i]) {
            vertices_connected_by_triangle[i].push_back(v);
        }
    }


    return;
}


// Data

Data::Data(double density, double triangle_density, double edge_density, double vertex_ratio)
    :density(density), triangle_density(triangle_density), edge_density(edge_density), vertex_ratio(vertex_ratio) {}

void Data::display(void) {
    std::cout << "degree density : " << 2 * density << std::endl;
    std::cout << "triangle density : " << 3 * triangle_density << std::endl;
    std::cout << "edge density : " << edge_density << std::endl;
    std::cout << "vertex ratio : " << vertex_ratio << std::endl;
    return;
}
