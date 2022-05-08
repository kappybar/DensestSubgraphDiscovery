#include "greedy.hpp"
#include <numeric>
#include <queue>
#include <vector>

#include <iostream>

// densest subgraph
// 1/2 approximation 
// unweighted
Data pq_greedy(const Graph &graph) {
    int n_vertices = graph.n_vertices;
    int n_edges = 0;

    auto compare = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        decltype(compare)> pq {compare};
    std::vector<bool> deleted(n_vertices, false);
    std::vector<int> degree(n_vertices, 0);

    for (int i = 0;i < graph.n_vertices; i++) {
        int degree_i = (int) graph.edges[i].size();
        n_edges += degree_i;
        degree[i] = degree_i;
        pq.emplace(degree_i, i);
    }

    n_edges /= 2;
    double maximum_degree_density = n_edges / (double) n_vertices;
    double vertex_ratio = 1;
    double edge_density = n_edges / ((double)n_vertices * (n_vertices - 1) / 2);

    // greedy peeling
    while (n_vertices > 1 && !pq.empty()) {
        auto [degree_i, i] = pq.top();
        pq.pop();
        if (deleted[i]) continue;
        deleted[i] = true;

        n_vertices --;

        for (const auto &edge : graph.edges[i]) {
            if (deleted[edge.to]) continue;
            degree[edge.to] --;
            n_edges --;
            pq.emplace(degree[edge.to], edge.to);
        }

        if (maximum_degree_density < n_edges / (double) n_vertices) {
            maximum_degree_density = n_edges / (double) n_vertices;
            vertex_ratio = n_vertices / (double) graph.n_vertices;
            edge_density = n_edges / ((double)n_vertices * (n_vertices - 1) / 2);
        }
        
    }

    Data result = Data(maximum_degree_density,
                       0,
                       edge_density,
                       vertex_ratio);

    return result;
}

// densest subgraph
// 1/2 approximation 
// weighted
double pq_weighted_greedy(const WeightedGraph &weighted_graph) {
    int n_vertices = weighted_graph.n_vertices;
    long long sum_edge_weight = 0;

    auto compare = [](const std::pair<long long, int> &a, const std::pair<long long, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<long long, int>, 
                        std::vector<std::pair<long long, int>>, 
                        decltype(compare)> pq {compare};
    std::vector<bool> deleted(n_vertices, false);
    std::vector<long long> weighted_degree(n_vertices, 0);

    for (int i = 0;i < weighted_graph.n_vertices; i++) {
        for (const auto &edge : weighted_graph.edges[i]) {
            weighted_degree[i] += edge.weight;
        }
        sum_edge_weight += weighted_degree[i];
        pq.emplace(weighted_degree[i], i);
    }
    
    sum_edge_weight /= 2;
    double maximum_weighted_edge_density = sum_edge_weight / (double) n_vertices;

    // greedy peeling
    while (n_vertices > 1 && !pq.empty()) {
        auto [weighted_degree_i, i] = pq.top();
        pq.pop();
        if (deleted[i]) continue;
        deleted[i] = true;

        n_vertices --;

        for (const auto &edge : weighted_graph.edges[i]) {
            if (deleted[edge.to]) continue;
            weighted_degree[edge.to] -= edge.weight;
            sum_edge_weight -= edge.weight;
            pq.emplace(weighted_degree[edge.to], edge.to);
        }

        maximum_weighted_edge_density = std::max(maximum_weighted_edge_density, sum_edge_weight / (double) n_vertices);
    }

    return maximum_weighted_edge_density;
}

// triangle densest subgraph
// 1/3 approximation
// unweighted
double pq_triangle_greedy(const Graph &graph) {
    int n_vertices = graph.n_vertices;
    int n_triangle = 0; 

    auto compare = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        decltype(compare)> pq {compare};
    std::vector<bool> deleted(n_vertices, false);
    std::vector<int> triangle(n_vertices, 0);

    for (int i = 0;i < graph.n_vertices; i++) {
        int triangle_i = (int) graph.triangles[i].size();
        n_triangle += triangle_i;
        triangle[i] = triangle_i;
        pq.emplace(triangle_i, i);
    }

    n_triangle /= 3;
    double maximum_triangle_density = n_triangle / (double) n_vertices;

    // greedy peeling
    while (n_vertices > 1 && !pq.empty()) {
        auto [triangle_i, i] = pq.top();
        pq.pop();
        if (deleted[i]) continue;
        deleted[i] = true;

        n_vertices --;

        for (const auto &triangle_edge : graph.triangles[i]) {
            if (deleted[triangle_edge.vertex1] || deleted[triangle_edge.vertex2]) continue;
            triangle[triangle_edge.vertex1] --;
            triangle[triangle_edge.vertex2] --;
            n_triangle --;
        }

        // improve?
        for (int v : graph.vertices_connected_by_triangle[i]) {
            if (deleted[v]) continue;
            pq.emplace(triangle[v], v);
        }

        maximum_triangle_density = std::max(maximum_triangle_density, n_triangle / (double) n_vertices);
    }

    return maximum_triangle_density;
}


// triangle densest triangle
// 1/3 approximation
// weighted
double pq_weighted_triangle_greedy(const WeightedGraph &weighted_graph) {
    int n_vertices = weighted_graph.n_vertices;
    long long sum_triangle_weight = 0; 

    auto compare = [](const std::pair<long long, int> &a, const std::pair<long long, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<long long, int>,
                        std::vector<std::pair<long long, int>>,
                        decltype(compare)> pq {compare};
    std::vector<bool> deleted(n_vertices, false);
    std::vector<long long> triangle_weight(n_vertices, 0);

    for (int i = 0;i < weighted_graph.n_vertices; i++) {
        for (const auto &triangle : weighted_graph.triangles[i]) {
            triangle_weight[i] += triangle.weight;
        }
        sum_triangle_weight += triangle_weight[i];
        pq.emplace(triangle_weight[i], i);
    }

    sum_triangle_weight /= 3;
    double maximum_weight_triangle_density = sum_triangle_weight / (double) n_vertices;

    // greedy peeling
    while (n_vertices > 1 && !pq.empty()) {
        auto [triangle_weight_i, i] = pq.top();
        pq.pop();
        if (deleted[i]) continue;
        deleted[i] = true;

        n_vertices --;

        for (const auto &triangle_edge : weighted_graph.triangles[i]) {
            if (deleted[triangle_edge.vertex1] || deleted[triangle_edge.vertex2]) continue;
            triangle_weight[triangle_edge.vertex1] -= triangle_edge.weight;
            triangle_weight[triangle_edge.vertex2] -= triangle_edge.weight;
            sum_triangle_weight -= triangle_edge.weight;
        }

        // improve?
        for (int v : weighted_graph.vertices_connected_by_triangle[i]) {
            if (deleted[v]) continue;
            pq.emplace(triangle_weight[v], v);
        }

        maximum_weight_triangle_density = std::max(maximum_weight_triangle_density, sum_triangle_weight / (double) n_vertices);
    }

    return maximum_weight_triangle_density;
}