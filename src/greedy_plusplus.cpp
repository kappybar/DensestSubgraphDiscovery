#include "greedy_plusplus.hpp"
#include <numeric>
#include <queue>
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>

// densest subgraph
// GREEDY++
// unweighted
Data pq_greedy_plusplus(const Graph &graph, int T) {

    auto compare = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        decltype(compare)> pq {compare};
    
    std::vector<int> load(graph.n_vertices, 0);

    double maximum_degree_density = 0;
    double vertex_ratio = 0;
    double edge_density = 0;


    // iter T times
    for (int t = 0;t < T; t++) {
        if (t % 100 == 0) {
            std::cout << "iter : " << t << std::endl;
            std::cout << "temporal maximum density : " << maximum_degree_density << std::endl;
        }

        int n_vertices = graph.n_vertices;
        int n_edges = 0;

        std::vector<bool> deleted(n_vertices, false);
        std::vector<int> degree(n_vertices, 0);

        for (int i = 0;i < graph.n_vertices; i++) {
            int degree_i = (int) graph.edges[i].size();
            n_edges += degree_i;
            degree[i] = degree_i;
            pq.emplace(degree_i + load[i], i);
        }

        n_edges /= 2;
        if (maximum_degree_density < n_edges / (double) n_vertices) {
            maximum_degree_density = n_edges / (double) n_vertices;
            vertex_ratio = n_vertices / (double) graph.n_vertices;
            edge_density = n_edges / ((double) n_vertices * (n_vertices - 1) / 2);
        }

        // greedy peeling
        while (n_vertices > 1 && !pq.empty()) {
            auto [weight_i, i] = pq.top();
            pq.pop();
            if (deleted[i]) continue;
            deleted[i] = true;

            n_vertices --;

            for (const auto &edge : graph.edges[i]) {
                if (deleted[edge.to]) continue;
                degree[edge.to] --;
                n_edges --;
                pq.emplace(degree[edge.to] + load[edge.to], edge.to);

                load[i] ++;
            }

            if (maximum_degree_density < n_edges / (double) n_vertices) {
                maximum_degree_density = n_edges / (double) n_vertices;
                vertex_ratio = n_vertices / (double) graph.n_vertices;
                edge_density = n_edges / ((double) n_vertices * (n_vertices - 1) / 2);
            }
        }
    }

    Data result = Data(maximum_degree_density,
                       0,
                       edge_density,
                       vertex_ratio);

    return result;
}



// densest subgraph
// GREEDY++
// weighted
double pq_weighted_greedy_plusplus(const WeightedGraph &weighted_graph, int T) {

    auto compare = [](const std::pair<long long, int> &a, const std::pair<long long, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<long long, int>, 
                        std::vector<std::pair<long long, int>>, 
                        decltype(compare)> pq {compare};

    std::vector<long long> load(weighted_graph.n_vertices, 0);
    double maximum_weighted_edge_density = 0;

    // iter T times
    for (int t = 0;t < T; t++) {
        if (t % 100 == 0) {
            std::cout << "iter : " << t << std::endl;
            std::cout << "temporal maximum density : " << maximum_weighted_edge_density << std::endl;
        }
        int n_vertices = weighted_graph.n_vertices;
        long long sum_edge_weight = 0;
        
        std::vector<bool> deleted(n_vertices, false);
        std::vector<long long> weighted_degree(n_vertices, 0);

        for (int i = 0;i < weighted_graph.n_vertices; i++) {
            for (const auto &edge : weighted_graph.edges[i]) {
                weighted_degree[i] += edge.weight;
            }
            sum_edge_weight += weighted_degree[i];
            pq.emplace(weighted_degree[i] + load[i], i);
        }
        
        sum_edge_weight /= 2;
        maximum_weighted_edge_density = std::max(maximum_weighted_edge_density, sum_edge_weight / (double) n_vertices);

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
                pq.emplace(weighted_degree[edge.to] + load[edge.to], edge.to);

                load[i] += edge.weight;
            }

            maximum_weighted_edge_density = std::max(maximum_weighted_edge_density, sum_edge_weight / (double) n_vertices);
        }
    }

    return maximum_weighted_edge_density;
}


// triangle densest subgraph
// Super-Greedy++
// unweighted
double pq_triangle_greedy_plusplus(const Graph &graph, int T) {

    auto compare = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        decltype(compare)> pq {compare};

    std::vector<int> load(graph.n_vertices, 0);
    double maximum_triangle_density = 0;

    for (int t = 0;t < T; t++) {
        if (t % 10 == 0) {
            std::cout << "iter :" << t << std::endl;
            std::cout << "temporal max density : " << maximum_triangle_density << std::endl;
        }
        int n_vertices = graph.n_vertices;
        int n_triangle = 0; 
        std::vector<bool> deleted(n_vertices, false);
        std::vector<int> triangle(n_vertices, 0);

        for (int i = 0;i < graph.n_vertices; i++) {
            int triangle_i = (int) graph.triangles[i].size();
            n_triangle += triangle_i;
            triangle[i] = triangle_i;
            pq.emplace(triangle_i + load[i], i);
        }

        n_triangle /= 3;
        maximum_triangle_density = std::max(maximum_triangle_density, n_triangle / (double) n_vertices);

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

                load[i] ++;
            }

            // improve?
            for (int v : graph.vertices_connected_by_triangle[i]) {
                if (deleted[v]) continue;
                pq.emplace(triangle[v] + load[v], v);
            }

            maximum_triangle_density = std::max(maximum_triangle_density, n_triangle / (double) n_vertices);
        }

    }

    return maximum_triangle_density;
}



// triangle densest triangle
// Super-Greedy++
// weighted
double pq_weighted_triangle_greedy_plusplus(const WeightedGraph &weighted_graph, int T) {

    auto compare = [](const std::pair<long long, int> &a, const std::pair<long long, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<long long, int>,
                        std::vector<std::pair<long long, int>>,
                        decltype(compare)> pq {compare};

    std::vector<long long> load(weighted_graph.n_vertices, 0);
    double maximum_weight_triangle_density = 0;

    for (int t = 0;t < T; t++) {
        if (t % 10 == 0) {
            std::cout << "iter : " << t << std::endl;
            std::cout << "temporal maximum : " << maximum_weight_triangle_density << std::endl;
        }
        int n_vertices = weighted_graph.n_vertices;
        long long sum_triangle_weight = 0; 
        std::vector<bool> deleted(n_vertices, false);
        std::vector<long long> triangle_weight(n_vertices, 0);

        for (int i = 0;i < weighted_graph.n_vertices; i++) {
            for (const auto &triangle : weighted_graph.triangles[i]) {
                triangle_weight[i] += triangle.weight;
            }
            sum_triangle_weight += triangle_weight[i];
            pq.emplace(triangle_weight[i] + load[i], i);
        }

        sum_triangle_weight /= 3;
        maximum_weight_triangle_density = std::max(maximum_weight_triangle_density, sum_triangle_weight / (double) n_vertices);

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

                load[i] += triangle_edge.weight;
            }

            // improve?
            for (int v : weighted_graph.vertices_connected_by_triangle[i]) {
                if (deleted[v]) continue;
                pq.emplace(triangle_weight[v] + load[v], v);
            }

            maximum_weight_triangle_density = std::max(maximum_weight_triangle_density, sum_triangle_weight / (double) n_vertices);
        }
    }

    return maximum_weight_triangle_density;
}


// densest subgraph
// GREEDY++
// unweighted
Data my_greedy_plusplus(const Graph &graph, int T, int alpha) {

    auto compare = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        decltype(compare)> pq {compare};
    
    std::vector<int> load(graph.n_vertices, 0);

    double maximum_degree_density = 0;
    double vertex_ratio = 0;
    double edge_density = 0;


    // iter T times
    for (int t = 0;t < T; t++) {
        if (t % 100 == 0) {
            std::cout << "iter : " << t << std::endl;
            std::cout << "temporal maximum density : " << maximum_degree_density << std::endl;
        }

        int n_vertices = graph.n_vertices;
        int n_edges = 0;

        std::vector<bool> deleted(n_vertices, false);
        std::vector<int> degree(n_vertices, 0);

        for (int i = 0;i < graph.n_vertices; i++) {
            int degree_i = (int) graph.edges[i].size();
            n_edges += degree_i;
            degree[i] = degree_i;
            pq.emplace(degree_i + load[i], i);
        }

        n_edges /= 2;
        if (maximum_degree_density < n_edges / (double) n_vertices) {
            maximum_degree_density = n_edges / (double) n_vertices;
            vertex_ratio = n_vertices / (double) graph.n_vertices;
            edge_density = n_edges / ((double) n_vertices * (n_vertices - 1) / 2);
        }

        // greedy peeling
        while (n_vertices > 1 && !pq.empty()) {
            auto [weight_i, i] = pq.top();
            pq.pop();
            if (deleted[i]) continue;
            deleted[i] = true;

            n_vertices --;

            for (const auto &edge : graph.edges[i]) {
                if (deleted[edge.to]) continue;
                degree[edge.to] --;
                n_edges --;
                pq.emplace(degree[edge.to] + load[edge.to], edge.to);

                load[i] += alpha; //best
            }

            if (maximum_degree_density < n_edges / (double) n_vertices) {
                maximum_degree_density = n_edges / (double) n_vertices;
                vertex_ratio = n_vertices / (double) graph.n_vertices;
                edge_density = n_edges / ((double) n_vertices * (n_vertices - 1) / 2);
            }
        }
    }

    Data result = Data(maximum_degree_density,
                       0,
                       edge_density,
                       vertex_ratio);

    return result;
}



// triangle densest subgraph
// Super-Greedy++
// unweighted
double my_triangle_greedy_plusplus(const Graph &graph, int T, int alpha) {

    auto compare = [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<int, int>, 
                        std::vector<std::pair<int, int>>, 
                        decltype(compare)> pq {compare};

    std::vector<int> load(graph.n_vertices, 0);
    double maximum_triangle_density = 0;

    for (int t = 0;t < T; t++) {
        if (t % 10 == 0) {
            std::cout << "iter :" << t << std::endl;
            std::cout << "temporal max density : " << maximum_triangle_density << std::endl;
        }
        int n_vertices = graph.n_vertices;
        int n_triangle = 0; 
        std::vector<bool> deleted(n_vertices, false);
        std::vector<int> triangle(n_vertices, 0);

        for (int i = 0;i < graph.n_vertices; i++) {
            int triangle_i = (int) graph.triangles[i].size();
            n_triangle += triangle_i;
            triangle[i] = triangle_i;
            pq.emplace(triangle_i + load[i], i);
        }

        n_triangle /= 3;
        maximum_triangle_density = std::max(maximum_triangle_density, n_triangle / (double) n_vertices);

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

                load[i] += alpha; // faster
            }

            // improve?
            for (int v : graph.vertices_connected_by_triangle[i]) {
                if (deleted[v]) continue;
                pq.emplace(triangle[v] + load[v], v);
            }

            maximum_triangle_density = std::max(maximum_triangle_density, n_triangle / (double) n_vertices);
        }

    }

    return maximum_triangle_density;
}


// weighted_graph have only k-top wieghted triangle
// let S := subgraph that maximize k-top-triangle-density
// let a := k-top-triangle-densty(S)
// let b := triangle-density(S)
// return (a, b)
std::pair<double, double> ktop_triangle_greedy_plusplus(const WeightedGraph &weighted_graph, int T) {
    auto compare = [](const std::pair<long long, int> &a, const std::pair<long long, int> &b) {
        return a.first > b.first;
    };

    std::priority_queue<std::pair<long long, int>,
                        std::vector<std::pair<long long, int>>,
                        decltype(compare)> pq {compare};

    std::vector<long long> load(weighted_graph.n_vertices, 0);
    double maximum_weight_triangle_density = 0;
    std::vector<bool> S(weighted_graph.n_vertices, 0);

    for (int t = 0;t < T; t++) {
        if (t % 10 == 0) {
            std::cout << "iter : " << t << std::endl;
            std::cout << "temporal maximum : " << maximum_weight_triangle_density << std::endl;
        }
        int n_vertices = weighted_graph.n_vertices;
        long long sum_triangle_weight = 0; 
        std::vector<bool> deleted(n_vertices, false);
        std::vector<long long> triangle_weight(n_vertices, 0);

        for (int i = 0;i < weighted_graph.n_vertices; i++) {
            for (const auto &triangle : weighted_graph.triangles[i]) {
                triangle_weight[i] += triangle.weight;
            }
            sum_triangle_weight += triangle_weight[i];
            pq.emplace(triangle_weight[i] + load[i], i);
        }

        sum_triangle_weight /= 3;
        std::vector<bool> tmp_S(n_vertices, 1);
        if (maximum_weight_triangle_density < sum_triangle_weight / (double) n_vertices) {
            maximum_weight_triangle_density =  sum_triangle_weight / (double) n_vertices;
            S = tmp_S;
        }

        // greedy peeling
        while (n_vertices > 1 && !pq.empty()) {
            auto [triangle_weight_i, i] = pq.top();
            pq.pop();
            if (deleted[i]) continue;
            deleted[i] = true;
            tmp_S[i] = 0;

            n_vertices --;

            for (const auto &triangle_edge : weighted_graph.triangles[i]) {
                if (deleted[triangle_edge.vertex1] || deleted[triangle_edge.vertex2]) continue;
                triangle_weight[triangle_edge.vertex1] -= triangle_edge.weight;
                triangle_weight[triangle_edge.vertex2] -= triangle_edge.weight;
                sum_triangle_weight -= triangle_edge.weight;

                // load[i] += triangle_edge.weight;
                load[i] += triangle_edge.weight * 10;
            }

            // improve?
            for (int v : weighted_graph.vertices_connected_by_triangle[i]) {
                if (deleted[v]) continue;
                pq.emplace(triangle_weight[v] + load[v], v);
            }

            if (maximum_weight_triangle_density < sum_triangle_weight / (double) n_vertices) {
                maximum_weight_triangle_density = sum_triangle_weight / (double) n_vertices;
                S = tmp_S;
            }
        }
    }

    // I need this section because I have to evaluate actual weighted triangle density (not only k-top weighted triangle)

    std::map<std::pair<int, int>, long long> edge2weight;
    for (int i = 0;i < weighted_graph.n_vertices; i++) {
        if (S[i]) {
            for (const auto &edge : weighted_graph.edges[i]) {
                if (S[edge.to] && i < edge.to) {
                    edge2weight[std::make_pair(i, edge.to)] = edge.weight;
                }
            }
        }
    }

    int n_vertices_of_S = 0;
    long long sum_weight_of_S = 0;
    for (int i = 0;i < weighted_graph.n_vertices; i++) {
        if (S[i]) {
            n_vertices_of_S ++;
            for (const auto &edge1 : weighted_graph.edges[i]) {
                for (const auto &edge2 : weighted_graph.edges[i]) {
                    if (edge1.to == edge2.to) break;
                    auto other_edge = std::minmax(edge1.to, edge2.to);
                    if (S[edge1.to] && S[edge2.to] && i < edge1.to && i < edge2.to && edge2weight.count(other_edge) > 0) {
                        sum_weight_of_S += edge1.weight + edge2.weight + edge2weight[other_edge];
                    }
                }
            }
        }
    }

    return std::make_pair(maximum_weight_triangle_density, sum_weight_of_S / (double) n_vertices_of_S);
}


