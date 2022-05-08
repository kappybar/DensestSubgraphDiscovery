#include "graph_reader.hpp"
#include <fstream>
#include <set>
#include <utility>
#include <algorithm>

Graph read_unweighted_graph(const std::string &file_name) {
    std::ifstream input_file(file_name);

    int n_vertices, n_edges;
    input_file >> n_vertices >> n_edges;

    Graph unweighted_graph(n_vertices);

    std::set<std::pair<int, int>> edges;
    for (int i = 0;i < n_edges; i++) {
        int u, v;
        input_file >> u >> v;

        // delete self-loop, multiple edge
        if (u == v) continue; 
        if (edges.count(std::minmax(u, v)) > 0) continue;

        unweighted_graph.add_edge(u, v);
        edges.emplace(std::minmax(u, v));
    }

    input_file.close();

    return unweighted_graph;
}

WeightedGraph read_weighted_graph(const std::string &file_name) {
    std::ifstream input_file(file_name);

    int n_vertices, n_edges;
    input_file >> n_vertices >> n_edges;

    WeightedGraph weighted_graph(n_vertices);

    std::set<std::pair<int, int>> edges;
    for (int i = 0;i < n_edges; i++) {
        int u, v;
        double weight;
        input_file >> u >> v >> weight;

        // delete self-loop, multiple edge
        if (u == v) continue; 
        if (edges.count(std::minmax(u, v)) > 0) continue;

        weighted_graph.add_edge(u, v, weight);
        edges.emplace(std::minmax(u, v));
    }

    input_file.close();

    return weighted_graph;
}