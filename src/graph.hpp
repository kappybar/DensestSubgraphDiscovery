#pragma once
#include <vector>


struct WeightedGraph {
    struct WeightedEdge {
        int to;
        long long weight;
        WeightedEdge(int to, long long weight);
    };
    struct WeightedTriangleEdge {
        int vertex1;
        int vertex2;
        long long weight;
        WeightedTriangleEdge(int vertex1, int vertex2, long long weight);
    };
    int n_vertices;
    std::vector<std::vector<WeightedEdge>> edges;
    std::vector<std::vector<WeightedTriangleEdge>> triangles;
    std::vector<std::vector<int>> vertices_connected_by_triangle;

    WeightedGraph(int n_vertices);
    void add_edge(int u, int v, long long weight);
    void extract_triangle(void);
    double greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps);
    double triangle_greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps);
    void add_triangle(int a, int b, int c, long long weight);   
    
};




struct Graph {
    struct Edge {
        int to;
        Edge(int to);
    };
    struct TriangleEdge {
        int vertex1;
        int vertex2;
        TriangleEdge(int vertex1, int vertex2);
    };
    int n_vertices;
    std::vector<std::vector<Edge>> edges;
    std::vector<std::vector<TriangleEdge>> triangles;
    std::vector<std::vector<int>> vertices_connected_by_triangle;

    Graph(int n_vertices);
    void add_edge(int u, int v);
    WeightedGraph put_uniform_weight(void); 
    void extract_triangle(void);
    double greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps);
    double triangle_greedy_plusplus_convergence_order(double lamda_star_lower_bound, double eps);
    void add_triangle(int a, int b, int c);

};

struct Data {
    double density;         // e(S) / |S|
    double triangle_density; // t(S) / |S|
    double edge_density;    // e(S) / (|S|, 2)
    double vertex_ratio;    // |S| / |V|

    Data(double density, double triangle_density, double edge_density, double vertex_ratio);
    void display(void);
};