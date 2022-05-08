#pragma once
#include "graph.hpp"

Data pq_greedy_plusplus(const Graph &graph, int T);
double pq_weighted_greedy_plusplus(const WeightedGraph &weighted_graph, int T);
double pq_triangle_greedy_plusplus(const Graph &graph, int T);
double pq_weighted_triangle_greedy_plusplus(const WeightedGraph &weighted_graph, int T);

Data my_greedy_plusplus(const Graph &graph, int T, int alpha=10);
double my_triangle_greedy_plusplus(const Graph &graph, int T, int alpha=10);

// it is not used in experiment
std::pair<double, double> ktop_triangle_greedy_plusplus(const WeightedGraph &weighted_graph, int T);
