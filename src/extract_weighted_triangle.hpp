#pragma once
#include "graph.hpp"
#include <utility>
#include <vector>

// it is not used in experiment
void dynamic_heavy_light(WeightedGraph &graph, int k, double alpha = 1.25);

std::vector<std::pair<double, double>> estimate_ktop_triangle_greedy_plusplus(WeightedGraph &weighted_graph, long long interval_k, long long last_k, int T);