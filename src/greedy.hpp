#pragma once

#include "graph.hpp"
#include <utility>

Data pq_greedy(const Graph &graph);
double pq_weighted_greedy(const WeightedGraph &weighted_graph);
double pq_triangle_greedy(const Graph &graph);
double pq_weighted_triangle_greedy(const WeightedGraph &weighted_graph);