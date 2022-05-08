#pragma once
#include "graph.hpp"

double exact_dsp_using_flow(const Graph &graph);
double exact_weighted_dsp_using_flow(const WeightedGraph &graph);
double exact_tdsp_using_flow(const Graph &graph);
double exact_weighted_tdsp_using_flow(const WeightedGraph &graph);

