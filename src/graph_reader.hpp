#pragma once
#include "graph.hpp"
#include <string>

Graph read_unweighted_graph(const std::string &file_name);
WeightedGraph read_weighted_graph(const std::string &file_name);