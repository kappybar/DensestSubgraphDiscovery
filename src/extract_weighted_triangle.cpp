
#include "extract_weighted_triangle.hpp"
#include <limits>
#include <set>
#include <algorithm>
#include <cmath>
#include <map>
#include <iostream>
#include <queue>


struct WeightedTriangle {
    int a, b, c;
    long long weight;
    WeightedTriangle(int a, int b, int c, long long weight):a(a), b(b), c(c), weight(weight) {}
    inline bool operator<(const WeightedTriangle& t) const {
        return weight > t.weight;
    }
};

struct WeightedEdge {
    int a, b;
    long long weight;
    WeightedEdge(int a, int b, long long weight):a(a), b(b), weight(weight) {}
};

// aproximately extract top-k weighted triangle
// assume graph.triangle is not caluculated when this function is called
// add these triangles to 1st argument "graph"
void dynamic_heavy_light(WeightedGraph &graph, int k, double alpha) {

    // sort edge in decreasing order of weight
    std::vector<WeightedEdge> edges;
    edges.reserve(10 * graph.n_vertices);
    for (int i = 0;i < graph.n_vertices; i++) {
        for (const auto &edge : graph.edges[i]) {
            if (i < edge.to) {
                edges.emplace_back(i, edge.to, edge.weight);
            }
        }
    }
    std::sort(edges.begin(), edges.end(), [](const auto &left, const auto &right) {
        return left.weight > right.weight;
    });

    std::vector<std::map<int, long long>> exists(graph.n_vertices);
    std::vector<long long> vert_to_wt(graph.n_vertices);
    std::vector<bool> computed(graph.n_vertices);
    std::vector<std::set<int>> deleted(graph.n_vertices); // edge deleted from light

    std::set<WeightedTriangle> counter, topk;
    long long num_tris = 0;

    std::vector<std::vector<WeightedEdge>> Gh; // heavy Graph
    int hi = 0, hj = 0;
    long long threshold = std::numeric_limits<long long>::max();
    counter.insert(WeightedTriangle(-1, -1, -1, threshold));
    auto curr = counter.begin();

    auto compute_exists_per_node = [&](int u) {
        if (!computed[u]) {
            computed[u] = true;
            for (const auto& edge : graph.edges[u]) {
                if (!deleted[u].count(edge.to)) {
                    exists[u][edge.to] = edge.weight;
                    exists[edge.to][u] = edge.weight;
                }
            }
        } 
    };

    auto search = [&](int u, int v) {
        if (graph.edges[u].size() > graph.edges[v].size()) std::swap(u, v);
        if (deleted[u].count(v)) return 0LL;
        if (exists[u].count(v)) return exists[u][v];
        for (const auto& edge : graph.edges[u]) {
            if (edge.to == v) return (long long) edge.weight;
        }
        return 0LL;
    };

    auto insert_triangle = [&](long long weight, long long threshold, WeightedTriangle t) {
        if (weight >= threshold) {
            topk.insert(t);
        } else {
            auto it = counter.insert(t).first;
            if (*it < *curr) {
                curr = it;
            }
        }
    };

    while ((int) topk.size() < k+1 && hj < (int) edges.size() ) {
        auto ei = edges[hi];
        auto ej = edges[hj];
        Gh.resize(std::max(ej.b + 1, (int) Gh.size()));
        threshold = 2 * ej.weight + ei.weight;

        bool advance_j = std::pow(ej.weight, alpha) >= ei.weight && hj < (int) edges.size() - 1;

        // always hi <= hj
        if (hi == hj || advance_j) {

            // from ej.a get (e : heavy edge)
            for (const auto& e : Gh[ej.a]) {
                long long weight = search(e.b, ej.b);
                if (weight > 0) {
                    WeightedTriangle t(e.b, ej.b, ej.a, e.weight + ej.weight + weight);
                    insert_triangle(t.weight, threshold, t);
                    num_tris ++;
                }
            }

            // from ej.b get (e : heavy edge)
            for (const auto& e : Gh[ej.b]) {
                long long weight = search(e.b, ej.a);
                if (weight > 0) {
                    WeightedTriangle t(e.b, ej.b, ej.a, e.weight + ej.weight + weight);
                    insert_triangle(t.weight, threshold, t);
                    num_tris ++;
                }
            }

            // from ej.a get heavy edge &&
            // from ej.b get heavy edge
            for (const auto& e : Gh[ej.a]) {
                vert_to_wt[e.b] = e.weight;
            }
            for (const auto& e : Gh[ej.b]) {
                if (vert_to_wt[e.b] > 0) {
                    WeightedTriangle t(e.b, ej.a, ej.b, ej.weight + e.weight + vert_to_wt[e.b]);
                    insert_triangle(t.weight, threshold, t);
                    num_tris++;
                }
            }
            for (const auto& e : Gh[ej.a]) {
                vert_to_wt[e.b] = 0;
            }

            hj++;

            deleted[ej.a].insert(ej.b);
            deleted[ej.b].insert(ej.a);
            Gh[ej.a].emplace_back(ej.a, ej.b, ej.weight);
            Gh[ej.b].emplace_back(ej.b, ej.a, ej.weight);
        } else {

            // 1 heavy edge
            compute_exists_per_node(ei.a);
            compute_exists_per_node(ei.b);
            for (const auto& [v, w] : exists[ei.a]) {
                vert_to_wt[v] = w;
            }
            for (const auto& [v, w] : exists[ei.b]) {
                if (vert_to_wt[v]) {
                    long long weight = ei.weight + w + vert_to_wt[v];
                    WeightedTriangle t(v, ei.a, ei.b, weight);
                    insert_triangle(weight, threshold, t);
                    num_tris ++;
                }
            }
            for (const auto& [v, w] : exists[ei.a]) {
                vert_to_wt[v] = 0;
            }
            hi++;

        }

        while (curr != counter.end() && curr->weight >= threshold) {
            topk.insert(*curr);
            auto prev = curr;
            curr++;
            if (prev != counter.begin()) {
                counter.erase(prev);
            }
        }
        if (curr == counter.end()) {
            curr--;
        }

    }

    counter.erase(counter.begin());
    topk.erase(topk.begin());

    std::cout << "Found" << num_tris << " triangles" << std::endl;
    std::cout << "Out of these, the top" << (int) topk.size() << " are found for sure." << std::endl;
    if (topk.size()) std::cout << "maximum weight : " << topk.begin()->weight << std::endl;

    for (const auto &triangle : topk) {
        auto [a, b, c, w] = triangle;
        graph.triangles[a].emplace_back(b, c, w);
        graph.triangles[b].emplace_back(a, c, w);
        graph.triangles[c].emplace_back(a, b, w);
    }

    return;
}


// weighted_graph have only k-top wieghted triangle
// let S := subgraph that maximize k-top-triangle-density
// let a := k-top-triangle-densty(S)
// let b := triangle-density(S)
// return (a, b)
std::vector<std::pair<double, double>> estimate_ktop_triangle_greedy_plusplus(WeightedGraph &weighted_graph, long long interval_k, long long last_k, int T) {

    // collect all triangle first, but don't add triangle to graph
    std::vector<WeightedTriangle> all_triangle;
    std::map<std::pair<int, int>, long long> all_edges;
    for (int i = 0;i < weighted_graph.n_vertices; i++) {
        for (const auto &edge : weighted_graph.edges[i]) {
            if (i < edge.to) {
                all_edges[std::make_pair(i, edge.to)] = edge.weight;
            }
        }
    }
    for (int i = 0;i < weighted_graph.n_vertices; i++) {
        for (const auto &edge1 : weighted_graph.edges[i]) {
            for (const auto &edge2 : weighted_graph.edges[i]) {
                if (edge1.to == edge2.to) break;
                auto other_edge = std::minmax(edge1.to, edge2.to);
                if (i < edge1.to && i < edge2.to && all_edges.count(other_edge) > 0) {
                    all_triangle.emplace_back(i, other_edge.first, other_edge.second, edge1.weight + edge2.weight + all_edges[other_edge]);
                }
            }
        }
    }
    std::sort(all_triangle.begin(), all_triangle.end());

    auto compare = [](const std::pair<long long, int> &a, const std::pair<long long, int> &b) {
        return a.first > b.first;
    };

    std::vector<std::pair<double, double>> results;

    int index = 0;
    for (int k = interval_k;k <= last_k;k += interval_k) {
        int ik = index;
        for (;ik < k && ik < (int) all_triangle.size(); ik++) {
            auto [a, b, c, w] = all_triangle[ik];
            weighted_graph.add_triangle(a, b, c, w);
        }
        index = ik;

        std::priority_queue<std::pair<long long, int>,
                            std::vector<std::pair<long long, int>>,
                            decltype(compare)> pq {compare};

        std::vector<long long> load(weighted_graph.n_vertices, 0);
        double maximum_weight_triangle_density = 0;
        std::vector<bool> S(weighted_graph.n_vertices, 0);

        for (int t = 0;t < T; t++) {
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

                    load[i] += triangle_edge.weight;
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

        results.emplace_back(index, sum_weight_of_S / (double) n_vertices_of_S);
    }

    return results;
}