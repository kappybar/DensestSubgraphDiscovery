
#include <vector>
#include <queue>
#include <limits>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cassert>
#include "exact_maxflow.hpp"

double EPS = 1e-4;


template<typename T> 
struct MaxFlowGraph {
private:
    struct Edge {
        int to, rev;
        T cap;
        Edge(int to, int rev, T cap):to(to), rev(rev), cap(cap) {}
    };
    int n;
    std::vector<std::vector<Edge>> G;
    std::vector<int> level;
    std::vector<int> iter;

    void bfs(int s) {
        fill(level.begin(), level.end(), -1);
        std::queue<int> q;
        level[s] = 0;
        q.push(s);
        while(!q.empty()) {
            int v = q.front();
            q.pop();
            for (const Edge &e : G[v]) {
                if (e.cap <= EPS || level[e.to] >= 0) continue;
                level[e.to] = level[v] + 1;
                q.push(e.to);
            }
        }
        return;
    }

    T dfs(int v, int t, T f) {
        if (v == t) return f;

        for(int &i = iter[v];i < (int)G[v].size(); i++) {
            Edge &e = G[v][i];
            if (e.cap <= EPS || level[v] >= level[e.to]) continue;
            T d = dfs(e.to, t, std::min(f, e.cap));
            if (d > EPS) {
                e.cap -= d;
                G[e.to][e.rev].cap += d;
                return d;
            }
        }

        return 0;
    }

public:
    MaxFlowGraph(int n):n(n), G(n), level(n), iter(n) {}

    void add_edge(int from, int to, T cap) {
        int rev_index_to = (int) G[to].size();
        int rev_index_from = (int) G[from].size();
        G[from].emplace_back(to, rev_index_to, cap);
        G[to].emplace_back(from, rev_index_from, 0);
        return;
    }

    bool reachable_from_s(int s) {
        for (const Edge &edge : G[s]) {
            if (edge.cap > EPS) return true;
        }
        return false;
    }

    T max_flow (int s, int t) {
        T flow = 0;
        T inf = std::numeric_limits<T>::max();
        while (1) {
            bfs(s);
            if(level[t] < 0) return flow;
            fill(iter.begin(), iter.end(), 0);
            T f;
            while((f = dfs(s, t, inf)) > EPS) {
                flow += f;
            }
        }
    }
};

double exact_dsp_using_flow(const Graph &graph) {
    int n_vertices = graph.n_vertices;
    int n_edges = 0;

    for (int i = 0;i < n_vertices; i++) {
        n_edges += (int) graph.edges[i].size();
    }
    n_edges /= 2;

    std::map<std::pair<int, int>, int> edgeid;
    int id = n_vertices;
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &edge : graph.edges[i]) {
            if (i < edge.to) {
                edgeid[std::make_pair(i, edge.to)] = id++;
            }
        }
    }

    double l = 0, u = (n_vertices - 1) / 2.0;
    while ((u - l) * n_vertices * (n_vertices - 1) >= 1) {
        double alpha = (u + l) / 2;
        MaxFlowGraph<double> max_flow_graph(n_vertices + n_edges + 2);

        int s = n_vertices + n_edges;
        int t = n_vertices + n_edges + 1;
        for (int i = 0;i < n_vertices; i++) {
            max_flow_graph.add_edge(s, i, (int)graph.edges[i].size());
            max_flow_graph.add_edge(i, t, 2 * alpha);
            for (const auto &edge : graph.edges[i]) {
                int eid = edgeid[std::minmax(i, edge.to)];
                max_flow_graph.add_edge(i, eid, 1);
                max_flow_graph.add_edge(eid, i, 1);
            }
        }

        max_flow_graph.max_flow(s, t);
        if (max_flow_graph.reachable_from_s(s)) {
            l = alpha;
        } else {
            u = alpha;
        }
        std::cout << l << " < opt < " << u << std::endl;

    }

    return l;
}

double exact_weighted_dsp_using_flow(const WeightedGraph &graph) {
    int n_vertices = graph.n_vertices;
    int n_edges = 0;

    for (int i = 0;i < n_vertices; i++) {
        n_edges += (int) graph.edges[i].size();
    }
    n_edges /= 2;

    std::map<std::pair<int, int>, int> edgeid;
    int id = n_vertices;
    long long weight_sum = 0;
    std::vector<long long> weighted_degree(n_vertices, 0);
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &edge : graph.edges[i]) {
            if (i < edge.to) {
                edgeid[std::make_pair(i, edge.to)] = id++;
            }
            weighted_degree[i] += edge.weight;
        }
        weight_sum += weighted_degree[i];
    }
    weight_sum /= 2;


    double l = 0, u = weight_sum;
    while ((u - l) * n_vertices * (n_vertices - 1) >= 1) {
        double alpha = (u + l) / 2;
        MaxFlowGraph<double> max_flow_graph(n_vertices + n_edges + 2);

        int s = n_vertices + n_edges;
        int t = n_vertices + n_edges + 1;
        for (int i = 0;i < n_vertices; i++) {
            max_flow_graph.add_edge(s, i, weighted_degree[i]);
            max_flow_graph.add_edge(i, t, 2 * alpha);
            for (const auto &edge : graph.edges[i]) {
                int eid = edgeid[std::minmax(i, edge.to)];
                max_flow_graph.add_edge(i, eid, edge.weight);
                max_flow_graph.add_edge(eid, i, edge.weight);
            }
        }

        max_flow_graph.max_flow(s, t);
        if (max_flow_graph.reachable_from_s(s)) {
            l = alpha;
        } else {
            u = alpha;
        }
        std::cout << l << " < opt < " << u << std::endl;

    }

    return l;
}

double exact_tdsp_using_flow(const Graph &graph) {
    int n_vertices = graph.n_vertices;
    int n_triangles = 0;

    for (int i = 0;i < n_vertices; i++) {
        n_triangles += (int) graph.triangles[i].size();
    }
    n_triangles /= 3;

    std::map<std::tuple<int, int, int>, int> triangleid;
    int id = n_vertices;
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &triangle : graph.triangles[i]) {
            if (i < triangle.vertex1 && i < triangle.vertex2) {
                triangleid[std::make_tuple(i, triangle.vertex1, triangle.vertex2)] = id++;
            }
        }
    }

    std::vector<std::vector<int>> connected_triangle(n_triangles);
    for (const auto &[triangle, id] : triangleid) {
        auto [a, b, c] = triangle;
        connected_triangle[a].push_back(id);
        connected_triangle[b].push_back(id);
        connected_triangle[c].push_back(id);
    }

    double l = 0, u = (n_vertices - 1) * (n_vertices - 2) / 6.0;
    while ((u - l) * n_vertices * (n_vertices - 1) >= 1) {
        double alpha = (u + l) / 2;
        MaxFlowGraph<double> max_flow_graph(n_vertices + n_triangles + 2);

        int s = n_vertices + n_triangles;
        int t = n_vertices + n_triangles + 1;
        for (int i = 0;i < n_vertices; i++) {
            max_flow_graph.add_edge(s, i, (int)graph.triangles[i].size());
            max_flow_graph.add_edge(i, t, 3 * alpha);
            for (int tid : connected_triangle[i]) {
                max_flow_graph.add_edge(i, tid, 1);
                max_flow_graph.add_edge(tid, i, 2);
            }
        }

        max_flow_graph.max_flow(s, t);
        if (max_flow_graph.reachable_from_s(s)) {
            l = alpha;
        } else {
            u = alpha;
        }
        std::cout << l << " < opt < " << u << std::endl;
    }

    return l;
}


double exact_weighted_tdsp_using_flow(const WeightedGraph &graph) {
    int n_vertices = graph.n_vertices;
    int n_triangles = 0;

    for (int i = 0;i < n_vertices; i++) {
        n_triangles += (int) graph.triangles[i].size();
    }
    n_triangles /= 3;

    std::map<std::tuple<int, int, int>, int> triangleid;
    int id = n_vertices;
    long long weight_sum = 0;
    std::vector<long long> weighted_triangle(n_vertices, 0);
    for (int i = 0;i < n_vertices; i++) {
        for (const auto &triangle : graph.triangles[i]) {
            if (i < triangle.vertex1 && i < triangle.vertex2) {
                if (triangle.vertex1 < triangle.vertex2) {
                    triangleid[std::make_tuple(i, triangle.vertex1, triangle.vertex2)] = id++;
                } else {
                    triangleid[std::make_tuple(i, triangle.vertex2, triangle.vertex1)] = id++;
                }
            }
            weighted_triangle[i] += triangle.weight;
        }
        weight_sum += weighted_triangle[i];
    }
    weight_sum /= 3;

    double l = 0, u = weight_sum;
    while ((u - l) * n_vertices * (n_vertices - 1) >= 1) {
        double alpha = (u + l) / 2;
        MaxFlowGraph<double> max_flow_graph(n_vertices + n_triangles + 2);

        int s = n_vertices + n_triangles;
        int t = n_vertices + n_triangles + 1;
        for (int i = 0;i < n_vertices; i++) {
            max_flow_graph.add_edge(s, i, weighted_triangle[i]);
            max_flow_graph.add_edge(i, t, 3 * alpha);
            for (const auto& triangle : graph.triangles[i]) {
                int vs[] = {i, triangle.vertex1, triangle.vertex2};
                std::sort(vs, vs + 3);
                int tid = triangleid[std::make_tuple(vs[0], vs[1], vs[2])];
                max_flow_graph.add_edge(i, tid, triangle.weight);
                max_flow_graph.add_edge(tid, i, 2 * triangle.weight);
            }
        }

        max_flow_graph.max_flow(s, t);
        if (max_flow_graph.reachable_from_s(s)) {
            l = alpha;
        } else {
            u = alpha;
        }
        std::cout << l << " < opt < " << u << std::endl;
    }

    return l;
}


