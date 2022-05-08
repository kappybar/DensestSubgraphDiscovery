#include "graph.hpp"
#include "graph_reader.hpp"
#include "greedy.hpp"
#include "greedy_plusplus.hpp"
#include "exact_maxflow.hpp"
#include "extract_weighted_triangle.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <fstream>

// ./main <input-file> <mode>
int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "specify file or mode!" << std::endl;
        return 0;
    }
    std::string file_name(argv[1]);
    std::string mode(argv[2]);

    std::cout << std::setprecision(10);
    std::cout << "file-name : " << file_name << std::endl;

    // decide weighted or unweighted
    bool weighted = file_name.find("weighted") != std::string::npos;

    // speed up
    int alpha = 10;

    Graph g(0); 
    WeightedGraph wg(0);
    if (weighted) {
        wg = read_weighted_graph(file_name);
    } else {
        g = read_unweighted_graph(file_name);
        assert(mode.find("weighted") == std::string::npos);
    }

    double eps = 0.01;
    int iter_count = 100;

    // Triangle density 
    if (mode == "tdsp-greedy") {
        g.extract_triangle();
        double t = pq_triangle_greedy(g);
        std::cout << "T density :" << 3*t << std::endl; 
        double iter = g.triangle_greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to need 0.01 approximation :" << iter << std::endl << std::endl;
    } 
    else if (mode == "tdsp-greedy++") {
        g.extract_triangle();
        double t = pq_triangle_greedy_plusplus(g, iter_count);
        std::cout << "T density : " << 3*t << std::endl;
        double iter = g.triangle_greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to need 0.01 approximation :" << iter << std::endl << std::endl;
    }
    else if (mode == "tdsp-exact") {
        g.extract_triangle();
        double t = exact_tdsp_using_flow(g);
        std::cout << "T density : " << 3*t << std::endl;
    }
    // Weighted Triangle density
    else if (mode == "weighted-tdsp-greedy") {
        wg.extract_triangle();
        double t = pq_weighted_triangle_greedy(wg);
        std::cout << "T density : " << t << std::endl;
        double iter = wg.triangle_greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to nned 0.01 approximation : " << iter << std::endl << std::endl;
    }
    else if (mode == "weighted-tdsp-greedy++") {
        wg.extract_triangle();
        double t = pq_weighted_triangle_greedy_plusplus(wg, iter_count);
        std::cout << "T density : " << t << std::endl;
        double iter = wg.triangle_greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to nned 0.01 approximation : " << iter << std::endl << std::endl;
    }
    else if (mode == "weighted-tdsp-exact") {
        wg.extract_triangle();
        double t = exact_weighted_tdsp_using_flow(wg);
        std::cout << "T density : " << t << std::endl;
        double iter = wg.triangle_greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to nned 0.01 approximation : " << iter << std::endl << std::endl;
    }
    // Only top-k weighted triangle density
    else if (mode == "tdsp-topk-weighted") {
        std::ofstream ofs("output.txt");
        auto results = estimate_ktop_triangle_greedy_plusplus(wg, 50000, 1700000, 10);
        for (auto [k, t_density] : results) {
            ofs << k << " " << t_density << std::endl;
        }
    }
    // Density
    else if (mode == "dsp-greedy") {
        Data result = pq_greedy(g);
        result.display();
        double iter = g.greedy_plusplus_convergence_order(result.density, eps);
        std::cout << "iter to need 0.01 approximation :" << iter << std::endl << std::endl;
    } 
    else if (mode == "dsp-greedy++") {
        Data result = pq_greedy_plusplus(g, iter_count);
        result.display();
        double iter = g.greedy_plusplus_convergence_order(result.density, eps);
        std::cout << "iter to need 0.01 approximation :" << iter << std::endl << std::endl;
    }
    else if (mode == "dsp-exact") {
        double t = exact_dsp_using_flow(g);
        std::cout << "density : " << 2*t << std::endl;
        double iter = g.greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to need 0.01 approximation :" << iter << std::endl << std::endl;
    }
    // Weighted Density
    else if (mode == "weighted-dsp-greedy") {
        double t = pq_weighted_greedy(wg);
        std::cout << "density : " << t << std::endl;
        double iter = wg.greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to need 0.01 approximation : " << iter << std::endl << std::endl;
    }
    else if (mode == "weighted-dsp-greedy++") {
        double t = pq_weighted_greedy_plusplus(wg, iter_count);
        std::cout << "density : " << t << std::endl;
        double iter = wg.greedy_plusplus_convergence_order(t, eps);
        std::cout << "iter to need 0.01 approximation : " << iter << std::endl << std::endl;
    }
    else if (mode == "weighted-dsp-exact") {
        double t = exact_weighted_dsp_using_flow(wg);
        std::cout << "density :" << t << std::endl;
    }
    // add alpha * deg
    else if (mode == "dsp-mygreedy") {
        Data result = my_greedy_plusplus(g, iter_count, alpha);
        result.display();
        double iter = g.greedy_plusplus_convergence_order(result.density, eps);
        std::cout << "iter to need 0.01 approximation :" << iter << std::endl << std::endl;
    }
    else if (mode == "tdsp-mygreedy") {
        g.extract_triangle();
        double t = my_triangle_greedy_plusplus(g, iter_count, alpha);
        std::cout << "T density : " << 3*t << std::endl;
    }
   
    else {
        std::cerr << "Not Allowed mode!" << std::endl;
    }


    return 0;
}