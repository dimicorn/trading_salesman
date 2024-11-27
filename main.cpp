#include <format>
#include <chrono>
#include "graph.hpp"


std::vector<std::pair<double, double> > parse(const std::string &input_file) {
    double x, y;
    std::vector<std::pair<double, double> > dots;
    std::ifstream input(input_file);
    while (input >> x >> y)
        dots.push_back(std::make_pair(x, y));
    input.close();
    return dots;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "No input file!";
        return 1;
    }
    std::string input_file = argv[1];
    std::vector<std::pair<double, double> > dots = parse(input_file);
    int vertices = dots.size();
    std::cout << "Number of vertices: " << vertices << '\n';
    std::vector<Edge> edges;
    for (int i = 0; i < vertices; ++i) {
        for (int j = i+1; j < vertices; ++j) {
            double ax = dots[i].first, ay = dots[i].second;
            double bx = dots[j].first, by = dots[j].second;
            double w = std::sqrt(pow(ax-bx, 2) + pow(ay-by, 2));
            Edge e = Edge(i, j, w);
            edges.push_back(e);
        }
    }

    Graph g = Graph(edges, vertices);
    const auto start1 = std::chrono::high_resolution_clock::now();
    auto [length, cycle] = g.ant_colony();
    const auto end1 = std::chrono::high_resolution_clock::now();
    std::cout << "Ant Colony Optimization Results:\n";
    std::cout << std::format("Optimal cycle length: {:.2f}\n", length);
    std::cout << "Optimal cycle path: ";
    for (const auto &it: cycle)
        std::cout << it << ' ';
    std::cout << '\n';
    const std::chrono::duration<double> diff1 = end1 - start1;
    std::cout << std::format("Time used: {}\n", diff1);

    if (vertices <= 25) {
        const auto start2 = std::chrono::high_resolution_clock::now();
        std::tie(length, cycle) = g.held_karp();
        const auto end2 = std::chrono::high_resolution_clock::now();
        std::cout << "\nHeld-Karp Algorithm Results:\n";
        std::cout << std::format("Optimal cycle length: {:.2f}\n", length);
        std::cout << "Optimal cycle path: ";
        for (const auto &it: cycle)
            std::cout << it << ' ';
        std::cout << '\n';
        const std::chrono::duration<double> diff2 = end2 - start2;
        std::cout << std::format("Time used: {}\n", diff2);
    }

    return 0;
}
