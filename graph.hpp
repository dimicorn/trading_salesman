#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <string>
#include <cmath>
#include <random>
#include <utility>

struct Edge {
    int a, b;
    double w;
    Edge(const int &a, const int &b, const double &w):
        a(a), b(b), w(w) {}
};

class Graph {
public:
    Graph(const std::vector<Edge> &edges, const int &vertices) {
        std::vector<std::vector<double> > matrix(vertices, std::vector<double> (vertices, 0));
        double mean = 0;
        for (const auto &e: edges) {
            matrix[e.a][e.b] = e.w;
            matrix[e.b][e.a] = e.w;
            mean += 2 * e.w;
        }
        mean = mean / (vertices * vertices) * 10;
        std::vector<std::vector<double> > mean_matrix(vertices, std::vector<double> (vertices, mean));
        intensity = std::move(mean_matrix);
        adj_matrix = std::move(matrix);
        num_vertices = vertices;
    }
    std::pair<double, std::vector<int> > held_karp() {
        min_path.resize(1 << num_vertices, std::vector<double>(num_vertices, MAX));
        min_path_last_vertex.resize(1 << num_vertices, std::vector<int>(num_vertices, -1));
        min_path[0][0] = 0;  // Our starting point
        
        // Iterating over sets this way guarantees we explore smaller sets first
        for (int explored_set = 0; explored_set < (1U << num_vertices); ++explored_set) {
            for (int last_vertex = 0; last_vertex < num_vertices; ++last_vertex) {
                if (min_path[explored_set][last_vertex] != MAX) {  // Valid state?
                    for (int next_vertex = 0; next_vertex < num_vertices; ++next_vertex) {
                        if (adj_matrix[last_vertex][next_vertex] != MAX // Edge exists
                            && (explored_set & (1U << next_vertex)) == 0) {  // We don't go back to already explored vertex
                            // Perform the traversal and see if we find a new shortest path/cycle to the destination state
                            const int next_set = explored_set + (1U << next_vertex);
                            const double next_cost = min_path[explored_set][last_vertex] + adj_matrix[last_vertex][next_vertex];
                            if (next_cost < min_path[next_set][next_vertex]) {
                                min_path[next_set][next_vertex] = next_cost;
                                min_path_last_vertex[next_set][next_vertex] = last_vertex;
                            }
                        }
                    }
                }
            }
        }
        double result_cost = min_path[(1 << num_vertices) - 1][0];
        std::vector<int> cycle = {0};
        if (result_cost != MAX) {
            int cur_state = (1 << num_vertices) - 1;
            int cur_vertex = 0;
            
            while (cur_state > 0) {  // Traverse back the minimum cycle
                int prev_vertex = min_path_last_vertex[cur_state][cur_vertex];
                cur_state -= (1 << cur_vertex);
                cur_vertex = prev_vertex;
                cycle.push_back(cur_vertex);
            }
            std::reverse(cycle.begin(), cycle.end());
        }
        return std::make_pair(result_cost, cycle);
    }
    std::pair<double, std::vector<int> > ant_colony(
        const int &iterations = 100, const int &ants_per_iter = 50,
        const int &q = 10, const double &degradation_factor = 0.9
    ) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_int_distribution<int> uni(0, num_vertices-1);
        std::vector<int> best_cycle(num_vertices, -1);
        auto best_length = MAX;

        for (int i = 0; i < iterations; ++i) {
            std::vector<std::pair<double, std::vector<int> > > cycles;
            for (int ant = 0; ant < ants_per_iter; ++ant) {
                auto cycle = traverse(uni(rng));
                cycles.push_back(cycle);
            }
            std::sort(cycles.begin(), cycles.end(), [](const auto &a, const auto &b) {
                    return a.first < b.first;
            });
            cycles = std::move(std::vector<std::pair<double, std::vector<int> > >(cycles.begin(), cycles.begin() + ants_per_iter/2));

            for (const auto &it: cycles) {
                double total_length = it.first;
                std::vector<int> cycle = it.second;
                if (total_length < best_length) {
                    best_length = total_length;
                    best_cycle = cycle;
                }
                double delta = q / total_length;
                int n = cycle.size();
                for (int i = 0; i < n-1; ++i)
                    intensity[cycle[i]][cycle[i+1]] += delta;
                intensity[cycle[n-1]][cycle[0]] += delta;
                std::for_each(intensity.begin(), intensity.end(), [&](auto &v) {
                    std::for_each(v.begin(), v.end(), [&](auto &el) {
                        el *= degradation_factor;
                    });
                });
            }
        }
        return std::make_pair(best_length, best_cycle);
    }
    std::pair<double, std::vector<int> > traverse(const int &source_node = 0) {
        std::vector<int> visited(num_vertices, 1);
        visited[source_node] = 0;
        std::vector<int> cycle = {source_node};
        int steps = 0, current = source_node;
        double total_length = 0;
        std::default_random_engine generator;
        while (steps < num_vertices - 1) {
            std::vector<double> jump_neighbors, jump_values;
            for (int node = 0; node < num_vertices; ++node) {
                if (visited[node] != 0) {
                    double pheromone_level = fmax(intensity[current][node], 1e-5);
                    double value = pow(pheromone_level, alpha) / pow(adj_matrix[current][node], beta);
                    jump_neighbors.push_back(node);
                    jump_values.push_back(value);
                }
            }
            std::discrete_distribution<int> distribution(jump_values.begin(), jump_values.end());
            int next_node = distribution(generator);
            visited[jump_neighbors[next_node]] = 0;
            current = jump_neighbors[next_node];
            cycle.push_back(current);
            ++steps;
        }
        total_length = cycle_length(cycle);
        return std::make_pair(total_length, cycle);
    }
    double cycle_length(const std::vector<int> &cycle) {
        double length = 0;
        int n = cycle.size();
        for (int i = 0; i < n - 1; ++i)
            length += adj_matrix[cycle[i]][cycle[i+1]];
        length += adj_matrix[cycle[n-1]][cycle[0]];
        return length;
    }

private:
    std::vector<std::vector<double> > adj_matrix, min_path, intensity;
    std::vector<std::vector<int> > min_path_last_vertex;
    int num_vertices;
    double alpha = 0.9, beta = 1.5;
    const double MAX = std::numeric_limits<double>::max();
};