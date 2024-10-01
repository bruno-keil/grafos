#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "Node.hpp"
#include "defines.hpp"

class Graph
{
public:
    /*Assinatura dos métodos básicos para o funcionamento da classe*/

    Graph(std::ifstream& instance, bool directed, bool weighted_edges, bool weighted_nodes);
    Graph(bool directed, bool weighted_edges, bool weighted_nodes);
    ~Graph();

    void remove_node(size_t node_id);
    void remove_edge(size_t node_id_1, size_t node_id_2);
    void add_node(size_t node_id, float weight = 0);
    void add_edge(size_t node_id_1, size_t node_id_2, float weight = 0);
    void print_dot(std::ofstream& output_file);
    void print_graph();
    void dfs(size_t start_id);
    void dfsHelper(size_t node_id, std::unordered_set<size_t>& visited);
    void indirect_transitive_closure_helper(size_t node_id, std::unordered_set<size_t>& visited);
    std::unordered_set<size_t> indirect_transitive_closure(size_t start_id);
    bool isDirected();

    void dijkstra_shortest_path(size_t start_id, size_t end_id);
    void floyd_warshall_shortest_path(size_t start_id, size_t end_id);

    std::unordered_set<size_t> direct_transitive_closure(size_t start_id);
    void direct_transitive_closure_helper(size_t node_id, std::unordered_set<size_t>& visited);

    int connected(size_t node_id_1, size_t node_id_2);

    Graph prim_minimum_spanning_tree(const std::set<size_t>& subset_of_vertices);
    Graph kruskalMST(const std::set<size_t>& subset_of_vertices);

    void calculateGraphMetrics(Graph& graph);

    std::vector<std::set<size_t>> partition_graph(size_t p);
    void greedy_partition(size_t p);
    void grasp_partition(size_t p, size_t iterations);
    void reactive_grasp_partition(size_t p, size_t iterations);
    void print_partition(const std::vector<std::set<size_t>>& partition);
    size_t get_random_node(const std::unordered_set<size_t>& unvisited);
    float calculate_total_gap(const std::vector<std::set<size_t>>& partition);
    void update_performance(std::vector<float>& performance, float randomness, float gap);
    float adapt_randomness(const std::vector<float>& performance);
    std::vector<std::set<size_t>> local_search_refine(const std::vector<std::set<size_t>>& partition);
    std::vector<std::set<size_t>> randomized_greedy_partition(size_t p);
    std::vector<std::set<size_t>> randomized_grasp(float randomness, size_t p);

    void print_partition_weight(const std::vector<std::set<size_t>>& partition);
    size_t get_best_node(const std::set<size_t>& subgraph, const std::unordered_set<size_t>& unvisited);


private:
    size_t _number_of_nodes;
    size_t _number_of_edges;
    bool   _directed;
    bool   _weighted_edges;
    bool   _weighted_nodes;
    Node  *_first;
    Node  *_last;

    Node *findNode(size_t node_id);
};

#endif  //GRAPH_HPP
