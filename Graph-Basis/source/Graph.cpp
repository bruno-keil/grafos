#include "../include/Graph.hpp"

Graph::Graph(std::ifstream& instance, bool directed, bool weighted_edges, bool weighted_nodes)
{
    if (!instance)
    {
        throw std::runtime_error("Invalid file stream");
    }

    _directed       = directed;
    _weighted_edges = weighted_edges;
    _weighted_nodes = weighted_nodes;

    std::string firstLine;
    getline(instance, firstLine);
    std::istringstream iss(firstLine);

    size_t number_of_nodes;
    iss >> number_of_nodes;

    _number_of_nodes = number_of_nodes;
    _number_of_edges = 0;
    _first           = nullptr;
    _last            = nullptr;

    std::string line;
    while (getline(instance, line))
    {
        std::istringstream edgeStream(line);
        size_t             source, target;
        float              weight;
        edgeStream >> source >> target >> weight;
        add_node(source, 0);
        add_node(target, 0);
        add_edge(source, target, weight);
    }
}

Graph::Graph(bool directed, bool weighted_edges, bool weighted_nodes)
{
    _number_of_nodes = 0;
    _number_of_edges = 0;
    _directed        = directed;
    _weighted_edges  = weighted_edges;
    _weighted_nodes  = weighted_nodes;
    _first           = nullptr;
    _last            = nullptr;
}

Graph::~Graph()
{
    Node *current_node = _first;
    while (current_node != nullptr)
    {
        Node *next_node    = current_node->_next_node;
        Edge *current_edge = current_node->_first_edge;
        while (current_edge != nullptr)
        {
            Edge *next_edge = current_edge->_next_edge;
            delete current_edge;
            current_edge = next_edge;
        }
        delete current_node;
        current_node = next_node;
    }
}

void Graph::remove_node(size_t node_id)
{
    Node *current_node = _first;
    Node *prev_node    = nullptr;

    while (current_node != nullptr && current_node->_id != node_id)
    {
        prev_node    = current_node;
        current_node = current_node->_next_node;
    }

    if (current_node == nullptr)
        return;

    if (prev_node != nullptr)
        prev_node->_next_node = current_node->_next_node;
    else
        _first = current_node->_next_node;

    Node *iterating_node = _first;
    while (iterating_node != nullptr)
    {
        Edge *current_edge = iterating_node->_first_edge;
        Edge *prev_edge    = nullptr;

        while (current_edge != nullptr)
        {
            if (current_edge->_target_id == node_id)
            {
                if (prev_edge != nullptr)
                    prev_edge->_next_edge = current_edge->_next_edge;
                else
                    iterating_node->_first_edge = current_edge->_next_edge;

                Edge *edge_to_delete = current_edge;
                current_edge         = current_edge->_next_edge;
                delete edge_to_delete;
            }
            else
            {
                prev_edge    = current_edge;
                current_edge = current_edge->_next_edge;
            }
        }

        iterating_node = iterating_node->_next_node;
    }

    delete current_node;
    _number_of_nodes--;
}

void Graph::remove_edge(size_t node_id_1, size_t node_id_2)
{
    Node *node1 = findNode(node_id_1);

    if (node1 == nullptr)
        return;

    Edge *current_edge = node1->_first_edge;
    Edge *prev_edge    = nullptr;

    while (current_edge != nullptr && current_edge->_target_id != node_id_2)
    {
        prev_edge    = current_edge;
        current_edge = current_edge->_next_edge;
    }

    if (current_edge == nullptr)
        return;

    if (prev_edge != nullptr)
        prev_edge->_next_edge = current_edge->_next_edge;
    else
        node1->_first_edge = current_edge->_next_edge;

    delete current_edge;
    node1->_number_of_edges--;

    if (!_directed)
    {
        Node *node2 = findNode(node_id_2);

        if (node2 == nullptr)
            return;

        current_edge = node2->_first_edge;
        prev_edge    = nullptr;

        while (current_edge != nullptr && current_edge->_target_id != node_id_1)
        {
            prev_edge    = current_edge;
            current_edge = current_edge->_next_edge;
        }

        if (current_edge == nullptr)
            return;

        if (prev_edge != nullptr)
            prev_edge->_next_edge = current_edge->_next_edge;
        else
            node2->_first_edge = current_edge->_next_edge;

        delete current_edge;
        node2->_number_of_edges--;
    }

    _number_of_edges--;
}

void Graph::add_node(size_t node_id, float weight)
{
    Node *new_node             = new Node;
    new_node->_id              = node_id;
    new_node->_weight          = weight;
    new_node->_number_of_edges = 0;
    new_node->_first_edge      = nullptr;
    new_node->_next_node       = nullptr;
    new_node->_previous_node   = _last;

    if (_last != nullptr)
        _last->_next_node = new_node;

    _last = new_node;

    if (_first == nullptr)
        _first = new_node;

    _number_of_nodes++;
}

void Graph::add_edge(size_t node_id_1, size_t node_id_2, float weight) {
    Node *node1 = findNode(node_id_1);
    Node *node2 = findNode(node_id_2);

    // Ensure both nodes exist
    if (node1 == nullptr || node2 == nullptr) {
        std::cerr << "Error: One or both nodes do not exist.\n";
        return;
    }

    // Create the edge for node1 -> node2
    Edge *new_edge = new Edge;
    new_edge->_weight = weight;
    new_edge->_target_id = node_id_2;
    new_edge->_next_edge = node1->_first_edge;  // Link to the existing first edge
    node1->_first_edge = new_edge;  // Update the first edge to the new edge
    node1->_number_of_edges++;

    // If the graph is undirected, add a reverse edge node2 -> node1
    if (!_directed) {
        Edge *new_reverse_edge = new Edge;
        new_reverse_edge->_weight = weight;
        new_reverse_edge->_target_id = node_id_1;
        new_reverse_edge->_next_edge = node2->_first_edge;
        node2->_first_edge = new_reverse_edge;
        node2->_number_of_edges++;
    }

    _number_of_edges++;
}


void Graph::print_dot(std::ofstream& output_file)
{
    output_file << (_directed ? "digraph" : "graph") << " G {\n";
    
    Node *current_node = _first;
    
    while (current_node != nullptr)
    {
        output_file << "    " << current_node->_id;
        if (_weighted_nodes)
        {
            output_file << " [label=\"" << current_node->_id << " (" << current_node->_weight << ")\"]";
        }
        output_file << ";\n";

        Edge *current_edge = current_node->_first_edge;
        while (current_edge != nullptr)
        {
            output_file << "    " << current_node->_id 
                        << (_directed ? " -> " : " -- ") 
                        << current_edge->_target_id;

            if (_weighted_edges)
            {
                output_file << " [label=\"" << current_edge->_weight << "\"]";
            }
            output_file << ";\n";
            current_edge = current_edge->_next_edge;
        }
        current_node = current_node->_next_node;
    }
    
    output_file << "}\n";
}

void Graph::print_graph()
{
    std::cout << "Graph " << (_directed ? "Directed" : "Undirected") << ", " << (_weighted_edges ? "Weighted Edges" : "Unweighted Edges") << "\n";

    Node                      *current_node = _first;
    std::unordered_set<size_t> printedNodes;

    while (current_node != nullptr)
    {
        if (printedNodes.find(current_node->_id) == printedNodes.end())
        {
            printedNodes.insert(current_node->_id);

            std::cout << "Node " << current_node->_id;
            if (_weighted_nodes)
            {
                std::cout << " (Weight: " << current_node->_weight << ")";
            }
            std::cout << " -> ";

            Edge *current_edge = current_node->_first_edge;
            if (!current_edge)
            {
                std::cout << "No outgoing edges";
            }
            while (current_edge != nullptr)
            {
                std::cout << current_edge->_target_id;
                if (_weighted_edges)
                {
                    std::cout << " (Weight: " << current_edge->_weight << ")";
                }
                if (current_edge->_next_edge != nullptr)
                {
                    std::cout << ", ";
                }
                current_edge = current_edge->_next_edge;
            }
            std::cout << std::endl;
        }
        current_node = current_node->_next_node;
    }
}

int Graph::connected(size_t node_id_1, size_t node_id_2)
{
    std::queue<size_t>         q;
    std::unordered_set<size_t> visited;

    q.push(node_id_1);
    visited.insert(node_id_1);

    while (!q.empty())
    {
        size_t current = q.front();
        q.pop();

        if (current == node_id_2)
        {
            return 1;
        }

        Node *node = findNode(current);
        if (node)
        {
            Edge *edge = node->_first_edge;
            while (edge != nullptr)
            {
                if (visited.find(edge->_target_id) == visited.end())
                {
                    visited.insert(edge->_target_id);
                    q.push(edge->_target_id);
                }
                edge = edge->_next_edge;
            }
        }
    }
    return 0;
}

Node *Graph::findNode(size_t node_id)
{
    Node *current_node = _first;
    while (current_node != nullptr && current_node->_id != node_id)
    {
        current_node = current_node->_next_node;
    }
    return current_node;
}

bool Graph::isDirected()
{
    return _directed;
}

void Graph::dfs(size_t start_id)
{
    std::unordered_set<size_t> visited;
    dfsHelper(start_id, visited);
}

void Graph::dfsHelper(size_t node_id, std::unordered_set<size_t>& visited)
{
    visited.insert(node_id);
    Node *current_node = findNode(node_id);

    // Process edges from the current node
    for (Edge *edge = current_node->_first_edge; edge != nullptr; edge = edge->_next_edge)
    {
        size_t target_id = edge->_target_id;
        if (visited.find(target_id) == visited.end())
        {
            // This is a forward edge
            std::cout << "Edge: " << node_id << " -> " << target_id << std::endl;
            dfsHelper(target_id, visited);
        }
        else
        {
            // This is a back edge
            std::cout << "Back Edge: " << node_id << " -> " << target_id << std::endl;
        }
    }
}

std::unordered_set<size_t> Graph::direct_transitive_closure(size_t start_id)
{
    std::unordered_set<size_t> visited;
    direct_transitive_closure_helper(start_id, visited);
    return visited;
}

void Graph::direct_transitive_closure_helper(size_t node_id, std::unordered_set<size_t>& visited)
{
    // Mark the current node as visited
    if (visited.find(node_id) == visited.end())
    {
        visited.insert(node_id);

        // Recur for all vertices adjacent to this vertex
        Node *node = findNode(node_id);
        if (node)
        {
            Edge *edge = node->_first_edge;
            while (edge != nullptr)
            {
                if (visited.find(edge->_target_id) == visited.end())
                {
                    direct_transitive_closure_helper(edge->_target_id, visited);
                }
                edge = edge->_next_edge;
            }
        }
    }
}

std::unordered_set<size_t> Graph::indirect_transitive_closure(size_t start_id)
{
    std::unordered_set<size_t> visited;
    indirect_transitive_closure_helper(start_id, visited);
    return visited;
}

void Graph::indirect_transitive_closure_helper(size_t node_id, std::unordered_set<size_t>& visited)
{
    // Mark the current node as visited
    if (visited.find(node_id) == visited.end())
    {
        visited.insert(node_id);

        // Recur for all vertices pointing to this vertex
        Node *current_node = _first;
        while (current_node != nullptr)
        {
            Edge *current_edge = current_node->_first_edge;
            while (current_edge != nullptr)
            {
                if (current_edge->_target_id == node_id)
                {
                    indirect_transitive_closure_helper(current_node->_id, visited);
                }
                current_edge = current_edge->_next_edge;
            }
            current_node = current_node->_next_node;
        }
    }
}

void Graph::dijkstra_shortest_path(size_t start_id, size_t end_id)
{
    std::unordered_map<size_t, float> distances;
    std::unordered_map<size_t, size_t> previous;
    std::priority_queue<std::pair<float, size_t>, std::vector<std::pair<float, size_t> >, std::greater<std::pair<float, size_t> > > queue;

    for (Node *node = _first; node != nullptr; node = node->_next_node)
    {
        distances[node->_id] = (node->_id == start_id) ? 0 : std::numeric_limits<float>::infinity();
        queue.push(std::make_pair(distances[node->_id], node->_id));
    }

    while (!queue.empty())
    {
        size_t current_id = queue.top().second;
        queue.pop();

        Node *current_node = findNode(current_id);
        if (!current_node) continue;

        for (Edge *edge = current_node->_first_edge; edge != nullptr; edge = edge->_next_edge)
        {
            float alt_distance = distances[current_id] + edge->_weight;
            if (alt_distance < distances[edge->_target_id])
            {
                distances[edge->_target_id] = alt_distance;
                previous[edge->_target_id] = current_id;
                queue.push(std::make_pair(alt_distance, edge->_target_id));
            }
        }
    }

    // Print shortest path from start_id to end_id
    std::vector<size_t> path;
    for (size_t at = end_id; at != start_id; at = previous[at])
    {
        path.push_back(at);
    }
    path.push_back(start_id);  // Add start node
    std::reverse(path.begin(), path.end());

    // DFS-style output
    std::cout << "Dijkstra's Algorithm: Shortest path from " << start_id << " to " << end_id << ":\n";
    std::cout << "Path: ";
    for (size_t node : path)
    {
        std::cout << node << " ";
    }
    std::cout << "\nTotal Distance: " << distances[end_id] << "\n" << std::endl;
}


void Graph::floyd_warshall_shortest_path(size_t start_id, size_t end_id)
{
    size_t num_nodes = _number_of_nodes;
    std::vector<std::vector<float> > dist(num_nodes, std::vector<float>(num_nodes, std::numeric_limits<float>::infinity()));
    std::vector<std::vector<size_t> > next(num_nodes, std::vector<size_t>(num_nodes));

    // Initialize distances with edge weights and next nodes
    for (Node *node = _first; node != nullptr; node = node->_next_node)
    {
        dist[node->_id][node->_id] = 0;
        for (Edge *edge = node->_first_edge; edge != nullptr; edge = edge->_next_edge)
        {
            dist[node->_id][edge->_target_id] = edge->_weight;
            next[node->_id][edge->_target_id] = edge->_target_id;
        }
    }

    // Floyd-Warshall algorithm
    for (size_t k = 0; k < num_nodes; ++k)
    {
        for (size_t i = 0; i < num_nodes; ++i)
        {
            for (size_t j = 0; j < num_nodes; ++j)
            {
                if (dist[i][k] < std::numeric_limits<float>::infinity() && dist[k][j] < std::numeric_limits<float>::infinity() &&
                    dist[i][j] > dist[i][k] + dist[k][j])
                {
                    dist[i][j] = dist[i][k] + dist[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }

    // Reconstruct shortest path from start_id to end_id
    std::vector<size_t> path;
    if (dist[start_id][end_id] == std::numeric_limits<float>::infinity())
    {
        std::cout << "Floyd-Warshall Algorithm: No path exists from " << start_id << " to " << end_id << std::endl;
        return;  // No path exists
    }
    for (size_t at = start_id; at != end_id; at = next[at][end_id])
    {
        path.push_back(at);
    }
    path.push_back(end_id);

    // DFS-style output
    std::cout << "Floyd-Warshall Algorithm: Shortest path from " << start_id << " to " << end_id << ":\n";
    std::cout << "Path: ";
    for (size_t node : path)
    {
        std::cout << node << " ";
    }
    std::cout << "\nTotal Distance: " << dist[start_id][end_id] << "\n" << std::endl;
}

Graph Graph::prim_minimum_spanning_tree(const std::set<size_t>& subset_of_vertices)
{
    // Initialize a new graph for the minimum spanning tree
    Graph minimum_spanning_tree(_directed, _weighted_edges, _weighted_nodes);

    // Define a comparator for the priority queue that uses edges
    struct Compare
    {
        bool operator()(const std::pair<float, std::pair<size_t, size_t> >& a, const std::pair<float, std::pair<size_t, size_t> >& b) const
        {
            return a.first > b.first;
        }
    };

    // Initialize a priority queue to store edges based on their weights
    std::priority_queue<std::pair<float, std::pair<size_t, size_t> >, std::vector<std::pair<float, std::pair<size_t, size_t> > >, Compare> pq;

    // Initialize a set to keep track of visited vertices
    std::unordered_set<size_t> visited;

    // Add the starting vertex to the visited set
    if (!subset_of_vertices.empty())
        visited.insert(*subset_of_vertices.begin());

    // Add all edges incident to the starting vertex to the priority queue
    for (std::set<size_t>::const_iterator it = subset_of_vertices.begin(); it != subset_of_vertices.end(); ++it)
    {
        size_t vertex = *it;
        Node  *node   = findNode(vertex);
        if (node)
        {
            Edge *edge = node->_first_edge;
            while (edge != nullptr)
            {
                pq.push(std::make_pair(edge->_weight, std::make_pair(vertex, edge->_target_id)));
                edge = edge->_next_edge;
            }
        }
    }

    // Loop until all vertices are visited or the priority queue is empty
    while (!pq.empty() && visited.size() < subset_of_vertices.size())
    {
        // Get the minimum weighted edge from the priority queue
        std::pair<float, std::pair<size_t, size_t> > top    = pq.top();
        float                                       weight = top.first;
        size_t                                      source = top.second.first;
        size_t                                      target = top.second.second;
        pq.pop();

        // Check if the target vertex has not been visited yet
        if (visited.find(target) == visited.end())
        {
            // Add the edge to the minimum spanning tree
            minimum_spanning_tree.add_node(source);
            minimum_spanning_tree.add_node(target);
            minimum_spanning_tree.add_edge(source, target, weight);

            // Mark the target vertex as visited
            visited.insert(target);

            // Add all edges incident to the target vertex to the priority queue
            Node *node = findNode(target);
            if (node)
            {
                Edge *edge = node->_first_edge;
                while (edge != nullptr)
                {
                    pq.push(std::make_pair(edge->_weight, std::make_pair(target, edge->_target_id)));
                    edge = edge->_next_edge;
                }
            }
        }
    }

    return minimum_spanning_tree;
}

class DisjointSets
{
public:
    DisjointSets(size_t size) : _parent(size), _rank(size, 0)
    {
        for (size_t i = 0; i < size; ++i)
            _parent[i] = i;
    }

    size_t find(size_t u)
    {
        if (_parent[u] != u)
            _parent[u] = find(_parent[u]);
        return _parent[u];
    }

    void merge(size_t u, size_t v)
    {
        size_t u_root = find(u);
        size_t v_root = find(v);

        if (u_root == v_root)
            return;

        if (_rank[u_root] < _rank[v_root])
            _parent[u_root] = v_root;
        else if (_rank[u_root] > _rank[v_root])
            _parent[v_root] = u_root;
        else
        {
            _parent[v_root] = u_root;
            _rank[u_root]++;
        }
    }

private:
    std::vector<size_t> _parent;
    std::vector<size_t> _rank;
};

Graph Graph::kruskalMST(const std::set<size_t>& subset_of_vertices)
{
    // Create a graph for the MST
    Graph MST(_directed, _weighted_edges, _weighted_nodes);

    // Container to store all edges (weight, u, v)
    std::vector<std::tuple<float, size_t, size_t> > edges;

    // Traverse all nodes and their edges
    Node* curr_node = _first;
    while (curr_node != nullptr)
    {
        Edge* curr_edge = curr_node->_first_edge;
        while (curr_edge != nullptr)
        {
            size_t u = curr_node->_id;
            size_t v = curr_edge->_target_id;

            // Consider only edges where both vertices are in the subset
            if (subset_of_vertices.count(u) && subset_of_vertices.count(v))
            {
                edges.emplace_back(curr_edge->_weight, u, v);
            }

            curr_edge = curr_edge->_next_edge;
        }
        curr_node = curr_node->_next_node;
    }

    // Sort edges by weight
    std::sort(edges.begin(), edges.end());

    // Initialize the DisjointSets structure
    DisjointSets ds(_number_of_nodes);
    float total_weight = 0.0;

    // Process each edge
    for (const auto& edge : edges)
    {
        float weight;
        size_t u, v;
        std::tie(weight, u, v) = edge;

        // If u and v are in different sets, add the edge to MST
        if (ds.find(u) != ds.find(v))
        {
            MST.add_edge(u, v, weight);
            ds.merge(u, v);
            total_weight += weight;

            std::cout << u << " -- " << v << " == " << weight << std::endl;
        }
    }

    std::cout << "Minimum Cost Spanning Tree: " << total_weight << std::endl;

    return MST;
}


void Graph::calculateGraphMetrics(Graph& graph)
{
    size_t                     radius   = std::numeric_limits<size_t>::max();
    size_t                     diameter = 0;
    std::unordered_set<size_t> center;
    std::unordered_set<size_t> periphery;

    for (Node *node = graph._first; node != nullptr; node = node->_next_node)
    {
        // Calculate eccentricity for each vertex
        std::unordered_map<size_t, float>                                                                                            distances;
        std::priority_queue<std::pair<float, size_t>, std::vector<std::pair<float, size_t> >, std::greater<std::pair<float, size_t> > > queue;

        for (Node *n = graph._first; n != nullptr; n = n->_next_node)
        {
            if (n->_id == node->_id)
            {
                distances[n->_id] = 0;
            }
            else
            {
                distances[n->_id] = std::numeric_limits<float>::infinity();
            }
            queue.push(std::make_pair(distances[n->_id], n->_id));
        }

        while (!queue.empty())
        {
            size_t current_id = queue.top().second;
            queue.pop();

            Node *current_node = graph.findNode(current_id);
            if (current_node == nullptr)
                continue;

            for (Edge *edge = current_node->_first_edge; edge != nullptr; edge = edge->_next_edge)
            {
                float alt_distance = distances[current_id] + edge->_weight;
                if (alt_distance < distances[edge->_target_id])
                {
                    distances[edge->_target_id] = alt_distance;
                    queue.push(std::make_pair(alt_distance, edge->_target_id));
                }
            }
        }

        // Find eccentricity for current vertex
        size_t max_distance = 0;
        for (const auto& pair : distances)
        {
            if (pair.second > max_distance && pair.second != std::numeric_limits<float>::infinity())
            {
                max_distance = pair.second;
            }
        }

        // Update radius and diameter
        if (max_distance < radius)
        {
            radius = max_distance;
            center.clear();
            center.insert(node->_id);
        }
        else if (max_distance == radius)
        {
            center.insert(node->_id);
        }

        if (max_distance > diameter)
        {
            diameter = max_distance;
            periphery.clear();
            periphery.insert(node->_id);
        }
        else if (max_distance == diameter)
        {
            periphery.insert(node->_id);
        }
    }

    // Print results
    std::cout << "Radius: " << radius << std::endl;
    std::cout << "Diameter: " << diameter << std::endl;
    std::cout << "Center: ";
    for (size_t id : center)
    {
        std::cout << id << " ";
    }
    std::cout << std::endl;
    std::cout << "Periphery: ";
    for (size_t id : periphery)
    {
        std::cout << id << " ";
    }
    std::cout << std::endl;
}

void Graph::greedy_partition(size_t p) {
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::set<size_t>> subgraphs(p);
    std::unordered_set<size_t> unvisited;

    // Initialize all nodes as unvisited
    for (Node* node = _first; node != nullptr; node = node->_next_node) {
        unvisited.insert(node->_id);
    }

    // Step 1: Initial greedy assignment to subgraphs
    for (size_t i = 0; i < p; i++) {
        if (unvisited.empty()) break;

        size_t node_id = *unvisited.begin();  // Get the first available node
        subgraphs[i].insert(node_id);
        unvisited.erase(node_id);

        // Try to add neighbors to this subgraph
        for (auto edge = findNode(node_id)->_first_edge; edge != nullptr; edge = edge->_next_edge) {
            size_t target_id = edge->_target_id;
            if (unvisited.find(target_id) != unvisited.end() && connected(node_id, target_id)) {
                subgraphs[i].insert(target_id);
                unvisited.erase(target_id);
            }
            // Stop when subgraph has at least two nodes
            if (subgraphs[i].size() >= 2) break;
        }

        // Ensure subgraph has at least 2 vertices
        while (subgraphs[i].size() < 2 && !unvisited.empty()) {
            size_t extra_node_id = *unvisited.begin();
            subgraphs[i].insert(extra_node_id);
            unvisited.erase(extra_node_id);
        }
    }

    // Step 2: Assign remaining unvisited nodes
    while (!unvisited.empty()) {
        for (auto& subgraph : subgraphs) {
            if (unvisited.empty()) break;

            size_t node_id = *unvisited.begin();
            subgraph.insert(node_id);
            unvisited.erase(node_id);
        }
    }

    // Final check to ensure all subgraphs have at least 2 nodes
    for (auto& subgraph : subgraphs) {
        while (subgraph.size() < 2 && !unvisited.empty()) {
            size_t node_id = *unvisited.begin();
            subgraph.insert(node_id);
            unvisited.erase(node_id);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;

    float total_gap = calculate_total_gap(subgraphs);

    std::cout << "Greedy Partition Quality: " << total_gap << "\n";
    std::cout << "Greedy Partition Time: " << duration.count() << " seconds\n";

    //print_partition(subgraphs);
    //print_partition_weight(subgraphs);
}

// GRASP partition method
void Graph::grasp_partition(size_t p, size_t iterations) {
    auto start = std::chrono::high_resolution_clock::now();

    std::srand(std::time(nullptr));
    std::vector<std::set<size_t>> best_partition;
    float best_gap = std::numeric_limits<float>::max();

    for (size_t i = 0; i < iterations; i++) {
        std::vector<std::set<size_t>> current_partition = randomized_greedy_partition(p);
        current_partition = local_search_refine(current_partition);
        float current_gap = calculate_total_gap(current_partition);

        if (current_gap < best_gap) {
            best_gap = current_gap;
            best_partition = current_partition;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;

    std::cout << "GRASP Partition Quality: " << best_gap << "\n";
    std::cout << "GRASP Partition Time: " << duration.count() << " seconds\n";

    //print_partition(best_partition);
    //print_partition_weight(best_partition);
}

// Reactive GRASP
void Graph::reactive_grasp_partition(size_t p, size_t iterations) {
    auto start = std::chrono::high_resolution_clock::now();

    std::srand(std::time(nullptr));
    std::vector<std::set<size_t>> best_partition;
    float best_gap = std::numeric_limits<float>::max();
    std::vector<float> randomness_levels = {0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<float> performance(randomness_levels.size(), 0);

    for (size_t i = 0; i < iterations; i++) {
        float randomness = adapt_randomness(performance);
        std::vector<std::set<size_t>> current_partition = randomized_grasp(randomness, p);
        current_partition = local_search_refine(current_partition);
        float current_gap = calculate_total_gap(current_partition);

        if (current_gap < best_gap) {
            best_gap = current_gap;
            best_partition = current_partition;
        }

        update_performance(performance, randomness, current_gap);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;

    std::cout << "Reactive GRASP Partition Quality: " << best_gap << "\n";
    std::cout << "Reactive GRASP Partition Time: " << duration.count() << " seconds\n";

    //print_partition(best_partition);
    //print_partition_weight(best_partition);
}

/* Randomized Greedy Partition for GRASP */
std::vector<std::set<size_t>> Graph::randomized_greedy_partition(size_t p) {
    std::vector<std::set<size_t>> subgraphs(p);
    std::unordered_set<size_t> unvisited;

    // Initialize all nodes as unvisited
    for (Node* node = _first; node != nullptr; node = node->_next_node) {
        unvisited.insert(node->_id);
    }

    // Randomized start for each subgraph
    for (size_t i = 0; i < p; i++) {
        if (unvisited.empty()) break;

        size_t node_id = get_random_node(unvisited);
        subgraphs[i].insert(node_id);
        unvisited.erase(node_id);

        // Randomly add neighbors to the subgraph
        for (auto edge = findNode(node_id)->_first_edge; edge != nullptr; edge = edge->_next_edge) {
            if (unvisited.find(edge->_target_id) != unvisited.end()) {
                if ((std::rand() / static_cast<float>(RAND_MAX)) < 0.5) { // Adjusted for testing
                    subgraphs[i].insert(edge->_target_id);
                    unvisited.erase(edge->_target_id);
                }
            }
            // Stop when subgraph has at least two nodes
            if (subgraphs[i].size() >= 2) break;
        }

        // Ensure subgraph has at least 2 vertices
        while (subgraphs[i].size() < 2 && !unvisited.empty()) {
            size_t extra_node_id = *unvisited.begin();
            subgraphs[i].insert(extra_node_id);
            unvisited.erase(extra_node_id);
        }
    }

    // Assign remaining unvisited nodes randomly to subgraphs
    while (!unvisited.empty()) {
        for (auto& subgraph : subgraphs) {
            if (unvisited.empty()) break;

            size_t node_id = get_random_node(unvisited);
            subgraph.insert(node_id);
            unvisited.erase(node_id);
        }
    }

    return subgraphs;
}

/* Local Search Refinement */
std::vector<std::set<size_t>> Graph::local_search_refine(const std::vector<std::set<size_t>>& partition) {
    std::vector<std::set<size_t>> refined_partition = partition;

    // Simple local search: try moving nodes between subgraphs
    for (size_t i = 0; i < refined_partition.size(); i++) {
        std::set<size_t> current_subgraph = refined_partition[i];
        for (auto it = current_subgraph.begin(); it != current_subgraph.end(); ++it) {
            size_t node_id = *it;
            float best_gap = calculate_total_gap(refined_partition);

            // Try moving node to a different subgraph
            for (size_t j = 0; j < refined_partition.size(); j++) {
                if (i == j) continue;

                refined_partition[i].erase(node_id);
                refined_partition[j].insert(node_id);

                float new_gap = calculate_total_gap(refined_partition);
                if (new_gap < best_gap) {
                    best_gap = new_gap;
                } else {
                    // Revert move
                    refined_partition[j].erase(node_id);
                    refined_partition[i].insert(node_id);
                }
            }
        }
    }

    return refined_partition;
}

/* Adapt Randomness Based on Performance */
float Graph::adapt_randomness(const std::vector<float>& performance) {
    // Choose randomness level based on performance
    float best_performance = *std::min_element(performance.begin(), performance.end());
    size_t best_index = std::distance(performance.begin(), std::find(performance.begin(), performance.end(), best_performance));
    return (0.1f + (best_index * 0.2f));  // Returns randomness level based on index
}

/* Update Performance of Randomness Levels */
void Graph::update_performance(std::vector<float>& performance, float randomness, float gap) {
    size_t index = (randomness - 0.1f) / 0.2f;
    performance[index] += gap;
}

/* Helper Functions */
float Graph::calculate_total_gap(const std::vector<std::set<size_t>>& partition) {
    float total_gap = 0;

    for (const auto& subgraph : partition) {
        float min_weight = std::numeric_limits<float>::max();
        float max_weight = std::numeric_limits<float>::lowest();

        for (auto node_id : subgraph) {
            Node* node = findNode(node_id);
            if (node) {
                min_weight = std::min(min_weight, node->_weight);
                max_weight = std::max(max_weight, node->_weight);
            }
        }

        total_gap += (max_weight - min_weight);
    }

    return total_gap;
}

size_t Graph::get_random_node(const std::unordered_set<size_t>& unvisited) {
    auto it = unvisited.begin();
    std::advance(it, std::rand() % unvisited.size());
    return *it;
}

void Graph::print_partition_weight(const std::vector<std::set<size_t>>& partition) {
    for (size_t i = 0; i < partition.size(); i++) {
        std::cout << "Subgraph " << i + 1 << ": ";
        for (auto node : partition[i]) {
            std::cout << findNode(node)->_weight << " ";
        }
        std::cout << std::endl;
    }
}
void Graph::print_partition(const std::vector<std::set<size_t>>& partition) {
    for (size_t i = 0; i < partition.size(); i++) {
        std::cout << "Subgraph " << i + 1 << ": ";
        for (auto node_id : partition[i]) {
            Node* node = findNode(node_id);
            if (node) {
                std::cout << node->_id << " ";
            }
        }
        std::cout << std::endl;
    }
}

std::vector<std::set<size_t>> Graph::randomized_grasp(float randomness, size_t p) {
    std::vector<std::set<size_t>> subgraphs(p);
    std::unordered_set<size_t> unvisited;

    // Initialize all nodes as unvisited
    for (Node* node = _first; node != nullptr; node = node->_next_node) {
        unvisited.insert(node->_id);
    }

    // Randomized start for each subgraph
    for (size_t i = 0; i < p; i++) {
        if (unvisited.empty()) break;

        size_t node_id = get_random_node(unvisited);
        subgraphs[i].insert(node_id);
        unvisited.erase(node_id);

        // Randomly add neighbors to minimize the gap
        for (auto edge = findNode(node_id)->_first_edge; edge != nullptr; edge = edge->_next_edge) {
            if (unvisited.find(edge->_target_id) != unvisited.end()) {
                // Use randomness to decide whether to add the neighbor
                if ((std::rand() / static_cast<float>(RAND_MAX)) < randomness) {
                    subgraphs[i].insert(edge->_target_id);
                    unvisited.erase(edge->_target_id);
                }
            }
            // Stop when subgraph has at least two nodes
            if (subgraphs[i].size() >= 2) break;
        }

        // Ensure subgraph has at least 2 vertices
        while (subgraphs[i].size() < 2 && !unvisited.empty()) {
            size_t extra_node_id = *unvisited.begin();
            subgraphs[i].insert(extra_node_id);
            unvisited.erase(extra_node_id);
        }
    }

    // Assign remaining unvisited nodes randomly to subgraphs
    while (!unvisited.empty()) {
        for (auto& subgraph : subgraphs) {
            if (unvisited.empty()) break;

            size_t node_id = get_random_node(unvisited);
            subgraph.insert(node_id);
            unvisited.erase(node_id);
        }
    }

    return subgraphs;
}