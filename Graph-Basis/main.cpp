#include "include/Graph.hpp"
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

// Function to ask user to save output to a file
void ask_to_save(const std::string& output)
{
    char choice;
    std::cout << "Do you want to save the output to a file? (y/n): ";
    std::cin >> choice;
    if (choice == 'y' || choice == 'Y')
    {
        std::string filename;
        std::cout << "Enter the filename: ";
        std::cin >> filename;
        std::ofstream outFile(filename);
        if (outFile)
        {
            outFile << output;
            outFile.close();
            std::cout << "Output saved to " << filename << " successfully.\n";
        }
        else
        {
            std::cerr << "Error opening file '" << filename << "'\n";
        }
    }
}

void save_to_file(const std::string& output) {
    std::ofstream outFile("results.txt", std::ios::app);  // Open file in append mode
    if (outFile) {
        outFile << output << "\n";  // Write output and add a newline
        outFile.close();
        std::cout << "Output saved to results.txt successfully.\n";
    } else {
        std::cerr << "Error opening file 'results.txt'\n";
    }
}

void read_edges(const std::string& inputFile, std::set<std::pair<size_t, size_t>>& E)
{
    std::ifstream file(inputFile);
    if (!file)
    {
        std::cerr << "Error opening file '" << inputFile << "'\n";
        return;
    }

    std::string line;
    bool edgeSectionFound = false;

    // Find the "set E :=" line and start reading edges
    while (std::getline(file, line))
    {
        if (line.find("set E :=") != std::string::npos) {
            edgeSectionFound = true;
            break;
        }
    }

    // Parse edges if the section was found
    if (edgeSectionFound) {
        std::string edgeLine;
        while (std::getline(file, edgeLine) && edgeLine.find(";") == std::string::npos) {
            std::istringstream edgeStream(edgeLine);
            char ch;
            size_t v1, v2;
            // Read edges in the format (v1, v2)
            while (edgeStream >> ch >> v1 >> ch >> v2 >> ch) {
                E.insert({v1, v2});
            }
        }
    } else {
        std::cerr << "Error: Edge section not found in the input file.\n";
    }
}

size_t get_p_value(const std::string& inputFile) {
    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "Error opening file '" << inputFile << "'\n";
        return 0;  // or an appropriate value indicating error
    }

    std::string line;
    size_t p = 0;  // Initialize p

    // Read through the file to find the parameter p
    while (std::getline(file, line)) {
        if (line.find("param p :=") != std::string::npos) {
            std::istringstream iss(line);
            std::string key;

            // Skip the first three tokens
            iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
            iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
            iss.ignore(std::numeric_limits<std::streamsize>::max(), ' ');

            // Extract the value of p
            iss >> p;

            // Skip any remaining tokens
            iss.ignore(std::numeric_limits<std::streamsize>::max());

            break;  // Stop reading once p is found
        }
    }

    return p;  // Return the value of p found in the file
}
// Capture output of a function into a string
std::string capture_output(std::function<void()> func)
{
    std::ostringstream oss;
    std::streambuf    *orig_buf = std::cout.rdbuf(oss.rdbuf());
    func();
    std::cout.rdbuf(orig_buf);  // Restore original buffer
    return oss.str();
}

int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <Op_Directed> <Op_WeightedEdges> <Op_WeightedNodes>\n";
        return 1;
    }

    std::string inputFile     = argv[1];
    std::string outputFile    = argv[2];
    bool        directed      = (std::stoi(argv[3]) != 0);
    bool        weightedEdges = (std::stoi(argv[4]) != 0);
    bool        weightedNodes = (std::stoi(argv[5]) != 0);

    std::ifstream file(inputFile);
    if (!file)
    {
        std::cerr << "Error opening file '" << inputFile << "'\n";
        return 1;
    }

    try {
    size_t                              p;  // number of subgraphs
    std::set<size_t>                    V;  // vertices set
    std::map<size_t, float>             w;  // weight vector
    std::set<std::pair<size_t, size_t>> E;  // edges set

    size_t v1, v2;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        
        // Parse vertices (V)
        if (line.find("set V :=") != std::string::npos) {
            size_t vertex;
            std::getline(file, line);  // Move to the next line with vertices
            std::istringstream vStream(line);
            while (vStream >> vertex) {
                V.insert(vertex);
            }
        }
        // Parse weights (w)
        else if (line.find("param w :=") != std::string::npos) {
            size_t vertex;
            float  weight;
            while (file >> vertex >> weight) {
                if (vertex == 0)
                    break;  // Stop when 0 is reached
                w[vertex] = weight;
            }
        }
    }

    // Initialize graph based on parsed data
    Graph graph(directed, weightedEdges, weightedNodes);

    // Add vertices and their weights to the graph
    for (const auto& vertex : V) {
        graph.add_node(vertex, w[vertex]);
    }

    p = get_p_value(inputFile);
    read_edges(inputFile, E);
    // Add edges to the graph
    for (const auto& edge : E) {
        graph.add_edge(edge.first, edge.second);
    }

        // Main menu for operations
        int option;
        do
        {
            std::cout << "Menu:\n"
                         "1. Graph Operations\n"
                         "2. Graph Functions\n"
                         "3. Algorithm Tests\n"
                         "0. Exit\n"
                         "Enter your choice: ";
            std::cin >> option;

            switch (option)
            {
                case 1:
                {
                    // Graph Operations submenu
                    int opOption;
                    do
                    {
                        std::cout << "Graph Operations:\n"
                                     "1. Print the graph\n"
                                     "2. Add a node\n"
                                     "3. Add an edge\n"
                                     "4. Remove a node\n"
                                     "5. Remove an edge\n"
                                     "6. Check if nodes are connected\n"
                                     "7. Export graph to DOT file\n"
                                     "0. Back to main menu\n"
                                     "Enter your choice: ";
                        std::cin >> opOption;

                        switch (opOption)
                        {
                            case 1:
                            {
                                std::string output = capture_output([&]() { graph.print_graph(); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 2:
                            {
                                size_t nodeId;
                                float  weight;
                                std::cout << "Enter node id and weight: ";
                                std::cin >> nodeId >> weight;
                                graph.add_node(nodeId, weight);
                                break;
                            }
                            case 3:
                            {
                                size_t sourceId, targetId;
                                float  weight;
                                std::cout << "Enter source node id, target node id, and weight: ";
                                std::cin >> sourceId >> targetId >> weight;
                                graph.add_edge(sourceId, targetId, weight);
                                break;
                            }
                            case 4:
                            {
                                size_t nodeId;
                                std::cout << "Enter node id to remove: ";
                                std::cin >> nodeId;
                                graph.remove_node(nodeId);
                                break;
                            }
                            case 5:
                            {
                                size_t sourceId, targetId;
                                std::cout << "Enter source node id and target node id to remove edge: ";
                                std::cin >> sourceId >> targetId;
                                graph.remove_edge(sourceId, targetId);
                                break;
                            }
                            case 6:
                            {
                                size_t nodeId1, nodeId2;
                                std::cout << "Enter two node ids to check connectivity: ";
                                std::cin >> nodeId1 >> nodeId2;
                                int         connected = graph.connected(nodeId1, nodeId2);
                                std::string output    = "Nodes " + std::to_string(nodeId1) + " and " + std::to_string(nodeId2) + " are " +
                                                     (connected ? "connected" : "not connected") + "\n";
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 7:  // New case for exporting to DOT file
                            {
                                std::ofstream dotFile(outputFile);
                                if (!dotFile)
                                {
                                    std::cerr << "Error opening output file '" << outputFile << "'\n";
                                    break;
                                }
                                graph.print_dot(dotFile);
                                dotFile.close();
                                std::cout << "Graph exported to " << outputFile << " successfully.\n";
                                break;
                            }
                            case 0:
                                std::cout << "Returning to main menu...\n";
                                break;
                            default:
                                std::cout << "Invalid choice. Please try again.\n";
                        }

                    } while (opOption != 0);
                    break;
                }
                case 2:
                {
                    // Graph Functions submenu
                    int fnOption;
                    do
                    {
                        std::cout << "Graph Functions:\n"
                                     "1. Compute direct transitive closure of a vertex\n"
                                     "2. Compute indirect transitive closure of a vertex\n"
                                     "3. Dijkstra's algorithm\n"
                                     "4. Floyd-Warshall algorithm\n"
                                     "5. Minimum Spanning Tree using Prim's algorithm\n"
                                     "6. Minimum Spanning Tree using Kruskal's algorithm\n"
                                     "7. DFS with back edges\n"
                                     "8. Graph Metrics\n"
                                     "0. Back to main menu\n"
                                     "Enter your choice: ";
                        std::cin >> fnOption;

                        switch (fnOption)
                        {
                            case 1:
                                if (!graph.isDirected())
                                {
                                    std::cout << "Graph is not directed." << std::endl;
                                    break;
                                }
                                {
                                    size_t vertexId;
                                    std::cout << "Enter a vertex id to compute direct transitive closure: ";
                                    std::cin >> vertexId;
                                    std::string output = capture_output([&]() { graph.direct_transitive_closure(vertexId); });
                                    std::cout << output;
                                    ask_to_save(output);
                                    break;
                                }
                            case 2:
                                if (!graph.isDirected())
                                {
                                    std::cout << "Graph is not directed." << std::endl;
                                    break;
                                }
                                {
                                    size_t vertexId;
                                    std::cout << "Enter a vertex id to compute indirect transitive closure: ";
                                    std::cin >> vertexId;
                                    std::string output = capture_output([&]() { graph.indirect_transitive_closure(vertexId); });
                                    std::cout << output;
                                    ask_to_save(output);
                                    break;
                                }
                            case 3:
                            {
                                size_t startId, endId;
                                std::cout << "Enter start and end node ids for Dijkstra's algorithm: ";
                                std::cin >> startId >> endId;
                                std::string output = capture_output([&]() { graph.dijkstra_shortest_path(startId, endId); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 4:
                            {
                                size_t startId, endId;
                                std::cout << "Enter start and end node ids for Floyd-Warshall algorithm: ";
                                std::cin >> startId >> endId;
                                std::string output = capture_output([&]() { graph.floyd_warshall_shortest_path(startId, endId); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 5:
                            {
                                size_t numVertices;
                                std::cout << "Enter the number of vertices: ";
                                std::cin >> numVertices;
                                std::set<size_t> vertexSet;
                                for (size_t i = 0; i < numVertices; ++i)
                                {
                                    size_t vertex;
                                    std::cout << "Enter vertex " << i + 1 << ": ";
                                    std::cin >> vertex;
                                    vertexSet.insert(vertex);
                                }

                                Graph       mst    = graph.prim_minimum_spanning_tree(vertexSet);
                                std::string output = capture_output([&]() { mst.print_graph(); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 6:
                            {
                                size_t numVertices;
                                std::cout << "Enter the number of vertices: ";
                                std::cin >> numVertices;
                                std::set<size_t> vertexSet;
                                for (size_t i = 0; i < numVertices; ++i)
                                {
                                    size_t vertex;
                                    std::cout << "Enter vertex " << i + 1 << ": ";
                                    std::cin >> vertex;
                                    vertexSet.insert(vertex);
                                }

                                Graph       mst    = graph.kruskalMST(vertexSet);
                                std::string output = capture_output([&]() { mst.print_graph(); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 7:
                            {
                                size_t startVertex;
                                std::cout << "Enter the starting vertex ID: ";
                                std::cin >> startVertex;

                                std::string output = capture_output([&]() { graph.dfs(startVertex); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 8:
                            {
                                std::string output = capture_output([&]() { graph.calculateGraphMetrics(graph); });
                                std::cout << output;
                                ask_to_save(output);
                                break;
                            }
                            case 0:
                                std::cout << "Returning to main menu...\n";
                                break;
                            default:
                                std::cout << "Invalid choice. Please try again.\n";
                        }

                    } while (fnOption != 0);
                    break;
                }
                case 3:
                {
                    // Algorithm Tests submenu
                    int algoOption;
                    do
                    {
                        std::cout << "Algorithm Tests:\n"
                                     "1. Test Greedy Algorithm\n"
                                     "2. Test GRASP Algorithm\n"
                                     "3. Test Reactive GRASP Algorithm\n"
                                     "0. Back to main menu\n"
                                     "Enter your choice: ";
                        std::cin >> algoOption;

                        switch (algoOption)
                        {
                            case 1:
                            {  // Greedy Algorithm
                                /*size_t p;
                                std::cout << "Enter the number of partitions (p): ";
                                std::cin >> p;*/
                                int x = 0;
                                while (x < 30)
                                {
                                    std::string output = capture_output([&]() { graph.greedy_partition(p); });
                                    std::cout << output;
                                    save_to_file(output);
                                    x++;
                                }
                                //ask_to_save(output);
                                break;
                            }
                            case 2:
                            {  // GRASP Algorithm
                                /*size_t p, iterations;
                                std::cout << "Enter the number of partitions (p): ";
                                std::cin >> p;
                                std::cout << "Enter the number of iterations: ";
                                std::cin >> iterations;*/

                                int x = 0;
                                while (x < 30)
                                {
                                    std::string output = capture_output([&]() { graph.grasp_partition(p, 500); });
                                    std::cout << output;
                                    save_to_file(output);
                                    x++;
                                }
                                //ask_to_save(output);
                                break;
                            }
                            case 3:
                            {  // Reactive GRASP Algorithm
                                /*size_t p, iterations;
                                std::cout << "Enter the number of partitions (p): ";
                                std::cin >> p;
                                std::cout << "Enter the number of iterations: ";
                                std::cin >> iterations;*/

                                int x = 0;
                                while (x < 30)
                                {
                                    std::string output = capture_output([&]() { graph.reactive_grasp_partition(p, 500); });
                                    std::cout << output;
                                    save_to_file(output);
                                    x++;
                                }
                                //ask_to_save(output);
                                break;
                            }
                            case 0:
                                std::cout << "Returning to main menu...\n";
                                break;
                            default:
                                std::cout << "Invalid choice. Please try again.\n";
                        }
                    } while (algoOption != 0);
                    break;
                }
                case 0:
                    std::cout << "Exiting program...\n";
                    break;
                default:
                    std::cout << "Invalid choice. Please try again.\n";
            }
        } while (option != 0);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        file.close();
        return 1;
    }

    file.close();
    return 0;
}
