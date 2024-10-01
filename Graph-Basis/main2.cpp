#include "include/Graph.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

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
// Function to read vertices and weights
void read_input_file(const std::string& inputFile, std::set<size_t>& V, 
                     std::map<size_t, float>& w) {
    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "Error opening file '" << inputFile << "'\n";
        return;
    }

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
}

// New function to read the edges
void read_edges(const std::string& inputFile, std::set<std::pair<size_t, size_t>>& E) {
    std::ifstream file(inputFile);
    if (!file) {
        std::cerr << "Error opening file '" << inputFile << "'\n";
        return;
    }

    std::string line;
    bool edgeSectionFound = false;

    // Find the "set E :=" line and start reading edges
    while (std::getline(file, line)) {
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

#include "include/Graph.hpp"
#include <chrono>
#include <filesystem>  // For directory iteration
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

// Function to capture output of graph operations
template <typename Func>
std::string capture_output(Func func) {
    std::ostringstream oss;
    func();
    return oss.str();
}

// Function to save output to a file
void save_to_file(const std::string& filename, const std::string& output) {
    std::ofstream outFile(filename, std::ios::app);
    if (!outFile) {
        std::cerr << "Error opening output file '" << filename << "'\n";
        return;
    }
    outFile << output;
}

// Function to process a single file
void process_file(const std::string& inputFile, const std::string& resultFile) {
    size_t p = get_p_value(inputFile);  // Get the value of p
    std::set<size_t> V;  // vertices set
    std::map<size_t, float> w;  // weight vector
    std::set<std::pair<size_t, size_t>> E;  // edges set

    // Read the input file for vertices and weights
    read_input_file(inputFile, V, w);

    // Read the edges separately
    read_edges(inputFile, E);

    // Initialize graph
    bool directed = false;  // Adjust if directed graph is needed
    bool weightedEdges = true;
    bool weightedNodes = true;
    Graph graph(directed, weightedEdges, weightedNodes);

    // Add vertices and their weights to the graph
    for (const auto& vertex : V) {
        graph.add_node(vertex, w[vertex]);
    }

    // Add edges to the graph
    for (const auto& edge : E) {
        graph.add_edge(edge.first, edge.second);
    }

    // Run each algorithm 30 times and save output to the result file
    for (int i = 0; i < 30; ++i) {
        // Greedy Algorithm
        auto start = std::chrono::high_resolution_clock::now();
        std::string greedyOutput = capture_output([&]() { graph.greedy_partition(p); });
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        save_to_file(resultFile, "Greedy Partition Quality: " + greedyOutput + "\n");
        save_to_file(resultFile, "Greedy Partition Time: " + std::to_string(elapsed.count()) + " seconds\n\n");

        // GRASP Algorithm
        size_t iterations = 100;  // Example number of iterations
        start = std::chrono::high_resolution_clock::now();
        std::string graspOutput = capture_output([&]() { graph.grasp_partition(p, iterations); });
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        save_to_file(resultFile, "GRASP Partition Quality: " + graspOutput + "\n");
        save_to_file(resultFile, "GRASP Partition Time: " + std::to_string(elapsed.count()) + " seconds\n\n");

        // Reactive GRASP Algorithm
        start = std::chrono::high_resolution_clock::now();
        std::string reactiveGraspOutput = capture_output([&]() { graph.reactive_grasp_partition(p, iterations); });
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        save_to_file(resultFile, "Reactive GRASP Partition Quality: " + reactiveGraspOutput + "\n");
        save_to_file(resultFile, "Reactive GRASP Partition Time: " + std::to_string(elapsed.count()) + " seconds\n\n");
    }
}

int main() {
    std::string inputFolder = "instances_example";
    std::string outputFolder = "results";

    // Create results folder if it doesn't exist
    std::__fs::filesystem::create_directory(outputFolder);

    // Iterate over each file in the instances_example folder
    for (const auto& entry : std::__fs::filesystem::directory_iterator(inputFolder)) {
        std::string inputFile = entry.path().string();  // Full path of the input file
        std::string filename = entry.path().filename().string();  // Just the filename
        std::string resultFile = outputFolder + "/result_" + filename;  // Output file path

        std::cout << "Processing file: " << filename << std::endl;

        // Process the file and save the results
        process_file(inputFile, resultFile);
    }

    return 0;
}
