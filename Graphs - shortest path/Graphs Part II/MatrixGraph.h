#ifndef __MATRIXGRAPH_H__
#define __MATRIXGRAPH_H__
#define ODD 1
#define EVEN 2
#include <vector>
#include <list>
#include <iostream>
#include <queue>
#include <stack>
#include <bits/stdc++.h>
#include <algorithm>
#define INF 0
#define STARTING_EDGES_CAPACITY 2

class Edge {
    public:
        int src{0};
        int dest{0};
        int distance{0};
        Edge(int src = 0, int dest = 0, int distance = 0) : src(src), dest(dest), distance(distance) {};
};

class Edges {
    private:
        std::list<Edge> my_edges;

    public:
        Edges() = default;

        void addEdge(int src, int dest, int distance) {
            my_edges.push_back(Edge(src, dest, distance));
        }

        void removeEdge(int src, int dest) {
            my_edges.remove_if([src,dest](Edge element) { if(element.dest == dest && element.src == src) return true; else return false;});
        }
    
        ~Edges() = default;

        int getSize() {
            return my_edges.size();
        }

        Edge getEdgeByIndex(int index) {
            std::list<Edge>::iterator it = my_edges.begin();
            std::advance(it, index);
            return *it;
        }
};

/**
 * Matrix implementation.
 */
class MatrixGraph {
private:
    std::vector<std::vector<int>> adjacency_matrix;
    int size{0};
    Edges my_edges;
    int max_distance{0};
    std::list< std::pair<int, int> > *adj;
public:
    // Constructor
    MatrixGraph(int size)
    {
        while(size <= 0) {
            std::cerr << "Invalid size! Insert a valid size!\n";
            std::cin >> size;
        }

        this->size = size;

        // Create vector size
        adjacency_matrix.resize(this->size);
        for (int i = 0; i < this->size; ++i) {
            adjacency_matrix.at(i).resize(this->size);
        }

        for(int i = 0; i < adjacency_matrix.size(); ++i) {
            for(int j = 0; j < adjacency_matrix.size(); ++j)
                if(i != j)
                    adjacency_matrix.at(i).at(j) = INF;
                else
                    adjacency_matrix.at(i).at(j) = 0;
        }

        adj = new std::list< std::pair<int, int> >[size];
    }

    // Destructor
    ~MatrixGraph() {
        delete adj;
        for(auto row : adjacency_matrix) {
            row.clear();
        }
        adjacency_matrix.clear();
    }

    /**
     * Adds an edge between two existing nodes.
     */
    void addEdge(int src, int dst, int distance = 0)
    {
        adjacency_matrix[src][dst] = distance;
        adjacency_matrix[dst][src] = distance;
        my_edges.addEdge(src, dst, distance);
        if(distance > max_distance)
            max_distance = distance;

        adj[src].push_back(std::make_pair(dst, distance));
        adj[dst].push_back(std::make_pair(src, distance));
    }

    /**
     * Get the weight between two nodes
    */
    int getWeightBetweenNodes(int src, int dest) {
        return adjacency_matrix.at(src).at(dest);
    }

    auto getAdj() {
        return adj;
    }

    // Get the whole edges
    Edges get_whole_edges() {
        return my_edges;
    }

    /*.
     */
    void removeEdge(int src, int dst)
    {
        adjacency_matrix[src][dst] = INF;
        adjacency_matrix[dst][src] = INF;
        my_edges.removeEdge(src, dst);

        int weight = getWeightBetweenNodes(src, dst);
        std::pair<int, int> my_pair_to_remove = std::make_pair(dst,weight);

        adj[src].remove(my_pair_to_remove);
        my_pair_to_remove.first = src;
        adj[dst].remove(my_pair_to_remove);
    }

    /**
     * Checks if there is an edge between two existing nodes.
     */
    bool hasEdge(int src, int dst)
    {
        return adjacency_matrix[src][dst];
    }

    /**
     * Gets the vector of neighbors associated with the given node.
     */
    std::vector<int> getNeighbors(int node)
    {
        std::vector<int> real_neighbors;
        for (int i = 0; i < size; ++i) {
            if (adjacency_matrix[node][i] != INF && adjacency_matrix[node][i] != 0) {
                real_neighbors.push_back(i);
            }
        }
        return real_neighbors;
    }

    /**
     * Gets the graph size.
     */
    int getSize()
    {
        return size;
    }

    /**
     * Gets the edges number.
    */
    int getEdges()
    {
       return my_edges.getSize();
    }

    /**
     * Get the biggest weight
    */
    int getMaxWight() {
        return max_distance;
    } 

    Edge getEdgeByIndex(int index) {
        return my_edges.getEdgeByIndex(index);
    }

    // BFS algorithm
    // Add the first element in the queue and mark it as visited
    // While the queue is not empty, get the first element, pop it
    // and get the neightbors of the node
    // Add them to the queue if they are not visited
    // Repeat until we fully visited all nodes
    void print_BFS(int root)
    {
        std::cout << "\nBFS traversal: ";

        // Create queue with the root
        std::queue<int> my_queue;
        my_queue.push(root);

        // Create visited vector
        int* visited =  new int[getSize()];
        for (int i = 0; i < size; i++)
            visited[i] = 0;

        visited[root] = 1;
        while (!my_queue.empty()) {
            // Get all neighbors of the current element and display our element
            std::vector<int> my_neighbors = getNeighbors(my_queue.front());
            std::cout << my_queue.front() << " ";
            my_queue.pop();

            // For each neighbour, if it is not visited, visit it
            for (auto& element : my_neighbors) {
                if (visited[element] == 0) {
                    visited[element] = 1;
                    my_queue.push(element);
                }
            }
        }
        delete[] visited;
        std::cout << "\n";
    }

    // DFS algorithm
    // Start from a root
    // Use a stack to start with the root
    // Find the first not-visited node
    // If we can't, pop the stack
    // Repeat until all nodes have been visited
    void print_DFS(int root)
    {
        std::cout << "\nDFS traversal: ";

        // Create visited vector
        int* visited =  new int[getSize()];
        for (int i = 0; i < size; i++)
            visited[i] = 0;

        // Create stack of elements
        std::stack<int> my_stack;

        // Add root and mark as visited
        visited[root]++;
        my_stack.push(root);
        std::cout << my_stack.top() << " ";

        // Repeat until the stack is empty
        while (!my_stack.empty()) {
            // Get neighbors
            std::vector<int> my_neighbors = getNeighbors(my_stack.top());

            // Find the virst non-visited neighbour
            int my_nonvisited_neighbour = -1;
            for (int i = 0; i < my_neighbors.size(); ++i) {
                if (!visited[my_neighbors[i]]) {
                    my_nonvisited_neighbour = my_neighbors[i];
                    break;
                }
            }

            // Check the first neighbour
            if (my_nonvisited_neighbour != -1) {
                visited[my_nonvisited_neighbour]++;
                my_stack.push(my_nonvisited_neighbour);
                std::cout << my_stack.top() << " ";
            }
            else {
                my_stack.pop();
            }
        }
        delete[] visited;
        std::cout << "\n";
    }

    void topo_DFS(int root, int visited[], std::list<std::pair<int, int> >& my_nodes_finish_time)
    {
        // Initialise finish counter
        static int finish_time_counter = 1;

        // Create stack of elements
        std::stack<int> my_stack;

        // Add root and mark as visited
        visited[root]++;
        my_stack.push(root);

        // Repeat until the stack is empty
        while (!my_stack.empty()) {

            // Get neighbors
            std::vector<int> my_neighbors = getNeighbors(my_stack.top());

            // Find the virst non-visited neighbour
            int my_nonvisited_neighbour = -1;
            for (int i = 0; i < my_neighbors.size(); ++i) {
                if (!visited[my_neighbors[i]]) {
                    my_nonvisited_neighbour = my_neighbors[i];
                    break;
                }
            }

            // Check the first neighbour
            if (my_nonvisited_neighbour != -1) {
                visited[my_nonvisited_neighbour]++;
                my_stack.push(my_nonvisited_neighbour);
                ++finish_time_counter;
            }
            else { // if we have no neighbour to visit, pop the stack
                ++finish_time_counter;
                my_nodes_finish_time.push_front({ finish_time_counter, my_stack.top() });
                my_stack.pop();
            }
        }
    }

    // Check if graph is bipartite
    bool bipartite_BFS(int root)
    {
        int* levels =  new int[getSize()];
        for (int i = 0; i < this->size; ++i)
            levels[i] = 0;

        // Create queue with the root
        std::queue<int> my_queue;
        my_queue.push(root);
        levels[root] = ODD;

        // Create visited vector
        int* visited =  new int[getSize()];
        for (int i = 0; i < size; i++)
            visited[i] = 0;

        visited[root] = 1;
        while (!my_queue.empty()) {
            // Get all neighbors of the current element and display our element
            int parent = my_queue.front();
            std::vector<int> my_neighbors = getNeighbors(my_queue.front());
            my_queue.pop();

            // For each neighbour, if it is not visited, visit it
            for (auto& element : my_neighbors) {
                if (visited[element] == 0) {
                    visited[element] = 1;
                    my_queue.push(element);
                    levels[parent] == EVEN ? levels[element] = ODD : levels[element] = EVEN;
                }
                else { // check neighbors parity
                    std::vector<int> my_neighbors_of_element = getNeighbors(element);
                    for (auto& new_element : my_neighbors_of_element) {
                        if (levels[new_element] == levels[element])
                            return false;
                    }
                }
            }
        }
        
        delete[] visited;
        delete[] levels;

        return true;
    }

    void printMatrixGraph()
    {
        // Display top row
        std::cout << "   ";
        for (int i = 0; i < size; ++i) {
            std::cout << i << " ";
        }

        std::cout << "\n";

        // Display matrix
        for (int i = 0; i < size; ++i) {
            std::cout << i << ": ";


            for (int j = 0; j < size; ++j) {
                std::cout << adjacency_matrix.at(i).at(j) << " ";
            }
            std::cout << "\n";
        }
    }

    auto operator[](int line) -> std::vector<int> {
        return adjacency_matrix.at(line);
    }

    auto get_element(int line, int column) -> int {
        return adjacency_matrix.at(line).at(column);
    }
};

#endif // __MATRIXGRAPH_H__