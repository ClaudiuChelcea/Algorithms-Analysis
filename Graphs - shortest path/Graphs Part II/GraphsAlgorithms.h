#pragma once
#include <iostream>
#include <vector>
#include <limits.h>
#include <fstream>
#include "ListGraph.h"
using namespace std;

#define INFINIT 999999

auto miniDist(std::vector<int> distances, std::vector<int> visited, int nodes_nr) -> int // finding minimum distance
{
    int minimum = INT_MAX, ind;
              
    for(int i = 0; i < nodes_nr; ++i) 
    {
        if(!visited.at(i) && distances.at(i) <= minimum)      
        {
            minimum = distances.at(i);
            ind = i;
        }
    }
    return ind;
}

void Dijkstra_shortests_paths(MatrixGraph* my_graph, int node) // adjacency matrix 
{
    std::vector<int> distances;
	distances.resize(my_graph->getSize()); // // array to calculate the minimum distance for each node                             
    
	std::vector<int> visited;
	visited.resize(my_graph->getSize()); // boolean array to mark visited and unvisited for each node
    
    for(int i = 0; i < my_graph->getSize(); ++i)
    {
		distances.at(i) = INT32_MAX;
		visited.at(i) = false;
    }
    
    distances[node] = 0;   // Source vertex distance is set 0               
    
    for(int i = 0; i < my_graph->getSize(); ++i)                           
    {
        int m = miniDist(distances, visited, my_graph->getSize()); 
        visited.at(m) = true;
        for(int j = 0; j < my_graph->getSize(); ++j)                  
        {
            // updating the distance of neighbouring vertex
            if(!visited.at(j) && my_graph->get_element(m,j) && distances.at(m) !=INT_MAX && ((distances.at(m) + my_graph->get_element(m,j)) < distances[j]))
                distances.at(j) = distances.at(m) + my_graph->get_element(m,j);
        }
    }

    for(int i = 0; i < my_graph->getSize(); ++i)                      
    {  
        std:: cout <<distances.at(i)<<std::endl;
    }
	
	distances.clear();
	visited.clear();
}

void runDijkstra(std::string in, std::string out) {
	// Open input file
    std::ifstream fileIN;
    fileIN.open(in, std::ios::in);
    assert(fileIN.is_open() && "File couldn't open!");

	// Create nodes
	int num_nodes = 0;
	fileIN >> num_nodes;
	MatrixGraph* my_graph = new MatrixGraph(num_nodes);

	// Add edges
	int src = 0, dest = 0, distance = 1;
	while (fileIN >> src >> dest >> distance) {
		
		if (src >= num_nodes || dest >= num_nodes) {
			std::cerr << "Error! The source or the destination node doesn't exist!";
			break;
		}

		if (src == dest && dest == -1)
			break;

		my_graph->addEdge(src, dest, distance);
	}

	my_graph->printMatrixGraph();

	// Apply algo Dijkstra
	std::cout << "\nDijkstra:\n";
	int wanted_node{0};
	fileIN >> wanted_node;
	Dijkstra_shortests_paths(my_graph, wanted_node);

	std::cout << "\nRef for Dijkstra:\n";
	std::ifstream fileREF;
    fileREF.open(out, std::ios::in);
    assert(fileREF.is_open() && "File couldn't open!");
	int output;
	while(fileREF >> output) {
		std::cout << output << "\n";
	}

	// Close file
    fileIN.close();
	fileREF.close();
}

void BellmanFord(MatrixGraph* my_graph, int node)
{
    int V = my_graph->getSize();
    int E = my_graph->getEdges();
    int dist[V];
  
    // Step 1: Initialize distances from node to all other vertices
    // as INFINITE
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX;
    dist[node] = 0;
  
    // Step 2: Relax all edges |V| - 1 times. A simple shortest
    // path from node to any other vertex can have at-most |V| - 1
    // edges
    for (int i = 1; i < V; i++) {
        for (int j = 0; j < E; j++) {
            int u = my_graph->getEdgeByIndex(j).src; // sursa laturii 0
            int v = my_graph->getEdgeByIndex(j).dest; // destinatia laturii 0
            int weight = my_graph->getEdgeByIndex(j).distance;
            if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
                dist[v] = dist[u] + weight;
        }
    }
  
    // Step 3: check for negative-weight cycles.  The above step
    // guarantees shortest distances if my_graph doesn't contain
    // negative weight cycle.  If we get a shorter path, then there
    // is a cycle.
    for (int i = 0; i < E; i++) {
        int u = my_graph->getEdgeByIndex(i).src; // sursa laturii 0
            int v = my_graph->getEdgeByIndex(i).dest; // destinatia laturii 0
            int weight = my_graph->getEdgeByIndex(i).distance;
        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
            printf("my_graph contains negative weight cycle");
            return; // If negative cycle is detected, simply return
        }
    }
  
    for (int i = 0; i < V; ++i)
        printf("%d\n", dist[i]);
}

void runBellmanFord(std::string in, std::string out) {
	// Open input file
    std::ifstream fileIN;
    fileIN.open(in, std::ios::in);
    assert(fileIN.is_open() && "File couldn't open!");

	// Create nodes
	int num_nodes = 0;
	fileIN >> num_nodes;
	MatrixGraph* my_graph = new MatrixGraph(num_nodes);

	// Add edges
	int src = 0, dest = 0, distance = 1;
	while (fileIN >> src >> dest >> distance) {
		
		if (src >= num_nodes || dest >= num_nodes) {
			std::cerr << "Error! The source or the destination node doesn't exist!";
			break;
		}

		if (src == dest && dest == -1)
			break;

		my_graph->addEdge(src, dest, distance);
	}

	my_graph->printMatrixGraph();

	// Apply algo Dijkstra
	std::cout << "\nBellmanFord:\n";
	int wanted_node{0};
	fileIN >> wanted_node;
	BellmanFord(my_graph, wanted_node);

	std::cout << "\nRef for BellmanFord:\n";
	std::ifstream fileREF;
    fileREF.open(out, std::ios::in);
    assert(fileREF.is_open() && "File couldn't open!");
	int output;
	while(fileREF >> output) {
		std::cout << output << "\n";
	}

	// Close file
    fileIN.close();
	fileREF.close();
}

void Dijkstra_optimized(MatrixGraph* my_graph, int node) {
	int W = my_graph->getMaxWight();
	int V = my_graph->getSize();
	auto adj = my_graph->getAdj();

	/* With each distance, iterator to that vertex in
       its bucket is stored so that vertex can be deleted
       in O(1) at time of updation. So
    dist[i].first = distance of ith vertex from node vertex
    dits[i].second = iterator to vertex i in bucket number */
    std::vector<std::pair<int, std::list<int>::iterator>> dist(V);
  
    // Initialize all distances as infinite (INF)
    for (int i = 0; i < V; i++)
        dist[i].first = INFINIT;
  
    // Create buckets B[].
    // B[i] keep vertex of distance label i
    std::list<int> B[W * V + 1];
  
    B[0].push_back(node);
    dist[node].first = 0;
  
    //
    int idx = 0;
    while (1)
    {
        // Go sequentially through buckets till one non-empty
        // bucket is found
        while (B[idx].size() == 0 && idx < W*V)
            idx++;
  
        // If all buckets are empty, we are done.
        if (idx == W * V)
            break;
  
        // Take top vertex from bucket and pop it
        int u = B[idx].front();
        B[idx].pop_front();
  
        // Process all adjacents of extracted vertex 'u' and
        // update their distanced if required.
        for (auto i = adj[u].begin(); i != adj[u].end(); ++i)
        {
            int v = (*i).first;
            int weight = (*i).second;
  
            int du = dist[u].first;
            int dv = dist[v].first;
  
            // If there is shorted path to v through u.
            if (dv > du + weight)
            {
                // If dv is not INF then it must be in B[dv]
                // bucket, so erase its entry using iterator
                // in O(1)
                if (dv != INFINIT)
                    B[dv].erase(dist[v].second);
  
                //  updating the distance
                dist[v].first = du + weight;
                dv = dist[v].first;
  
                // pushing vertex v into updated distance's bucket
                B[dv].push_front(v);
  
                // storing updated iterator in dist[v].second
                dist[v].second = B[dv].begin();
            }
        }
    }
  
    // Print shortest distances stored in dist[]
    for (int i = 0; i < V; ++i)
        printf("%d\n", dist[i].first);
}

void runDijkstra_optimized(std::string in, std::string out) {
	// Open input file
    std::ifstream fileIN;
    fileIN.open(in, std::ios::in);
    assert(fileIN.is_open() && "File couldn't open!");

	// Create nodes
	int num_nodes = 0;
	fileIN >> num_nodes;
	MatrixGraph* my_graph = new MatrixGraph(num_nodes);

	// Add edges
	int src = 0, dest = 0, distance = 1;
	while (fileIN >> src >> dest >> distance) {
		
		if (src >= num_nodes || dest >= num_nodes) {
			std::cerr << "Error! The source or the destination node doesn't exist!";
			break;
		}

		if (src == dest && dest == -1)
			break;

		my_graph->addEdge(src, dest, distance);
	}

	my_graph->printMatrixGraph();

	// Apply algo Dijkstra
	std::cout << "\nDijkstra optimized:\n";
	int wanted_node{0};
	fileIN >> wanted_node;
	Dijkstra_optimized(my_graph, wanted_node);

	std::cout << "\nRef for Dijkstra optimized:\n";
	std::ifstream fileREF;
    fileREF.open(out, std::ios::in);
    assert(fileREF.is_open() && "File couldn't open!");
	int output;
	while(fileREF >> output) {
		std::cout << output << "\n";
	}

	// Close file
    fileIN.close();
	fileREF.close();
}

void topologicalSortUtil(int v, bool visited[], std::stack<int> &Stack, ListGraph* my_list_graph, MatrixGraph* my_graph)
{
    // Mark the current node as visited
    visited[v] = true;
 
    // Recur for all the vertices adjacent to this vertex
    auto neighbors = my_graph->getNeighbors(v);
    for(int i = 0; i < neighbors.size(); ++i) {
        if(!visited[neighbors.at(i)]) {
            topologicalSortUtil(neighbors.at(i), visited, Stack, my_list_graph, my_graph);
        }
    }
 
    // Push current vertex to stack which stores topological sort
    Stack.push(v);
}
 
// The function to find shortest paths from given vertex. It uses recursive
// topologicalSortUtil() to get topological sorting of given graph.
void ShortestPathAcyclicOrientedGraph(MatrixGraph* my_graph, int s)
{
    std::stack<int> Stack;
	int V = my_graph->getSize();
	ListGraph* my_list_graph = new ListGraph(my_graph);
    int dist[V];
    for(int i = 0; i < V; ++i)
        dist[i] = INFINIT;
 
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
 
    // Call the recursive helper function to store Topological Sort
    // starting from all vertices one by one
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            topologicalSortUtil(i, visited, Stack, my_list_graph, my_graph);
 
    // Initialize distances to all vertices as infinite and distance
    // to source as 0
    for (int i = 0; i < V; i++)
        dist[i] = INFINIT;
    dist[s] = 0;
 
    // Process vertices in topological order
    while (Stack.empty() == false)
    {
        // Get the next vertex from topological order
        int u = Stack.top();
        Stack.pop();
 
        // Update distances of all adjacent vertices
        if (dist[u] != INFINIT) {
            auto adj = my_graph->getNeighbors(u);
            for(int i = 0; i < adj.size(); ++i) {
                int weight = my_graph->getWeightBetweenNodes(u, adj[i]);
                if(dist[adj.at(i)] > dist[u] + weight)
                    dist[adj.at(i)] = dist[u] + weight;
            }
        }
    }
 
    // Print the calculated shortest distances
    for (int i = 0; i < V; i++)
        std::cout << dist[i] << " ";
}

void runShortestPathAcyclicOrientedGraph(std::string in, std::string out) {
	// Open input file
    std::ifstream fileIN;
    fileIN.open(in, std::ios::in);
    assert(fileIN.is_open() && "File couldn't open!");

	// Create nodes
	int num_nodes = 0;
	fileIN >> num_nodes;
	MatrixGraph* my_graph = new MatrixGraph(num_nodes);

	// Add edges
	int src = 0, dest = 0, distance = 1;
	while (fileIN >> src >> dest >> distance) {
		
		if (src >= num_nodes || dest >= num_nodes) {
			std::cerr << "Error! The source or the destination node doesn't exist!";
			break;
		}

		if (src == dest && dest == -1)
			break;

		my_graph->addEdge(src, dest, distance);
	}

	my_graph->printMatrixGraph();

	// Apply algo Dijkstra
	std::cout << "\nShortest Path Acyclic Oriented Graph:\n";
	int wanted_node{0};
	fileIN >> wanted_node;
	ShortestPathAcyclicOrientedGraph(my_graph, wanted_node);

	std::cout << "\nRef for Shortest Path Acyclic Oriented Graph:\n";
	std::ifstream fileREF;
    fileREF.open(out, std::ios::in);
    assert(fileREF.is_open() && "File couldn't open!");
	int output;
	while(fileREF >> output) {
		std::cout << output << "\n";
	}

	// Close file
    fileIN.close();
	fileREF.close();
}
