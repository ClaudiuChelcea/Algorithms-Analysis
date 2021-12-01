#include "MatrixGraph.h"
#include "GraphsAlgorithms.h"
#include <iostream>
#include <cassert>
#define INF 0

int main(void)
{
	// Choose algorithm
	std::string algorithm;
	std::cout << "Choose algorithm from: Dijkstra | BellmanFord | Dijkstra_optimized | Shortest_path!\n";
	do {
		std::cout << "Insert the name as a perfect match to one of the options above: ";
		std::cin >> algorithm;
	} while(algorithm != "Dijkstra" && algorithm != "BellmanFord" && algorithm != "Dijkstra_optimized" && algorithm != "Shortest_path");

	// Choose test number
	char test_number{0};
	do {
		std::cout << "Insert test number (from 1 to 5): ";
		std::cin>> test_number;
		if(test_number - 48 >= 1 && test_number - 48 <= 5)
			break;
	} while(true);

	// Create path
	std::string path_in = "Tests/" + algorithm + "/test" + test_number;
	std::string path_out = "Refs/" + algorithm + "/test" + test_number+ "_ref";

	// Choose test and algorithm to run
	if(algorithm == "Dijkstra") {
		// Run Dijkstra
		runDijkstra(path_in, path_out);
	} else if(algorithm == "BellmanFord") {
		// run BellmanFord
		runBellmanFord(path_in, path_out);
	} else if(algorithm == "Dijkstra_optimized") {
		// run Dijkstra_optimized
		runDijkstra_optimized(path_in, path_out);
	} else if(algorithm == "Shortest_path") {
		// run Dijkstra_optimized
		runShortestPathAcyclicOrientedGraph(path_in, path_out);
	} else {
		std::cerr << "Couldn't recognize algorithm! TRY AGAIN!\n";
	}

    return 0;
}
