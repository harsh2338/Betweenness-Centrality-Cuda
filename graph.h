#include <bits/stdc++.h>
using namespace std;
class Graph {

public:

	int numNodes, numEdges;
	int *adjList, *adjListPointers;

	int getNumNodes() {
		return numNodes;
	}

	int getNumEdges() {
		return numEdges;
	}

	void readGraph() {
		cin >> numNodes >> numEdges;
		adjListPointers = new int[numNodes +1];
		adjList = new int[2 * numEdges +1];
		for(int i=0; i<(2 * numEdges); i++)
			cin >> adjList[i];
		for(int i=0; i<=numNodes; i++) 
			cin >> adjListPointers[i];

	}

	int *getadjList(int node) {
		return adjList;
	}

	int *getadjListPointers(int node) {
		return adjListPointers;
	}

};