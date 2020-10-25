#include <bits/stdc++.h>
using namespace std;
class Graph {

public:

	int numNodes, numEdges;
	int *adjList, *adjListPointers;
	int *edgeList1, *edgeList2;

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

	void convertToCOO() {
		edgeList2 = adjList;
		edgeList1 = new int[2 * numEdges +1];

		for(int i=0; i <numNodes; ++i) {
			for(int j=adjListPointers[i]; j<adjListPointers[i+1]; ++j){
				edgeList1[j] = i;
			}
		}
	}
	int *getadjList(int node) {
		return adjList;
	}

	int *getadjListPointers(int node) {
		return adjListPointers;
	}

};