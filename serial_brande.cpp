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
float *betweennessCentrality(Graph *graph)
{
    int numNodes = graph->getNumNodes();

    float *bc = new float[numNodes]();
    vector<int> *predecessor = new vector<int>[numNodes];

    double *dependency = new double[numNodes];
    int *sigma = new int[numNodes];
    int *distance = new int[numNodes];

    for (int s = 0; s < numNodes; s++)
    {
        stack<int> st;
        memset(distance, -1, numNodes * sizeof(int));
        memset(sigma, 0, numNodes * sizeof(int));
        memset(dependency, 0, numNodes * sizeof(double));

        distance[s] = 0;
        sigma[s] = 1;
        queue<int> q;
        q.push(s);
        while (!q.empty())
        {
            int v = q.front();
            q.pop();
            st.push(v);
            for (int i = graph->adjListPointers[v]; i < graph->adjListPointers[v + 1]; i++)
            {
                int w = graph->adjList[i];

                if (distance[w] ==-1)
                {
                    q.push(w);
                    distance[w] = distance[v] + 1;
                }

                if (distance[w] == distance[v] + 1)
                {
                    sigma[w] += sigma[v];
                    predecessor[w].push_back(v);
                }
            }
        }
        //4 1 3 2 0
        while (!st.empty())
        {
            int w = st.top();
            st.pop();

            for (const int &v : predecessor[w])
            {
                if (sigma[w] != 0)
                    dependency[v] += (sigma[v] * 1.0 / sigma[w]) * (1 + dependency[w]);
            }
            if (w != s)
                bc[w] += dependency[w] / 2;
            
        }
        for(int i=0; i<numNodes; ++i)
            predecessor[i].clear();
        
    }

    return bc;
}

int main()
{

    freopen("input.txt", "r", stdin);
    Graph *graph = new Graph();
    graph->readGraph();

    int numNodes = graph->getNumNodes();
    int numEdges = graph->getNumEdges();

    clock_t start, end;
    start = clock();

    float *bc = betweennessCentrality(graph);

    end = clock();
    float time_taken = 1000.0 * (end - start) / (float)CLOCKS_PER_SEC;
    double maxBetweenness = -1, maxBetweennessNode=-1;
    cout<<"Time taken :"<<time_taken<<"\n";
    //cout<<"Betweenness Centrality for each node:\n";
    for (int i = 0; i < numNodes; i++)
    {
        if(maxBetweenness<bc[i]){
            maxBetweenness = bc[i];
            maxBetweennessNode=i;
        }
        cout<<i<<" : "<<bc[i]<<"\n";
    }
    cout<<"Maximum betweenness centrality : ";
    cout<<maxBetweenness<<"\n";
    cout<<"Node with Maximum betweenness centrality : ";
    cout<<maxBetweennessNode<<"\n";
    return 0;
}