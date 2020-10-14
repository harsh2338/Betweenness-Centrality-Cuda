#include "graph.h"
using namespace std;

double *betweennessCentrality(Graph *graph)
{
    int numNodes = graph->getNumNodes();

    double *bc = new double[numNodes]();
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

    double *bc = betweennessCentrality(graph);

    end = clock();
    float time_taken = 1000.0 * (end - start) / (float)CLOCKS_PER_SEC;
    double maxBetweenness = -1;
    cout<<"Time taken :"<<time_taken<<"\n";
    cout<<"Betweenness Centrality for node:\n";
    for (int i = 0; i < numNodes; i++)
    {
        maxBetweenness = max(maxBetweenness, bc[i]);
        cout<<i<<" : "<<bc[i]<<"\n";
    }
    cout<<maxBetweenness;
    return 0;
}
