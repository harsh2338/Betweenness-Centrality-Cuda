// %%cu
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
#include <iostream>
#include <cuda.h>
#define MAX_THREAD_COUNT 1024
#define MAX_MEMORY ((long long)4e9)

void make_graph(int LIM){
        int x, y;
    int n=LIM;
    int m=LIM;
    int cur = 0;

    freopen("input.txt", "w", stdout);
    cout << n << " " << m << endl;

    vector<int> *v = new vector<int>[n+1];
    for (int i = 0; i < m; i++) {
        do {
            x = rand() % n;
            y = rand() % n;
        } while (x == y);

        v[x].push_back(y);
        v[y].push_back(x);
    }

    for(int i=0; i<n; i++)
        for(int node : v[i])
            cout<<node<<" ";
    cout<<"\n";

    cout<<"0 ";
    for(int i=0; i<n; i++) {
        cur += v[i].size();
        cout<<cur<<" ";
    }
    cout<<"\n";
}

using namespace std;

__global__ void betweennessCentralityKernel(Graph *graph, float *bc, int nodeCount,
            int *sigma, int *distance, float *dependency, int *Q, int *Qpointers) {
    
    int idx = threadIdx.x;
    if(idx >= nodeCount)
        return;
    
    __shared__ int s;
    __shared__ int Q_len;
    __shared__ int Qpointers_len;
    __shared__ int noOfBlocks;

    if(idx == 0) {
        s = blockIdx.x - gridDim.x;
        noOfBlocks = gridDim.x;
    }
    __syncthreads();
    
    while(s < nodeCount - noOfBlocks)
    {
        if(idx == 0)
        {
            s += noOfBlocks;

            Q[0 + (blockIdx.x * nodeCount)] = s;
            Q_len = 1;
            Qpointers[0 + (blockIdx.x * nodeCount)] = 0;
            Qpointers[1 + (blockIdx.x * nodeCount)] = 1;
            Qpointers_len = 1;
        }
        __syncthreads();

        for(int v=idx; v<nodeCount; v+=blockDim.x)
        {
            if(v == s)
            {
                distance[v + (blockIdx.x * nodeCount)] = 0;
                sigma[v + (blockIdx.x * nodeCount)] = 1;
            }
            else
            {
                distance[v + (blockIdx.x * nodeCount)] = INT_MAX;
                sigma[v + (blockIdx.x * nodeCount)] = 0;
            }
            dependency[v + (blockIdx.x * nodeCount)] = 0.0;
        }
        __syncthreads();
        
        while(true)
        {
            __syncthreads();
            for(int k=idx; k<Qpointers[Qpointers_len + (blockIdx.x * nodeCount)]; k+=blockDim.x) 
            {
                if(k < Qpointers[Qpointers_len -1 + (blockIdx.x * nodeCount)])
                    continue;

                int v = Q[k + (blockIdx.x * nodeCount)];
                for(int r = graph->adjListPointers[v]; r < graph->adjListPointers[v + 1]; r++)
                {
                    int w = graph->adjList[r];
                    if(atomicCAS(&distance[w + (blockIdx.x * nodeCount)], INT_MAX, distance[v + (blockIdx.x * nodeCount)] +1) == INT_MAX)
                    {
                        int t = atomicAdd(&Q_len, 1);
                        Q[t + (blockIdx.x * nodeCount)] = w;
                    }
                    if(distance[w + (blockIdx.x * nodeCount)] == (distance[v + (blockIdx.x * nodeCount)]+1))
                    {
                        atomicAdd(&sigma[w + (blockIdx.x * nodeCount)], sigma[v + (blockIdx.x * nodeCount)]);
                    }
                }
            }
            __syncthreads();

            if(Q_len == Qpointers[Qpointers_len + (blockIdx.x * nodeCount)])
                break;

            if(idx == 0)
            {
                Qpointers_len++;
                Qpointers[Qpointers_len + (blockIdx.x * nodeCount)] = Q_len;
            }
            __syncthreads();
        }
        __syncthreads();
        
        while(Qpointers_len > 0)
        {
            for(int k=idx; k < Qpointers[Qpointers_len + (blockIdx.x * nodeCount)]; k+=blockDim.x) 
            {
                if(k < Qpointers[Qpointers_len -1 + (blockIdx.x * nodeCount)])
                    continue;

                int v = Q[k + (blockIdx.x * nodeCount)];
                for(int r = graph->adjListPointers[v]; r < graph->adjListPointers[v + 1]; r++)
                {
                    int w = graph->adjList[r];
                    if(distance[w + (blockIdx.x * nodeCount)] == (distance[v + (blockIdx.x * nodeCount)] + 1))
                    {
                        if (sigma[w + (blockIdx.x * nodeCount)] != 0)
                            dependency[v + (blockIdx.x * nodeCount)] += (sigma[v + (blockIdx.x * nodeCount)] * 1.0 / sigma[w + (blockIdx.x * nodeCount)]) * (1 + dependency[w + (blockIdx.x * nodeCount)]);
                    }
                }
                if (v != s)
                {
                    atomicAdd(bc + v, dependency[v + (blockIdx.x * nodeCount)] / 2);
                }
            }
            __syncthreads();

            if(idx == 0)
                Qpointers_len--;

            __syncthreads();
        }
    }
}

float *betweennessCentrality(Graph *graph, int nodeCount)
{
    float *bc = new float[nodeCount]();
    float *device_bc, *dependency;
    int *sigma, *distance, *Q, *Qpointers;

    const int BLOCK_COUNT = MAX_MEMORY / (4 * 5 * nodeCount);

    cudaMalloc((void **)&device_bc, sizeof(float) * nodeCount);
    cudaMalloc((void **)&sigma, sizeof(int) * nodeCount * BLOCK_COUNT);
    cudaMalloc((void **)&distance, sizeof(int) * nodeCount * BLOCK_COUNT);
    cudaMalloc((void **)&Q, sizeof(int) * (nodeCount) * BLOCK_COUNT);
    cudaMalloc((void **)&Qpointers, sizeof(int) * (nodeCount) * BLOCK_COUNT);
    cudaMalloc((void **)&dependency, sizeof(float) * nodeCount * BLOCK_COUNT);
    cudaMemcpy(device_bc, bc, sizeof(float) * nodeCount, cudaMemcpyHostToDevice);

    betweennessCentralityKernel<<<BLOCK_COUNT, MAX_THREAD_COUNT>>>(graph, device_bc, nodeCount, sigma, distance, dependency, Q, Qpointers);
    cudaDeviceSynchronize();

    cudaMemcpy(bc, device_bc, sizeof(float) * nodeCount, cudaMemcpyDeviceToHost);
    cudaFree(device_bc);
    cudaFree(sigma);
    cudaFree(dependency);
    cudaFree(distance);
    cudaFree(Q);
    cudaFree(Qpointers);
    return bc;
}
double getWD(){
        freopen("input.txt", "r", stdin);

    Graph *graph = new Graph();
    Graph *cudaGraph;

    cudaMalloc((void **)&cudaGraph, sizeof(Graph));
    graph->readGraph();

    int numNodes = graph->getNumNodes();
    int numEdges = graph->getNumEdges();
    cudaMemcpy(cudaGraph, graph, sizeof(Graph), cudaMemcpyHostToDevice);

    int *adjList;
    cudaMalloc((void **)&adjList, sizeof(int) * (2 * numEdges + 1));
    cudaMemcpy(adjList, graph->adjList, sizeof(int) * (2 * numEdges + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(&(cudaGraph->adjList), &adjList, sizeof(int *), cudaMemcpyHostToDevice);

    int *adjListPointers;
    cudaMalloc((void **)&adjListPointers, sizeof(int) * (numNodes + 1));
    cudaMemcpy(adjListPointers, graph->adjListPointers, sizeof(int) * (numNodes + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(&(cudaGraph->adjListPointers), &adjListPointers, sizeof(int *), cudaMemcpyHostToDevice);
 
    clock_t start, end;
    start = clock();
    float *bc = betweennessCentrality(cudaGraph, numNodes);
    end = clock();
    double time_taken = 1000.0 * (end - start) / (double)CLOCKS_PER_SEC;
    double maxBetweenness = -1, maxBetweennessNode=-1;
    cout<<"Time taken :"<<time_taken<<"\n";
    // cout<<"Betweenness Centrality for each node:\n";
    for (int i = 0; i < numNodes; i++)
    {
        if(maxBetweenness<bc[i]){
            maxBetweenness = bc[i];
            maxBetweennessNode=i;
        }
        // cout<<i<<" : "<<bc[i]<<"\n";
    }
    cout<<"Maximum betweenness centrality : ";
    cout<<maxBetweenness<<"\n";
    cout<<"Node with Maximum betweenness centrality : ";
    cout<<maxBetweennessNode<<"\n";

    cudaFree(adjList);
    cudaFree(adjListPointers);
    cudaFree(cudaGraph);
    return time_taken;
}
int main(){
    vector<double>x,v,e,we,wd;
    for(int i=500;i<=1500;i+=500){
        make_graph(i);
        x.push_back(i);
        //v.push_back(getVertex());
        //e.push_back(getEdge());
        //we.push_back(getWE());
        wd.push_back(getWD());
    }
    for(int i:x)cout<<i<<" , ";
    cout<<"\n\n";
    for(int i:v)cout<<i<<" , ";
        cout<<"\n\n";


    for(int i:e)cout<<i<<" , ";

    cout<<"\n\n";

    for(int i:we)cout<<i<<" , ";
    cout<<"\n\n";

    for(int i:wd)cout<<i<<" , ";
    cout<<"\n\n";
    
    
}