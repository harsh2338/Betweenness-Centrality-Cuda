#include "graph.h"
#include <iostream>
#include <cuda.h>
#define MAX_THREAD_COUNT 1024

using namespace std;


__global__ void betweennessCentralityKernel(Graph *graph, float *bc, int numNodes,
                                            int *sigma, int *distance, float *dependency)
{

    int idx = threadIdx.x;
    if(idx >= max((2*(graph->numEdges)), numNodes))
        return;

    __shared__ int s;
    __shared__ int current_depth;
    __shared__ bool done;

    if(idx == 0) {
        s = -1;
    }
    __syncthreads();

    while(s < numNodes - 1)
    {
        if (idx == 0)
        {
            ++s;
            done = false;
            current_depth = -1;
        }
        __syncthreads();
    
        for(int v=idx; v<numNodes; v+=blockDim.x)
        {
            if(v == s)
            {
                distance[v] = 0;
                sigma[v] = 1;
            }
            else
            {
                distance[v] = INT_MAX;
                sigma[v] = 0;
            }
            dependency[v] = 0.0;
        }
        __syncthreads();


        while (!done)
        {
            __syncthreads();

            if (threadIdx.x == 0){
                current_depth++;
            }
            done = true;
            __syncthreads();
    
            
            for (int i = idx; i < (2*(graph->numEdges)); i += blockDim.x) 
            {
                int v = graph->edgeList1[i];
                if (distance[v] == current_depth)
                {    
                    int w = graph->edgeList2[i];
                    if (distance[w] == INT_MAX)
                    {
                        distance[w] = distance[v] + 1;
                        done = false;
                    }
                    if (distance[w] == (distance[v] + 1))
                    {
                        atomicAdd(&sigma[w], sigma[v]);
                    }
                }
            }
            __syncthreads();
        }
        
        __syncthreads();
        
        while(current_depth)
        {
            if(idx == 0){
                current_depth--;
            }
            __syncthreads();

            for (int i = idx; i < (2*(graph->numEdges)); i += blockDim.x) 
            {
                int v = graph->edgeList1[i];
                if(distance[v] == current_depth)
                {

                    int w = graph->edgeList2[i];
                    if(distance[w] == (distance[v] + 1))
                    {
                        if (sigma[w] != 0) {
                            atomicAdd(dependency + v, (sigma[v] * 1.0 / sigma[w]) * (1 + dependency[w]));
                        }
                       
                    }
  
                }
            }
            __syncthreads();
        }

        for(int v=idx; v<numNodes; v+=blockDim.x){
            if (v != s)
            {
                bc[v] += dependency[v] / 2;
            }
        }
        __syncthreads();
    }
}

float *betweennessCentrality(Graph *graph, int numNodes)
{
    float *bc = new float[numNodes]();
    float *cuda_bc, *dependency;
    int *sigma, *distance;

    cudaMalloc((void **)&cuda_bc, sizeof(float) * numNodes);
    cudaMalloc((void **)&sigma, sizeof(int) * numNodes);
    cudaMalloc((void **)&distance, sizeof(int) * numNodes);
    cudaMalloc((void **)&dependency, sizeof(float) * numNodes);
    cudaMemcpy(cuda_bc, bc, sizeof(float) * numNodes, cudaMemcpyHostToDevice);

    betweennessCentralityKernel<<<1, MAX_THREAD_COUNT>>>(graph, cuda_bc, numNodes, sigma, distance, dependency);
    cudaDeviceSynchronize();
    
    cudaMemcpy(bc, cuda_bc, sizeof(float) * numNodes, cudaMemcpyDeviceToHost);
    cudaFree(cuda_bc);
    cudaFree(sigma);
    cudaFree(dependency);
    cudaFree(distance);

    return bc;
}

int main()
{

    freopen("input.txt", "r", stdin);

    Graph *graph = new Graph();
    Graph *cudaGraph;

    cudaMalloc((void **)&cudaGraph, sizeof(Graph));
    graph->readGraph();
    graph->convertToCOO();

    int numNodes = graph->getNumNodes();
    int numEdges = graph->getNumEdges();
    cudaMemcpy(cudaGraph, graph, sizeof(Graph), cudaMemcpyHostToDevice);

    int *edgeList1;
    int *edgeList2;
    cudaMalloc((void **)&edgeList1, sizeof(int) * (2 * numEdges + 1));
    cudaMemcpy(edgeList1, graph->edgeList1, sizeof(int) * (2 * numEdges + 1), cudaMemcpyHostToDevice);

    cudaMalloc((void **)&edgeList2, sizeof(int) * (2 * numEdges + 1));
    cudaMemcpy(edgeList2, graph->edgeList2, sizeof(int) * (2 * numEdges + 1), cudaMemcpyHostToDevice);

    cudaMemcpy(&(cudaGraph->edgeList1), &edgeList1, sizeof(int *), cudaMemcpyHostToDevice);
    cudaMemcpy(&(cudaGraph->edgeList2), &edgeList2, sizeof(int *), cudaMemcpyHostToDevice);

 
    clock_t start, end;
    start = clock();
    float *bc = betweennessCentrality(cudaGraph, numNodes);
    end = clock();
    float time_taken = 1000.0 * (end - start) / (float)CLOCKS_PER_SEC;
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

    cudaFree(cudaGraph);
}