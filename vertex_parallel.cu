#include "graph.h"
#include <iostream>
#include <cuda.h>
#define MAX_THREAD_COUNT 1024

using namespace std;


__global__ void betweennessCentralityKernel(Graph *graph, double *bc, int numNodes,
            int *sigma, int *distance, double *dependency) {
    
    int idx = threadIdx.x;
    if(idx >= numNodes)
        return;
    
    __shared__ int s;
    __shared__ int current_depth;
    __shared__ bool done;

    if(idx == 0) {
        s = -1;
    }
    __syncthreads();

    while(s < numNodes -1)
    {    
        if(idx == 0)
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
        
        while(!done)
        {
            if(idx == 0){
                current_depth++;
            }
            done = true;
            __syncthreads();

            for(int v=idx; v<numNodes; v+=blockDim.x)
            {
                if(distance[v] == current_depth)
                {
                    for(int r = graph->adjListPointers[v]; r < graph->adjListPointers[v + 1]; r++)
                    {
                        int w = graph->adjList[r];
                        if(distance[w] == INT_MAX)
                        {
                            distance[w] = distance[v] + 1;
                            done = false;
                        }
                        if(distance[w] == (distance[v] + 1))
                        {
                            atomicAdd(&sigma[w], sigma[v]);
                        }
                    }
                }
            }
            __syncthreads();
        }

        while(current_depth)
        {
            if(idx == 0){
                current_depth--;
            }
            __syncthreads();

            for(int v=idx; v<numNodes; v+=blockDim.x) 
            {
                if(distance[v] == current_depth)
                {
                    for(int r = graph->adjListPointers[v]; r < graph->adjListPointers[v + 1]; r++)
                    {
                        int w = graph->adjList[r];
                        if(distance[w] == (distance[v] + 1))
                        {
                            if (sigma[w] != 0)
                                dependency[v] += (sigma[v] * 1.0 / sigma[w]) * (1 + dependency[w]);
                        }
                    }
                    if (v != s)
                    {
                        bc[v] += dependency[v] / 2;
                    }
                }
            }
            __syncthreads();
        }
    }
}

double *betweennessCentrality(Graph *graph, int numNodes)
{
    double *bc = new double[numNodes]();
    double *cuda_bc, *dependency;
    int *sigma, *distance;

    cudaMalloc((void **)&cuda_bc, sizeof(double) * numNodes);
    cudaMalloc((void **)&sigma, sizeof(int) * numNodes);
    cudaMalloc((void **)&distance, sizeof(int) * numNodes);
    cudaMalloc((void **)&dependency, sizeof(double) * numNodes);
    cudaMemcpy(cuda_bc, bc, sizeof(double) * numNodes, cudaMemcpyHostToDevice);


    betweennessCentralityKernel<<<1, MAX_THREAD_COUNT>>>(graph, cuda_bc, numNodes, sigma, distance, dependency);
    cudaDeviceSynchronize();

    cudaMemcpy(bc, cuda_bc, sizeof(double) * numNodes, cudaMemcpyDeviceToHost);
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
    double *bc = betweennessCentrality(cudaGraph, numNodes);
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
}