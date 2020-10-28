#include "graph.h"
#include <iostream>
#include <cuda.h>
#define MAX_THREAD_COUNT 1024

using namespace std;

__global__ void betweennessCentralityKernel(Graph *graph, double *bwCentrality, int nodeCount,
            int *sigma, int *distance, double *dependency, int *Q, int *Qpointers) {
    
    int idx = threadIdx.x;
    if(idx >= nodeCount)
        return;
    
    __shared__ int s;
    __shared__ int Q_len;
    __shared__ int Qpointers_len;

    if(idx == 0) {
        s = -1;
    }
    __syncthreads();

    while(s < nodeCount -1)
    {    
        if(idx == 0)
        {
            ++s;
            
            Q[0] = s;
            Q_len = 1;
            Qpointers[0] = 0;
            Qpointers[1] = 1;
            Qpointers_len = 1;
        }
        __syncthreads();

        for(int v=idx; v<nodeCount; v+=blockDim.x)
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
        
        while(true)
        {
            __syncthreads();
            for(int k=idx; k<Qpointers[Qpointers_len]; k+=blockDim.x)
            {
                if(k < Qpointers[Qpointers_len -1])
                    continue;

                int v = Q[k];
                for(int r = graph->adjListPointers[v]; r < graph->adjListPointers[v + 1]; r++)
                {
                    int w = graph->adjList[r];
                    if(atomicCAS(&distance[w], INT_MAX, distance[v] +1) == INT_MAX)
                    {
                        int t = atomicAdd(&Q_len, 1);
                        Q[t] = w;
                    }
                    if(distance[w] == (distance[v]+1))
                    {
                        atomicAdd(&sigma[w], sigma[v]);
                    }
                }
            }
            __syncthreads();

            if(Q_len == Qpointers[Qpointers_len])
                break;

            if(idx == 0)
            {
                Qpointers_len++;
                Qpointers[Qpointers_len] = Q_len;
            }
            __syncthreads();
        }
        __syncthreads();
        
        while(Qpointers_len > 0)
        {
            for(int k=idx; k < Qpointers[Qpointers_len]; k+=blockDim.x) 
            {
                if(k < Qpointers[Qpointers_len -1])
                    continue;

                int v = Q[k];
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
                    bwCentrality[v] += dependency[v] / 2;
                }
            }
            __syncthreads();

            if(idx == 0)
                Qpointers_len--;

            __syncthreads();
        }
    }
}

double *betweennessCentrality(Graph *graph, int nodeCount)
{
    double *bwCentrality = new double[nodeCount]();
    double *device_bwCentrality, *dependency;
    int *sigma, *distance, *Q, *Qpointers;

    cudaMalloc((void **)&device_bwCentrality, sizeof(double) * nodeCount);
    cudaMalloc((void **)&sigma, sizeof(int) * nodeCount);
    cudaMalloc((void **)&distance, sizeof(int) * nodeCount);
    cudaMalloc((void **)&Q, sizeof(int) * (nodeCount +1));
    cudaMalloc((void **)&Qpointers, sizeof(int) * (nodeCount +1));
    cudaMalloc((void **)&dependency, sizeof(double) * nodeCount);
    cudaMemcpy(device_bwCentrality, bwCentrality, sizeof(double) * nodeCount, cudaMemcpyHostToDevice);

    betweennessCentralityKernel<<<1, MAX_THREAD_COUNT>>>(graph, device_bwCentrality, nodeCount, sigma, distance, dependency, Q, Qpointers);
    cudaDeviceSynchronize();

    cudaMemcpy(bwCentrality, device_bwCentrality, sizeof(double) * nodeCount, cudaMemcpyDeviceToHost);
    cudaFree(device_bwCentrality);
    cudaFree(sigma);
    cudaFree(dependency);
    cudaFree(distance);
    cudaFree(Q);
    cudaFree(Qpointers);
    return bwCentrality;
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