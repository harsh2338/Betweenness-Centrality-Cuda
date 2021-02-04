#include<bits/stdc++.h>
#include "graph.h"
using namespace std;
int main(){
	freopen("input.txt", "r", stdin);
    Graph *graph = new Graph();
    graph->readGraph();
    graph->convertToCOO();
}