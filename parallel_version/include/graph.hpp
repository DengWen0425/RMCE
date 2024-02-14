#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <sys/time.h>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_hash_map.h>
#include <atomic>
#include <thread>

using namespace std;

class Graph
{
public:
    Graph();
    ~Graph();

    void from_path(string const &filename);

    void computeDegeneracyOrder();
    long computeDegeneracyWithGlobalReduction();
    long computeDegeneracyWithGlobalReductionParallel();

    int n;
    int m;
    int degeneracy=-1;
    int* vals;
    int* inds;
    int* seps;
    int* seps2;
    int* orders;
    int* o2v;
    int* orderedVals;
    // bool* bfsParent;
    tbb::concurrent_vector<int> Xpruned;
    // int* Xpruned;
    int offset = 0;
};

class SubGraph
{
public:
    SubGraph();
    ~SubGraph();

    void from_graph(Graph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X);
    void from_graph_reduced(Graph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X);

    tbb::concurrent_unordered_map<int, tbb::concurrent_unordered_set<int>> adj;
    // tbb::concurrent_unordered_map<int, tbb::concurrent_unordered_set<int>> adjInX;
    int anchor;
    // unordered_map<int, unordered_set<int>> adj;
    // tbb::concurrent_unordered_map<int, bool> bfsParent;
    // tbb::concurrent_unordered_map<int, int> neiCnt;
};

struct NeighborList
{
    NeighborList()
    : vertex(-1)
    , earlier()
    , later()
    , orderNumber(-1) {}

    int vertex; //!< the vertex that owns this neighbor list
    std::list<int> earlier; //!< a linked list of neighbors that come before this vertex in the ordering
    std::list<int> later; //!< a linked list of neighbors that come after this vertex in the ordering
    int orderNumber; //!< the position of this vertex in the ordering
};

typedef struct NeighborList NeighborList;

#endif 