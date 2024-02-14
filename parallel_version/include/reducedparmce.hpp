#ifndef _PARMCE_H_
#define _PARMCE_H_

#include <iostream>
#include <string>
#include <tbb/atomic.h>
#include "graph.hpp"

using namespace std;

inline atomic<long> cliqueCount(0);
inline vector<vector<bool>> allBfsParent;
inline vector<vector<int>> allDegCnt;
// inline tbb::concurrent_vector<tbb::atomic<int>> debugInfo;

void reduced_parmce_main(string file_path, int num_thread);

void reduced_parTTT(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X);

int reduced_parPivot(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X);

#endif 