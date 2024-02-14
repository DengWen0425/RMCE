#ifndef _PARMCE_H_
#define _PARMCE_H_

#include <iostream>
#include <string>
#include <tbb/atomic.h>
#include "graph.hpp"

using namespace std;

inline atomic<long> cliqueCount(0);

void parmce_main(string file_path, int num_thread);

void parTTT(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X);

int parPivot(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X);

#endif 