#include <iostream>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_hash_map.h>
#include <sys/time.h>
#include <atomic>
#include <thread>
#include <chrono>

#include "parmce.hpp"


void parmce_main(string file_name, int num_thread)
{
    /* Initialization */
    tbb::task_scheduler_init init(num_thread);
    cout << "Allocated thread number: " << num_thread << endl;

    /* Initilize CSR */
    Graph g;
    g.from_path(file_name);

    std::chrono::steady_clock::time_point start, mid, end;
    start = std::chrono::steady_clock::now();

    /* tail: smaller, Gamma: larger */
    tbb::concurrent_unordered_set<int> P;
    tbb::concurrent_unordered_set<int> X;
    tbb::concurrent_unordered_set<int> R;

    g.computeDegeneracyOrder();

    mid = std::chrono::steady_clock::now();

    /* Use tbb::global control */
    tbb::task_group tg0;
    tg0.run([&] {
        tbb::parallel_for_each(g.o2v+g.offset, g.o2v+g.n, 
            [&](int v) {
                tbb::concurrent_unordered_set<int> P_local;
                tbb::concurrent_unordered_set<int> X_local;

                tbb::task_group tg;
                tg.run([&] {
                tbb::parallel_for(g.inds[v], g.seps[v], 1,
                    [&](int i) {
                        X_local.insert(g.orderedVals[i]);
                    });
                });
                tg.wait();

                tg.run([&] {
                tbb::parallel_for(g.seps[v], g.inds[v+1], 1,
                    [&](int i) {
                        P_local.insert(g.orderedVals[i]);
                    });
                });
                tg.wait();

                SubGraph subg;
                subg.from_graph(g, P_local, X_local);
                parTTT(subg, P_local, X_local);
                
            });
    });
    tg0.wait();

    end = std::chrono::steady_clock::now();
    std::chrono::duration<double> time = (end-start);
    std::chrono::duration<double> time1 = (mid-start);
    cout << "Order compute time: " << time1.count() << " sec" << endl;
    cout << "Execution time: " << time.count() << " sec" << endl;

    cout << "MBE number: " << cliqueCount << endl;
}

void parTTT(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X)
{
    if (P.empty() && X.empty()) {
        cliqueCount++;
        return;
    }

    int pivot = parPivot(g, P, X);
    // cout << "pivot: " << pivot << endl;
    tbb::concurrent_vector<int> candidates_vec;

    tbb::task_group tg;
    {
        unordered_set<int> candidates;
        for (auto it = P.begin(); it != P.end(); it++) {
            if(g.adj[pivot].find(*it) == g.adj[pivot].end())
                candidates.insert(*it);
        }
        for (auto it = candidates.begin(); it != candidates.end(); it++) {
            candidates_vec.push_back(*it);
        }
    }    // free candidates

    tg.run([&] {
    tbb::parallel_for(0, (int)candidates_vec.size(), 1,
        [&](int i) {
            int v = candidates_vec[i];

            tbb::task_group tg1;
            tbb::concurrent_unordered_set<int> P_local;
            tbb::concurrent_unordered_set<int> X_local;

            for (auto it = g.adj[v].begin(); it != g.adj[v].end(); it++) {
                if (P.find(*it) != P.end()) {
                    P_local.insert(*it);
                }
            }
            for (auto it = X.begin(); it != X.end(); it++) {
                if (g.adj[*it].find(v) != g.adj[*it].end())
                    X_local.insert(*it);
            }

            for (int j = 0; j < i; j++) {
                P_local.unsafe_erase(candidates_vec[j]);
                if (g.adj[v].find(candidates_vec[j]) != g.adj[v].end())
                    X_local.insert(candidates_vec[j]);
            }

            parTTT(g, P_local, X_local);
        });
    });
    tg.wait();
}

int parPivot(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X)
{
    int max_size = -1;
    int pivot;
    tbb::concurrent_unordered_map<int, int> neiCnt;

    tbb::task_group tg;

    tg.run([&] {
        tbb::parallel_for_each(X.begin(), X.end(), 
            [&max_size, &g, &pivot, &P, &neiCnt](int v) {
                int count = 0;
                // sequential
                for (auto it = g.adj[v].begin(); it != g.adj[v].end(); it++) {
                    if (P.find(*it) != P.end()) {
                        count++;
                    }
                }
                neiCnt[v] = count;
            });
        });
    tg.wait();
    for (auto it = X.begin(); it != X.end(); it++){
        int v = *it;
        int count = neiCnt[v];
        if (count > max_size) {
            max_size = count;
            pivot = v;
        }
    }

    tg.run([&] {
        tbb::parallel_for_each(P.begin(), P.end(), 
            [&max_size, &g, &pivot, &P, &neiCnt](int v) {
                int count = 0;
                // sequential
                for (auto it = g.adj[v].begin(); it != g.adj[v].end(); it++) {
                    if (P.find(*it) != P.end()) {
                        count++;
                    }
                }
                neiCnt[v] = count;

            });
        });
    tg.wait();

    for (auto it = P.begin(); it != P.end(); it++){
        int v = *it;
        int count = neiCnt[v];
        if (count > max_size) {
            max_size = count;
            pivot = v;
        }
    }

    return pivot;
}

