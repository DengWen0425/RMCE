#include <iostream>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_hash_map.h>
#include <sys/time.h>
#include <atomic>
#include <thread>
#include <assert.h>
#include <chrono>

#include "reducedparmce.hpp"


void reduced_parmce_main(string file_name, int num_thread)
{
    tbb::task_scheduler_init init(num_thread);
    cout << "Allocated thread number: " << num_thread << endl;

    /* Initilize CSR */
    Graph g;
    g.from_path(file_name);
    tbb::concurrent_unordered_set<int> P;
    tbb::concurrent_unordered_set<int> X;

    tbb::task_group tg0;

    std::chrono::steady_clock::time_point start, mid, end;
    start  = std::chrono::steady_clock::now();

    cliqueCount = cliqueCount + g.computeDegeneracyWithGlobalReductionParallel();

    mid = std::chrono::steady_clock::now();

    /* Use tbb::global control */
    tbb::blocked_range<int> range(g.offset, g.n);
    tg0.run([&] {
        tbb::parallel_for(range,
            [&](const tbb::blocked_range<int>& r) {
                for (int i = r.begin(); i != r.end(); i++) {
                    int v = g.o2v[i];
                    tbb::concurrent_unordered_set<int> P_local;
                    tbb::concurrent_unordered_set<int> X_local;

                    tbb::task_group tg;
                    tg.run([&] {
                    tbb::parallel_for(g.inds[v], g.seps[v], 1,
                        [&](int i) {
                            X_local.insert(g.orderedVals[i]);
                        });
                    });
                    // tg.wait();

                    tg.run([&] {
                    tbb::parallel_for(g.seps2[v], g.inds[v+1], 1,
                        [&](int i) {
                            P_local.insert(g.orderedVals[i]);
                        });
                    });
                    tg.wait();

                    SubGraph subg;
                    subg.from_graph_reduced(g, P_local, X_local);
                    subg.anchor = v;
                    reduced_parTTT(subg, P_local, X_local);
                }
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

void reduced_parTTT(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X)
{
    if (P.empty() && X.empty()) {
        cliqueCount++;
        return;
    }

    int pivot = reduced_parPivot(g, P, X);
    if (pivot < 0) {
        return;
    }
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

            reduced_parTTT(g, P_local, X_local);
        });
    });
    tg.wait();
}

int reduced_parPivot(SubGraph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X)
{
    tbb::concurrent_unordered_map<int, int> neiCnt;
    int max_size = -1;
    int pivot = -3;

    tbb::task_group tg;

    tg.run([&] {
        tbb::parallel_for_each(X.begin(), X.end(), 
            [&](int v) {
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
        if (count == (int)P.size()) {
            max_size = count;
            pivot = -1;
            break;
        }
        if (count > max_size) {
            max_size = count;
            pivot = v;
        }
    }

    if (pivot == -1) {
        return pivot;
    }

    int ori_p_size = P.size();
    tbb::concurrent_unordered_set<int> toDelete;
    tbb::atomic<int> eCnt=0;
    tg.run([&] {
        tbb::parallel_for_each(P.begin(), P.end(), 
            [&](int v) {
                int count = 0;
                // sequential
                for (auto it = g.adj[v].begin(); it != g.adj[v].end(); it++) {
                    if (P.find(*it) != P.end()) {
                        count++;
                    }
                }
                if (count == 0){
                    bool find = false;
                    for (auto it = X.begin(); it != X.end(); it++) {
                        if (g.adj[*it].find(v) != g.adj[*it].end()){
                            find = true;
                            break;
                        }
                    }
                    if (!find){
                        cliqueCount++;
                    }
                    toDelete.insert(v);
                }
                else{
                    eCnt = eCnt + count;
                    neiCnt[v] = count;
                }
            });
        });
    tg.wait();
    
    if (toDelete.size() == 0 && eCnt == (int)P.size() * ((int)P.size()-1)) {
        cliqueCount++;
        return -1;
    }
    for (auto it = toDelete.begin(); it != toDelete.end(); it++) {
        P.unsafe_erase(*it);
    }
    if (eCnt == 0) {
        return -1;
    }

    for (auto it = P.begin(); it != P.end(); it++){
        int v = *it;
        int count = neiCnt[v];
        if (count > max_size) {
            max_size = count;
            pivot = v;
        }
    }

    toDelete.clear();
    tg.run([&] {
        tbb::parallel_for_each(P.begin(), P.end(), 
            [&](int v) {
                int count = neiCnt[v];
                if (count == (int)P.size()-1) {
                    toDelete.insert(v);
                }
            });
        });
    tg.wait();
    for (auto it = toDelete.begin(); it != toDelete.end(); it++) {
        P.unsafe_erase(*it);
    }

    if ((int)P.size() != ori_p_size) {
        // recompute the pivot 
        bool recomputePivot = (int)toDelete.size() > 0;
        if (recomputePivot){
            max_size = -1;
            pivot = -3;
        }
        tbb::concurrent_unordered_set<int> toDelete2;
        tg.run([&] {
            tbb::parallel_for_each(X.begin(), X.end(), 
                [&](int v) {
                    int count = 0;
                    for (auto it = toDelete.begin(); it != toDelete.end(); it++) {
                        if (g.adj[v].find(*it) != g.adj[v].end()){
                            count++;
                        }
                    }
                    if (count != (int)toDelete.size()){
                        toDelete2.insert(v);
                    }
                    else{
                        if (recomputePivot){
                            int count = 0;
                            // sequential
                            for (auto it = g.adj[v].begin(); it != g.adj[v].end(); it++) {
                                if (P.find(*it) != P.end()) {
                                    count++;
                                }
                            }
                            neiCnt[v] = count;
                        }
                    }
            });
        });
        tg.wait();
        for (auto it = toDelete2.begin(); it != toDelete2.end(); it++) {
            X.unsafe_erase(*it);
        }
        toDelete2.clear();

        for (auto it = X.begin(); it != X.end(); it++){
            int v = *it;
            int count = neiCnt[v];
            if (count > max_size) {
                max_size = count;
                pivot = v;
            }
        }

        if (P.size() == 0) {
            toDelete.clear();
            if (X.size() == 0){
                cliqueCount++;
            }
            return pivot;
        }

        if (pivot == -1) {
            return pivot;
        }
        if (recomputePivot){
            tg.run([&] {
                tbb::parallel_for_each(P.begin(), P.end(), 
                    [&](int v) {
                        int count = 0;
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
        }
    }
    toDelete.clear();
    return pivot;
}

