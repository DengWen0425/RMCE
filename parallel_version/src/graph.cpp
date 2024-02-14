#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include "graph.hpp"


Graph::Graph(){

}
Graph::~Graph() {
    if (vals != nullptr) {
        free(vals);
    }
    if (inds != nullptr) {
        free(inds);
    }
    if (seps != nullptr) {
        free(seps);
    }
    if (seps2 != nullptr) {
        free(seps2);
    }
    if (orders != nullptr) {
        free(orders);
    }
    if (o2v != nullptr) {
        free(o2v);
    }
    if (orderedVals != nullptr) {
        free(orderedVals);
    }
}

void Graph::from_path(string const &fileName){
    ifstream instream(fileName.c_str());
    if (instream.good() && !instream.eof())
        instream >> n;
    else {
        fprintf(stderr, "problem with line 1 in input file\n");
        exit(1);
    }
    if (instream.good() && !instream.eof())
        instream >> m;
    else {

        fprintf(stderr, "problem with line 2 in input file\n");
        exit(1);
    }
    vector<vector<int>> adjList(n);
    // vector<list<int>> adjList(n);
    int u, v; // endvertices, to read edges.
    int i = 0;
    while(i < m)
    {
        char comma;
        if (instream.good() && !instream.eof()) {
            instream >> u >> comma >> v;
        } else {
            fprintf(stderr, "problem with line %d in input file\n", i+2);
            exit(1);
        }
        assert(u < n && u > -1);
        assert(v < n && v > -1);
        if(u==v)
            fprintf(stderr, "Detected loop %d->%d\n", u, v);
        assert(u != v);
        adjList[u].emplace_back(v);
        i++;
    }

    int e_idx = 0;
    inds = (int*)malloc(sizeof(int) * (n+1));
    vals = (int*)malloc(sizeof(int) * m);
    // isEdgeDeleted = (int*)malloc(sizeof(int) * m);
    for (int i = 0; i < n; ++i) {
        inds[i] = e_idx;
        sort(adjList[i].begin(), adjList[i].end());
        for (auto u : adjList[i]) {
            // isEdgeDeleted[e_idx] = -1;
            vals[e_idx++] = u;
        }
    }
    inds[n] = m;
    //initialize orders, o2v, seps, orderedVals
    orders = (int*)malloc(sizeof(int) * n);
    o2v = (int*)malloc(sizeof(int) * n);
    seps = (int*)malloc(sizeof(int) * n);
    seps2 = (int*)malloc(sizeof(int) * 1);
    orderedVals = (int*)malloc(sizeof(int) * m);
}

void Graph::computeDegeneracyOrder() {

    vector<NeighborList> vOrdering(n);
    vector<list<int>> verticesByDegree(n);
    vector<list<int>::iterator> vertexLocator(n);
    vector<int> degree(n);

    for(int i=0; i<n; i++)
    {
        degree[i] = inds[i+1] - inds[i];
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }
    int currentDegree = 0;
    int numVerticesRemoved = 0;
    while(numVerticesRemoved < n)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();

            vOrdering[numVerticesRemoved].vertex = vertex;
            vOrdering[numVerticesRemoved].orderNumber = numVerticesRemoved;
            degree[vertex] = -1;

            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor]!=-1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    vOrdering[numVerticesRemoved].later.push_back(neighbor);

                    degree[neighbor]--;

                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
                else
                {
                    vOrdering[numVerticesRemoved].earlier.push_back(neighbor);
                }
            }
            numVerticesRemoved++;
            currentDegree = 0;
        }
        else
        {
            currentDegree++;
        }
    }
    for (int i = 0; i < n; i++){
        int order = vOrdering[i].orderNumber;
        int vertex = vOrdering[i].vertex;
        o2v[order] = vertex;
        orders[vertex] = order;
        int l = inds[vertex];
        for (int const nei:vOrdering[i].earlier){
            orderedVals[l++] = nei;
        }
        seps[vertex] = l;
        for (int const nei:vOrdering[i].later){
            orderedVals[l++] = nei;
        }
    }
}

long Graph::computeDegeneracyWithGlobalReduction(){
    free(seps2);
    seps2 = (int*)malloc(sizeof(int) * n);
    long currentCliques = 0;
    vector<list<int>> verticesByDegree(n);
    vector<list<int>::iterator> vertexLocator(n);
    vector<int> degree(n);
    int* isEdgeDeleted = (int*)malloc(sizeof(int) * m);
    // set all to -1
    fill_n(isEdgeDeleted, m, -1);

    // struct timeval start, mid, end;
    // gettimeofday(&start, NULL);

    for(int i=0; i<n; i++)
    {
        degree[i] = inds[i+1] - inds[i];
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }
    int numVerticesRemoved = 0;

    // non triangle edge removal
    for (int i = 0; i < n; i++)
    {
        if (degree[i] <= 0)
        {
            continue;
        }
        for (int j = inds[i]; j < inds[i+1]; j++)
        {
            int k = vals[j];
            if (i>k || degree[k] <= 0 || isEdgeDeleted[j] != -1)
            {
                continue;
            }
            int u_ind = inds[i], v_ind = inds[k];
            int reverse_j = -1;
            bool canDelete = true;
            while (u_ind < inds[i+1] && v_ind < inds[k+1]) {
                if (vals[v_ind] == i)
                {
                    reverse_j = v_ind;
                }
                if(degree[vals[u_ind]]<=0 || isEdgeDeleted[u_ind] == 1){
                    u_ind++;
                    continue;
                } 
                if(degree[vals[v_ind]]<=0 || isEdgeDeleted[v_ind] == 1){
                    v_ind++;
                    continue;
                } 
                if (vals[u_ind] == vals[v_ind]) {  
                    canDelete = false;  
                    if (isEdgeDeleted[u_ind] == -1 && isEdgeDeleted[v_ind] == -1)
                    {
                        int reverse1 = lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], i) - vals;
                        int reverse2 = lower_bound(vals+reverse1, vals+inds[vals[u_ind]+1], k) - vals;
                        isEdgeDeleted[u_ind] = 0;
                        isEdgeDeleted[v_ind] = 0;
                        isEdgeDeleted[reverse1] = 0;
                        isEdgeDeleted[reverse2] = 0;
                    }
                    else if (isEdgeDeleted[u_ind] == -1){
                        isEdgeDeleted[u_ind] = 0;
                        isEdgeDeleted[lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], k) - vals] = 0;
                    }
                    else if (isEdgeDeleted[v_ind] == -1){
                        isEdgeDeleted[v_ind] = 0;
                        isEdgeDeleted[lower_bound(vals+inds[vals[v_ind]], vals+inds[vals[v_ind]+1], i) - vals] = 0;
                    }
                    v_ind++;
                    u_ind++;
                    break;
                }
                else if (vals[u_ind] < vals[v_ind]) u_ind++;
                else v_ind++;
            }
            if (reverse_j==-1)
            {
                reverse_j = lower_bound(vals+v_ind, vals+inds[k+1], i) - vals; 
            }
            if (canDelete)
            {
                currentCliques++;
                isEdgeDeleted[j]=1;
                isEdgeDeleted[reverse_j]=1;
                verticesByDegree[degree[i]].erase(vertexLocator[i]);
                degree[i]--;
                verticesByDegree[degree[k]].erase(vertexLocator[k]);
                degree[k]--;
                verticesByDegree[degree[i]].push_front(i);
                vertexLocator[i] = verticesByDegree[degree[i]].begin();
                verticesByDegree[degree[k]].push_front(k);
                vertexLocator[k] = verticesByDegree[degree[k]].begin();
            }
            else{
                isEdgeDeleted[j] = 0;
                isEdgeDeleted[reverse_j] = 0;
            }
        }
    }

    // count how many edges are marked as 1, -1, 0 respectively, use i < j format
    int count1 = 0, count_1 = 0, count0 = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = inds[i]; j < inds[i+1]; j++)
        {
            // if (vals[j] < i)
            // {
            //     continue;
            // }
            if (isEdgeDeleted[j] == 1)
            {
                count1++;
            }
            else if (isEdgeDeleted[j] == -1)
            {
                count_1++;
            }
            else
            {
                count0++;
            }
        }
    }
    cout << "count1: " << count1 << endl;
    cout << "count_1: " << count_1 << endl;
    cout << "count0: " << count0 << endl;
    // show deg 0, 1, 2 vertices' number
    cout << "deg 0: " << verticesByDegree[0].size() << endl;
    cout << "deg 1: " << verticesByDegree[1].size() << endl;
    cout << "deg 2: " << verticesByDegree[2].size() << endl;

    // gettimeofday(&mid, NULL);

    cout << "check 2: " << numVerticesRemoved << " " << currentCliques << endl;

    // vertex reduction
    while (!verticesByDegree[1].empty() || !verticesByDegree[2].empty() || !verticesByDegree[0].empty()){

        while (!verticesByDegree[2].empty())
        {
            int const vertex = verticesByDegree[2].front();
            verticesByDegree[2].pop_front();
            int neighbor1=-1, neighbor2=-1;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor] > 0 && isEdgeDeleted[j]!=1)
                {
                    if (neighbor1 == -1)
                    {
                        neighbor1 = neighbor;
                    }
                    else{
                        neighbor2 = neighbor;
                        break;
                    }
                }
            }
            int* p = lower_bound(vals+inds[neighbor1], vals+inds[neighbor1+1], neighbor2);
            bool found = (p != vals+inds[neighbor1+1] && !(neighbor2<*p));
            if (!found)
            {
                verticesByDegree[degree[neighbor1]].erase(vertexLocator[neighbor1]);
                degree[neighbor1]--;
                verticesByDegree[degree[neighbor2]].erase(vertexLocator[neighbor2]);
                degree[neighbor2]--;
                currentCliques = currentCliques + 2;
                verticesByDegree[degree[neighbor1]].push_front(neighbor1);
                vertexLocator[neighbor1] = verticesByDegree[degree[neighbor1]].begin();
                verticesByDegree[degree[neighbor2]].push_front(neighbor2);
                vertexLocator[neighbor2] = verticesByDegree[degree[neighbor2]].begin();
                verticesByDegree[2].erase(vertexLocator[vertex]);
                orders[vertex] = numVerticesRemoved;
                o2v[numVerticesRemoved] = -1;
                degree[vertex] = -2;
                numVerticesRemoved++;
                continue;
            }
            else{ // remove all deg 2 vertices
                verticesByDegree[degree[neighbor1]].erase(vertexLocator[neighbor1]);
                verticesByDegree[degree[neighbor2]].erase(vertexLocator[neighbor2]);
                int u_ind = inds[neighbor1], v_ind = inds[neighbor2];
                bool canDelete = true;
                while (u_ind < inds[neighbor1+1] && v_ind < inds[neighbor2+1]) {
                    if(degree[vals[u_ind]]<=0 || vals[u_ind]==vertex || isEdgeDeleted[u_ind] == 1){
                        u_ind++;
                        continue;
                    } 
                    if(degree[vals[v_ind]]<=0 || vals[v_ind]==vertex || isEdgeDeleted[v_ind] == 1){
                        v_ind++;
                        continue;
                    } 
                    if (vals[u_ind] == vals[v_ind]) {
                        canDelete = false;
                        break;
                    }
                    else if (vals[u_ind] < vals[v_ind]) u_ind++;
                    else v_ind++;
                }
                if(canDelete){
                    isEdgeDeleted[p-vals] = 1;
                    isEdgeDeleted[lower_bound(vals+inds[neighbor2], vals+inds[neighbor2+1], neighbor1) -vals] = 1;
                    degree[neighbor1]--;
                    degree[neighbor2]--;
                }
                degree[neighbor1]--;
                degree[neighbor2]--;
                currentCliques++;
                verticesByDegree[degree[neighbor1]].push_front(neighbor1);
                vertexLocator[neighbor1] = verticesByDegree[degree[neighbor1]].begin();
                verticesByDegree[degree[neighbor2]].push_front(neighbor2);
                vertexLocator[neighbor2] = verticesByDegree[degree[neighbor2]].begin();
                orders[vertex] = numVerticesRemoved;
                o2v[numVerticesRemoved] = -1;
                degree[vertex] = -2;
                numVerticesRemoved++;
                continue;
            }
        }

        while (!verticesByDegree[1].empty())
        {
            currentCliques++;
            int const vertex = verticesByDegree[1].front();
            verticesByDegree[1].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = -1;
            degree[vertex] = -2;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor]!=-2 && isEdgeDeleted[j]!=1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    degree[neighbor]--;
                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }
            numVerticesRemoved++;
        }

        // degree zero removal
        while (!verticesByDegree[0].empty())
        {
            int const vertex = verticesByDegree[0].front();
            verticesByDegree[0].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = -1;
            degree[vertex] = -2;
            numVerticesRemoved++;
        }
    }

    offset = numVerticesRemoved;
    cout << "check 3: " << offset << " currentCliques: " << currentCliques << endl;

    int currentDegree = 0;
    while(numVerticesRemoved < n)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = vertex;
            degree[vertex] = -1;
            int l = inds[vertex], r = inds[vertex+1]-1;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                if(isEdgeDeleted[j] == 1){
                    // orderedVals[l++] = -1;
                    continue;
                }
                int neighbor = vals[j];
                if (degree[neighbor] == -2)
                {
                    continue; // mark deleted vertex as -1, which can not be accessed later
                }
                else if(degree[neighbor] ==-1)
                {
                    orderedVals[l++] = neighbor;
                }
                else
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    orderedVals[r--] = neighbor;
                    degree[neighbor]--;
                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }
            seps[vertex] = l;
            seps2[vertex] = r+1;
            numVerticesRemoved++;
            currentDegree = max(currentDegree-1, 0);
        }
        else
        {
            currentDegree++;
        }
    }
    // gettimeofday(&end, NULL);
    // show time 
    // double time1 = (mid.tv_sec - start.tv_sec) + (mid.tv_usec - start.tv_usec) / 1000000.0;
    // double time2 = (end.tv_sec - mid.tv_sec) + (end.tv_usec - mid.tv_usec) / 1000000.0;
    // cout << "time1: " << time1 << endl;
    // cout << "time2: " << time2 << endl;
    return currentCliques;
}

long Graph::computeDegeneracyWithGlobalReductionParallel(){
    free(seps2);
    seps2 = (int*)malloc(sizeof(int) * n);
    tbb::atomic<long> currentCliques = 0;
    // tbb::concurrent_vector<tbb::concurrent_unordered_set<int>> verticesByDegree(n);
    vector<list<int>> verticesByDegree(n);
    vector<list<int>::iterator> vertexLocator(n);
    tbb::concurrent_vector<tbb::atomic<int>> degree(n);
    tbb::concurrent_vector<int> isEdgeDeleted(m);
    // set all element in isEdgeDeleted to -1
    tbb::task_group tg1;
    tg1.run([&] {
        tbb::parallel_for(0, m, 1, 
        [&](int i) {
            isEdgeDeleted[i] = -1;
        });
    });
    // tg1.wait();

    // tbb::task_group tg2;
    tg1.run([&] {
        tbb::parallel_for(0, n, 1, 
        [&](int i) {
            degree[i] = inds[i+1] - inds[i];
            // verticesByDegree[degree[i]].insert(i);
        });
    });
    tg1.wait();

    // for(int i=0; i<n; i++)
    // {
    //     degree[i] = inds[i+1] - inds[i];
    //     // verticesByDegree[degree[i]].insert(i);
    // }
    int numVerticesRemoved = 0;

    // cout << "check 1: " << n << endl;

    // parallel the non triangle edge removal
    tbb::task_group tg3;
    tg3.run([&] {
        tbb::parallel_for(0, n, 1, 
        [&](int i) {
        // cout << "i: " << i << "-- " << endl;
        if (degree[i] <= 0)
        {
            return;
        }
        tbb::task_group tg4;
        tg4.run([&] {
            tbb::parallel_for(inds[i], inds[i+1], 1, 
                [&](int j) {
                // for (int j = inds[i]; j < inds[i+1]; j++){
                        int k = vals[j];
                        if (i>k || degree[k] <= 0 || isEdgeDeleted[j] != -1)
                        {
                            return;
                        }
                        // cout << "i, k: " << i << " " << k << endl;
                        int u_ind = inds[i], v_ind = inds[k];
                        int reverse_j = -1;
                        bool canDelete = true;
                        while (u_ind < inds[i+1] && v_ind < inds[k+1]) {
                            if (vals[v_ind] == i)
                            {
                                reverse_j = v_ind;
                            }
                            if(degree[vals[u_ind]]<=0 || isEdgeDeleted[u_ind] == 1){
                                u_ind++;
                                continue;
                            } 
                            if(degree[vals[v_ind]]<=0 || isEdgeDeleted[v_ind] == 1){
                                v_ind++;
                                continue;
                            } 
                            if (vals[u_ind] == vals[v_ind]) {  
                                canDelete = false;  
                                if (isEdgeDeleted[u_ind] == -1 && isEdgeDeleted[v_ind] == -1)
                                {
                                    int reverse1 = lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], i) - vals;
                                    int reverse2 = lower_bound(vals+reverse1, vals+inds[vals[u_ind]+1], k) - vals;
                                    isEdgeDeleted[u_ind] = 0;
                                    isEdgeDeleted[v_ind] = 0;
                                    isEdgeDeleted[reverse1] = 0;
                                    isEdgeDeleted[reverse2] = 0;
                                }
                                else if (isEdgeDeleted[u_ind] == -1){
                                    isEdgeDeleted[u_ind] = 0;
                                    isEdgeDeleted[lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], k) - vals] = 0;
                                }
                                else if (isEdgeDeleted[v_ind] == -1){
                                    isEdgeDeleted[v_ind] = 0;
                                    isEdgeDeleted[lower_bound(vals+inds[vals[v_ind]], vals+inds[vals[v_ind]+1], i) - vals] = 0;
                                }
                                v_ind++;
                                u_ind++;
                                break;
                            }
                            else if (vals[u_ind] < vals[v_ind]) u_ind++;
                            else v_ind++;
                        }
                        if (reverse_j==-1)
                        {
                            reverse_j = lower_bound(vals+v_ind, vals+inds[k+1], i) - vals; 
                        }
                        if (canDelete)
                        {
                            currentCliques++;
                            isEdgeDeleted[j]=1;
                            isEdgeDeleted[reverse_j]=1;
                            degree[i]--;
                            degree[k]--;
                        }
                        else{
                            isEdgeDeleted[j] = 0;
                            isEdgeDeleted[reverse_j] = 0;
                        }
                });
            });
            tg4.wait();
        });
    });
    tg3.wait();

    // count how many edges are marked as 1, -1, 0 respectively
    // int count1 = 0, count_1 = 0, count0 = 0;
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = inds[i]; j < inds[i+1]; j++)
    //     {
    //         // if (vals[j] < i)
    //         // {
    //         //     continue;
    //         // }
    //         if (isEdgeDeleted[j] == 1)
    //         {
    //             count1++;
    //             // check whether i is marked 1 in vals[j]
    //             int k = vals[j];
    //             int reverse_j = lower_bound(vals+inds[k], vals+inds[k+1], i) - vals;
    //             if (isEdgeDeleted[reverse_j] != 1)
    //             {
    //                 cout << "i, k: " << i << " " << k << endl;
    //             }
    //         }
    //         else if (isEdgeDeleted[j] == -1)
    //         {
    //             count_1++;
    //         }
    //         else
    //         {
    //             count0++;
    //         }
    //         // cout << " i, j, isEdgeDeleted[j]: " << i << " " << vals[j] << " " << isEdgeDeleted[j] << endl;
    //     }
    // }
    // cout << "count1: " << count1 << endl;
    // cout << "count_1: " << count_1 << endl;
    // cout << "count0: " << count0 << endl;

    // build verticesByDegree and vertexLocator in sequential
    for(int i=0; i<n; i++)
    {
        verticesByDegree[degree[i]].push_front(i);
        vertexLocator[i] = verticesByDegree[degree[i]].begin();
    }
    // show deg 0, 1, 2 vertices' number
    // cout << "deg 0: " << verticesByDegree[0].size() << endl;
    // cout << "deg 1: " << verticesByDegree[1].size() << endl;
    // cout << "deg 2: " << verticesByDegree[2].size() << endl;

    // cout << "check 2: " << numVerticesRemoved << " " << currentCliques << endl;

    while (!verticesByDegree[1].empty() || !verticesByDegree[2].empty() || !verticesByDegree[0].empty()){
        // cout << "check 2.1: " << numVerticesRemoved << " " << currentCliques << endl;
        while (!verticesByDegree[2].empty())
        {
            int const vertex = verticesByDegree[2].front();
            verticesByDegree[2].pop_front();
            int neighbor1=-1, neighbor2=-1;
            // cout << "vertex: " << vertex << endl;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor] > 0 && isEdgeDeleted[j]!=1)
                {
                    if (neighbor1 == -1)
                    {
                        neighbor1 = neighbor;
                    }
                    else{
                        neighbor2 = neighbor;
                        break;
                    }
                }
            }
            // cout << "neighbor1, neighbor2: " << neighbor1 << " " << neighbor2 << endl;
            int* p = lower_bound(vals+inds[neighbor1], vals+inds[neighbor1+1], neighbor2);
            bool found = (p != vals+inds[neighbor1+1] && !(neighbor2<*p));
            if (!found)
            {
                verticesByDegree[degree[neighbor1]].erase(vertexLocator[neighbor1]);
                degree[neighbor1]--;
                verticesByDegree[degree[neighbor2]].erase(vertexLocator[neighbor2]);
                degree[neighbor2]--;
                currentCliques = currentCliques + 2;
                verticesByDegree[degree[neighbor1]].push_front(neighbor1);
                vertexLocator[neighbor1] = verticesByDegree[degree[neighbor1]].begin();
                verticesByDegree[degree[neighbor2]].push_front(neighbor2);
                vertexLocator[neighbor2] = verticesByDegree[degree[neighbor2]].begin();
                verticesByDegree[2].erase(vertexLocator[vertex]);
                orders[vertex] = numVerticesRemoved;
                o2v[numVerticesRemoved] = -1;
                degree[vertex] = -2;
                numVerticesRemoved++;
                continue;
            }
            else{ // remove all deg 2 vertices
                verticesByDegree[degree[neighbor1]].erase(vertexLocator[neighbor1]);
                verticesByDegree[degree[neighbor2]].erase(vertexLocator[neighbor2]);
                int u_ind = inds[neighbor1], v_ind = inds[neighbor2];
                bool canDelete = true;
                while (u_ind < inds[neighbor1+1] && v_ind < inds[neighbor2+1]) {
                    if(degree[vals[u_ind]]<=0 || vals[u_ind]==vertex || isEdgeDeleted[u_ind] == 1){
                        u_ind++;
                        continue;
                    } 
                    if(degree[vals[v_ind]]<=0 || vals[v_ind]==vertex || isEdgeDeleted[v_ind] == 1){
                        v_ind++;
                        continue;
                    } 
                    if (vals[u_ind] == vals[v_ind]) {
                        canDelete = false;
                        break;
                    }
                    else if (vals[u_ind] < vals[v_ind]) u_ind++;
                    else v_ind++;
                }
                if(canDelete){
                    isEdgeDeleted[p-vals] = 1;
                    isEdgeDeleted[lower_bound(vals+inds[neighbor2], vals+inds[neighbor2+1], neighbor1) -vals] = 1;
                    degree[neighbor1]--;
                    degree[neighbor2]--;
                }
                degree[neighbor1]--;
                degree[neighbor2]--;
                currentCliques++;
                verticesByDegree[degree[neighbor1]].push_front(neighbor1);
                vertexLocator[neighbor1] = verticesByDegree[degree[neighbor1]].begin();
                verticesByDegree[degree[neighbor2]].push_front(neighbor2);
                vertexLocator[neighbor2] = verticesByDegree[degree[neighbor2]].begin();
                orders[vertex] = numVerticesRemoved;
                o2v[numVerticesRemoved] = -1;
                degree[vertex] = -2;
                numVerticesRemoved++;
                continue;
            }
        }
        // cout << "check 2.2: " << numVerticesRemoved << " " << currentCliques << endl;
        while (!verticesByDegree[1].empty())
        {
            currentCliques++;
            int const vertex = verticesByDegree[1].front();
            verticesByDegree[1].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = -1;
            degree[vertex] = -2;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                int neighbor = vals[j];
                if(degree[neighbor]!=-2 && isEdgeDeleted[j]!=1)
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    degree[neighbor]--;
                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }
            numVerticesRemoved++;
        }

        // cout << "check 2.3: " << numVerticesRemoved << " " << currentCliques << endl;
        while (!verticesByDegree[0].empty())
        {
            int const vertex = verticesByDegree[0].front();
            verticesByDegree[0].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = -1;
            degree[vertex] = -2;
            numVerticesRemoved++;
        }
    }

    offset = numVerticesRemoved;
    // cout << "check 3: " << offset << " currentCliques: " << currentCliques << endl;

    int currentDegree = 0;
    while(numVerticesRemoved < n)
    {
        if(!verticesByDegree[currentDegree].empty())
        {
            degeneracy = max(degeneracy,currentDegree);
            int const vertex = verticesByDegree[currentDegree].front();
            verticesByDegree[currentDegree].pop_front();
            orders[vertex] = numVerticesRemoved;
            o2v[numVerticesRemoved] = vertex;
            degree[vertex] = -1;
            int l = inds[vertex], r = inds[vertex+1]-1;
            for(int j=inds[vertex]; j < inds[vertex+1]; j++)
            {
                if(isEdgeDeleted[j] == 1){
                    // orderedVals[l++] = -1;
                    continue;
                }
                int neighbor = vals[j];
                if (degree[neighbor] == -2)
                {
                    continue; // mark deleted vertex as -1, which can not be accessed later
                }
                else if(degree[neighbor] ==-1)
                {
                    orderedVals[l++] = neighbor;
                }
                else
                {
                    verticesByDegree[degree[neighbor]].erase(vertexLocator[neighbor]);
                    orderedVals[r--] = neighbor;
                    degree[neighbor]--;
                    if(degree[neighbor] != -1)
                    {
                        verticesByDegree[degree[neighbor]].push_front(neighbor);
                        vertexLocator[neighbor] = verticesByDegree[degree[neighbor]].begin();
                    }
                }
            }
            seps[vertex] = l;
            seps2[vertex] = r+1;
            numVerticesRemoved++;
            currentDegree = max(currentDegree-1, 0);
        }
        else
        {
            currentDegree++;
        }
    }
    // gettimeofday(&end, NULL);
    // show time 
    // double time1 = (mid.tv_sec - start.tv_sec) + (mid.tv_usec - start.tv_usec) / 1000000.0;
    // double time2 = (end.tv_sec - mid.tv_sec) + (end.tv_usec - mid.tv_usec) / 1000000.0;
    // cout << "time1: " << time1 << endl;
    // cout << "time2: " << time2 << endl;
    // cout << "check 4 "  << endl;
    return currentCliques;
}

// long Graph::computeDegeneracyWithGlobalReductionParallel(tbb::task_group &tg){
//     free(seps2);
//     seps2 = (int*)malloc(sizeof(int) * n);
//     tbb::atomic<long> currentCliques = 0;
//     tbb::concurrent_vector<tbb::concurrent_unordered_set<int>> verticesByDegree(n);
//     // vector<list<int>> verticesByDegree(n);
//     // vector<list<int>::iterator> vertexLocator(n);
//     tbb::concurrent_vector<tbb::atomic<int>> degree(n);
//     tbb::concurrent_vector<int> isEdgeDeleted(m);
//     // set all element in isEdgeDeleted to -1
//     tg.run([&] {
//         tbb::parallel_for(0, m, 1, 
//         [&](int i) {
//             isEdgeDeleted[i] = -1;
//         });
//     });
//     tg.wait();
//     for(int i=0; i<n; i++)
//     {
//         degree[i] = inds[i+1] - inds[i];
//         // verticesByDegree[degree[i]].insert(i);
//     }
//     int numVerticesRemoved = 0;

//     // cout << "check 1: " << n << endl;

//     // parallel the non triangle edge removal
//     tg.run([&] {
//         tbb::parallel_for(0, n, 1, 
//         [&](int i) {
//         // cout << "i: " << i << "-- " << endl;
//         if (degree[i] <= 0)
//         {
//             return;
//         }
//         tbb::task_group tg1;
//         tg1.run([&] {
//             tbb::parallel_for(inds[i], inds[i+1], 1, 
//                 [&](int j) {
//                 // for (int j = inds[i]; j < inds[i+1]; j++){
//                         int k = vals[j];
//                         if (i>k || degree[k] <= 0 || isEdgeDeleted[j] != -1)
//                         {
//                             return;
//                         }
//                         // cout << "i, k: " << i << " " << k << endl;
//                         int u_ind = inds[i], v_ind = inds[k];
//                         int reverse_j = -1;
//                         bool canDelete = true;
//                         while (u_ind < inds[i+1] && v_ind < inds[k+1]) {
//                             if (vals[v_ind] == i)
//                             {
//                                 reverse_j = v_ind;
//                             }
//                             if(degree[vals[u_ind]]<=0 || isEdgeDeleted[u_ind] == 1){
//                                 u_ind++;
//                                 continue;
//                             } 
//                             if(degree[vals[v_ind]]<=0 || isEdgeDeleted[v_ind] == 1){
//                                 v_ind++;
//                                 continue;
//                             } 
//                             if (vals[u_ind] == vals[v_ind]) {  
//                                 canDelete = false;  
//                                 if (isEdgeDeleted[u_ind] == -1 && isEdgeDeleted[v_ind] == -1)
//                                 {
//                                     int reverse1 = lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], i) - vals;
//                                     int reverse2 = lower_bound(vals+reverse1, vals+inds[vals[u_ind]+1], k) - vals;
//                                     isEdgeDeleted[u_ind] = 0;
//                                     isEdgeDeleted[v_ind] = 0;
//                                     isEdgeDeleted[reverse1] = 0;
//                                     isEdgeDeleted[reverse2] = 0;
//                                 }
//                                 else if (isEdgeDeleted[u_ind] == -1){
//                                     isEdgeDeleted[u_ind] = 0;
//                                     isEdgeDeleted[lower_bound(vals+inds[vals[u_ind]], vals+inds[vals[u_ind]+1], k) - vals] = 0;
//                                 }
//                                 else if (isEdgeDeleted[v_ind] == -1){
//                                     isEdgeDeleted[v_ind] = 0;
//                                     isEdgeDeleted[lower_bound(vals+inds[vals[v_ind]], vals+inds[vals[v_ind]+1], i) - vals] = 0;
//                                 }
//                                 v_ind++;
//                                 u_ind++;
//                                 break;
//                             }
//                             else if (vals[u_ind] < vals[v_ind]) u_ind++;
//                             else v_ind++;
//                         }
//                         if (reverse_j==-1)
//                         {
//                             reverse_j = lower_bound(vals+v_ind, vals+inds[k+1], i) - vals; 
//                         }
//                         if (canDelete)
//                         {
//                             currentCliques++;
//                             isEdgeDeleted[j]=1;
//                             isEdgeDeleted[reverse_j]=1;
//                             if (degree[i] == 0 || degree[k] == 0)
//                             {
//                                 cout << "i, k: " << i << " " << k << endl;
//                             }
//                             degree[i]--;
//                             degree[k]--;
//                         }
//                         else{
//                             isEdgeDeleted[j] = 0;
//                             isEdgeDeleted[reverse_j] = 0;
//                         }
//                     // }
//                 });
//             });
//             tg1.wait();
//         });
//     });
//     tg.wait();

//     // count how many edges are marked as 1, -1, 0 respectively
//     // int count1 = 0, count_1 = 0, count0 = 0;
//     // for (int i = 0; i < n; i++)
//     // {
//     //     for (int j = inds[i]; j < inds[i+1]; j++)
//     //     {
//     //         // if (vals[j] < i)
//     //         // {
//     //         //     continue;
//     //         // }
//     //         if (isEdgeDeleted[j] == 1)
//     //         {
//     //             count1++;
//     //             // check whether i is marked 1 in vals[j]
//     //             int k = vals[j];
//     //             int reverse_j = lower_bound(vals+inds[k], vals+inds[k+1], i) - vals;
//     //             if (isEdgeDeleted[reverse_j] != 1)
//     //             {
//     //                 cout << "i, k: " << i << " " << k << endl;
//     //             }
//     //         }
//     //         else if (isEdgeDeleted[j] == -1)
//     //         {
//     //             count_1++;
//     //         }
//     //         else
//     //         {
//     //             count0++;
//     //         }
//     //         // cout << " i, j, isEdgeDeleted[j]: " << i << " " << vals[j] << " " << isEdgeDeleted[j] << endl;
//     //     }
//     // }
//     // cout << "count1: " << count1 << endl;
//     // cout << "count_1: " << count_1 << endl;
//     // cout << "count0: " << count0 << endl;

//     // build verticesByDegree and vertexLocator
//     tg.run([&] {
//         tbb::parallel_for(0, n, 1, 
//         [&](int i) {
//             if (degree[i] > 0)
//             {
//                 verticesByDegree[degree[i]].insert(i);
//             }
//         });
//     });
//     tg.wait();


//     while (!verticesByDegree[1].empty() || !verticesByDegree[2].empty() || !verticesByDegree[0].empty()){
//         // cout << "check 2.1: " << numVerticesRemoved << " " << currentCliques << endl;
//         while (!verticesByDegree[2].empty())
//         {
//             int const vertex = *verticesByDegree[2].begin();
//             verticesByDegree[2].unsafe_erase(vertex);
//             int neighbor1=-1, neighbor2=-1;
//             // cout << "vertex: " << vertex << endl;
//             for(int j=inds[vertex]; j < inds[vertex+1]; j++)
//             {
//                 int neighbor = vals[j];
//                 if(degree[neighbor] > 0 && isEdgeDeleted[j]!=1)
//                 {
//                     if (neighbor1 == -1)
//                     {
//                         neighbor1 = neighbor;
//                     }
//                     else{
//                         neighbor2 = neighbor;
//                         break;
//                     }
//                 }
//             }
//             // cout << "neighbor1, neighbor2: " << neighbor1 << " " << neighbor2 << endl;
//             int* p = lower_bound(vals+inds[neighbor1], vals+inds[neighbor1+1], neighbor2);
//             bool found = (p != vals+inds[neighbor1+1] && !(neighbor2<*p));
//             if (!found)
//             {
//                 verticesByDegree[degree[neighbor1]].unsafe_erase(neighbor1);
//                 degree[neighbor1]--;
//                 verticesByDegree[degree[neighbor2]].unsafe_erase(neighbor2);
//                 degree[neighbor2]--;
//                 currentCliques = currentCliques + 2;
//                 verticesByDegree[degree[neighbor1]].insert(neighbor1);
//                 verticesByDegree[degree[neighbor2]].insert(neighbor2);
//                 verticesByDegree[2].unsafe_erase(vertex);
//                 orders[vertex] = numVerticesRemoved;
//                 o2v[numVerticesRemoved] = -1;
//                 degree[vertex] = -2;
//                 numVerticesRemoved++;
//                 continue;
//             }
//             else{ // remove all deg 2 vertices
//                 verticesByDegree[degree[neighbor1]].unsafe_erase(neighbor1);
//                 verticesByDegree[degree[neighbor2]].unsafe_erase(neighbor2);
//                 int u_ind = inds[neighbor1], v_ind = inds[neighbor2];
//                 bool canDelete = true;
//                 while (u_ind < inds[neighbor1+1] && v_ind < inds[neighbor2+1]) {
//                     if(degree[vals[u_ind]]<=0 || vals[u_ind]==vertex || isEdgeDeleted[u_ind] == 1){
//                         u_ind++;
//                         continue;
//                     } 
//                     if(degree[vals[v_ind]]<=0 || vals[v_ind]==vertex || isEdgeDeleted[v_ind] == 1){
//                         v_ind++;
//                         continue;
//                     } 
//                     if (vals[u_ind] == vals[v_ind]) {
//                         canDelete = false;
//                         break;
//                     }
//                     else if (vals[u_ind] < vals[v_ind]) u_ind++;
//                     else v_ind++;
//                 }
//                 if(canDelete){
//                     isEdgeDeleted[p-vals] = 1;
//                     isEdgeDeleted[lower_bound(vals+inds[neighbor2], vals+inds[neighbor2+1], neighbor1) -vals] = 1;
//                     degree[neighbor1]--;
//                     degree[neighbor2]--;
//                 }
//                 degree[neighbor1]--;
//                 degree[neighbor2]--;
//                 currentCliques++;
//                 verticesByDegree[degree[neighbor1]].insert(neighbor1);
//                 verticesByDegree[degree[neighbor2]].insert(neighbor2);
//                 orders[vertex] = numVerticesRemoved;
//                 o2v[numVerticesRemoved] = -1;
//                 degree[vertex] = -2;
//                 numVerticesRemoved++;
//                 continue;
//             }
//         }
//         // cout << "check 2.2: " << numVerticesRemoved << " " << currentCliques << endl;
//         while (!verticesByDegree[1].empty())
//         {
//             currentCliques++;
//             int const vertex = *verticesByDegree[1].begin();
//             verticesByDegree[1].unsafe_erase(vertex);
//             orders[vertex] = numVerticesRemoved;
//             o2v[numVerticesRemoved] = -1;
//             degree[vertex] = -2;
//             for(int j=inds[vertex]; j < inds[vertex+1]; j++)
//             {
//                 int neighbor = vals[j];
//                 if(degree[neighbor]!=-2 && isEdgeDeleted[j]!=1)
//                 {
//                     verticesByDegree[degree[neighbor]].unsafe_erase(neighbor);
//                     degree[neighbor]--;
//                     if(degree[neighbor] != -1)
//                     {
//                         verticesByDegree[degree[neighbor]].insert(neighbor);
//                     }
//                 }
//             }
//             numVerticesRemoved++;
//         }

//         // cout << "check 2.3: " << numVerticesRemoved << " " << currentCliques << endl;
//         while (!verticesByDegree[0].empty())
//         {
//             int const vertex = *verticesByDegree[0].begin();
//             verticesByDegree[0].unsafe_erase(vertex);
//             orders[vertex] = numVerticesRemoved;
//             o2v[numVerticesRemoved] = -1;
//             degree[vertex] = -2;
//             numVerticesRemoved++;
//         }
//     }

//     offset = numVerticesRemoved;
//     // cout << "check 3: " << offset << " currentCliques: " << currentCliques << endl;

//     int currentDegree = 0;
//     while(numVerticesRemoved < n)
//     {
//         if(!verticesByDegree[currentDegree].empty())
//         {
//             degeneracy = max(degeneracy,currentDegree);
//             int const vertex = *verticesByDegree[currentDegree].begin();
//             verticesByDegree[currentDegree].unsafe_erase(vertex);
//             orders[vertex] = numVerticesRemoved;
//             o2v[numVerticesRemoved] = vertex;
//             degree[vertex] = -1;
//             int l = inds[vertex], r = inds[vertex+1]-1;
//             for(int j=inds[vertex]; j < inds[vertex+1]; j++)
//             {
//                 if(isEdgeDeleted[j] == 1){
//                     // orderedVals[l++] = -1;
//                     continue;
//                 }
//                 int neighbor = vals[j];
//                 if (degree[neighbor] == -2)
//                 {
//                     continue; // mark deleted vertex as -1, which can not be accessed later
//                 }
//                 else if(degree[neighbor] ==-1)
//                 {
//                     orderedVals[l++] = neighbor;
//                 }
//                 else
//                 {
//                     verticesByDegree[degree[neighbor]].unsafe_erase(neighbor);
//                     orderedVals[r--] = neighbor;
//                     degree[neighbor]--;
//                     if(degree[neighbor] != -1)
//                     {
//                         verticesByDegree[degree[neighbor]].insert(neighbor);
//                     }
//                 }
//             }
//             seps[vertex] = l;
//             seps2[vertex] = r+1;
//             numVerticesRemoved++;
//             currentDegree = max(currentDegree-1, 0);
//         }
//         else
//         {
//             currentDegree++;
//         }
//     }
//     // gettimeofday(&end, NULL);
//     // show time 
//     // double time1 = (mid.tv_sec - start.tv_sec) + (mid.tv_usec - start.tv_usec) / 1000000.0;
//     // double time2 = (end.tv_sec - mid.tv_sec) + (end.tv_usec - mid.tv_usec) / 1000000.0;
//     // cout << "time1: " << time1 << endl;
//     // cout << "time2: " << time2 << endl;
//     // cout << "check 4 "  << endl;
//     return currentCliques;
// }

SubGraph::SubGraph(){}
SubGraph::~SubGraph(){

}

void SubGraph::from_graph(Graph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X){
    adj.clear();
    tbb::task_group tg;
    tg.run([&] {
    tbb::parallel_for_each(X.begin(), X.end(),
        [&](int v) {
            tbb::task_group tg1;
            tg1.run([&] {
                tbb::parallel_for(g.seps[v], g.inds[v+1], 1, 
                [&](int i) {
                    if(P.find(g.orderedVals[i])!=P.end()){
                        adj[v].insert(g.orderedVals[i]);
                    }
                });
            });
            tg1.wait();
        });
    });
    tg.wait();
    tg.run([&] {
    tbb::parallel_for_each(P.begin(), P.end(),
        [&](int v) {
            tbb::task_group tg1;
            tg1.run([&] {
                tbb::parallel_for(g.seps[v], g.inds[v+1], 1, 
                [&](int i) {
                    if(P.find(g.orderedVals[i])!=P.end()){
                        adj[v].insert(g.orderedVals[i]);
                        adj[g.orderedVals[i]].insert(v);
                    }
                });
            });
            tg1.wait();
        });
    });
    tg.wait();
}

void SubGraph::from_graph_reduced(Graph &g, tbb::concurrent_unordered_set<int> &P, tbb::concurrent_unordered_set<int> &X){
    adj.clear();
    tbb::task_group tg;
    tg.run([&] {
    tbb::parallel_for_each(X.begin(), X.end(),
        [&](int v) {
            //sequential
            // for(int i=g.seps[v]; i < g.inds[v+1]; i++){
            //     if(P.find(g.orderedVals[i])!=P.end()){
            //         adj[v].insert(g.orderedVals[i]);
            //     }
            // }
            tbb::task_group tg1;
            tg1.run([&] {
                tbb::parallel_for(g.seps2[v], g.inds[v+1], 1, 
                [&](int i) {
                    if(P.find(g.orderedVals[i])!=P.end()){
                        adj[v].insert(g.orderedVals[i]);
                        // adjInX[g.orderedVals[i]].insert(v);
                    }
                });
            });
            tg1.wait();
        });
    });
    tg.wait();
    tg.run([&] {
    tbb::parallel_for_each(P.begin(), P.end(),
        [&](int v) {
            // sequential 
            // for (int i = g.seps2[v]; i < g.inds[v+1]; i++){
            //     if(P.find(g.orderedVals[i])!=P.end()){
            //         adj[v].insert(g.orderedVals[i]);
            //         adj[g.orderedVals[i]].insert(v);
            //     }
            // }
            tbb::task_group tg1;
            tg1.run([&] {
                tbb::parallel_for(g.seps2[v], g.inds[v+1], 1, 
                [&](int i) {
                    if(P.find(g.orderedVals[i])!=P.end()){
                        adj[v].insert(g.orderedVals[i]);
                        adj[g.orderedVals[i]].insert(v);
                    }
                });
            });
            tg1.wait();
        });
    });
    tg.wait();
    // set every vertex in P in variable bfsParent to false in sequential
    // for (auto it = P.begin(); it != P.end(); it++) {
    //     bfsParent[*it] = false;
    // }
}

