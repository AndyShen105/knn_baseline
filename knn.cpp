//
// Created by Hang Shen on 2018/3/15.
//

#include <iostream>
#include <iomanip>
#include <math.h>
#include <queue>

#include "knn.h"
#include "util.h"

/*
    data: audience pool
    query: audience seeds
    distance: similarty
    n: size of audience pool
    r: dimensions of audience
    q: numbers of queries
*/

/* Return cosine distance */
inline float query_distance(float *data[], float *query[], int i, int j, int r, int n, int q) {
    float sqDistance = 0.0f;
    float datanorm = 0.0f;
    float querynorm = 0.0f;
    for (int d = 0; d < r; d++) {
        double temp=(*data)[i*r+d] * (*query)[j*r+d];
        datanorm+=(*data)[i*r+d] * (*data)[i*r+d];
        querynorm+=(*query)[j*r+d] * (*query)[j*r+d];
        sqDistance+=temp;
    }
    return sqDistance/(sqrt(datanorm) * sqrt(querynorm));
}

/* Distances array size = n*q */
void serial_distances(float *data[], priority_queue<canducate_user> &top_k, int n, int r, float *query[], int q, int k) {
    for (int i = 0; i < n; i++) { // For each object
        float max_sim = 0.0;
        for (int d = 0; d < q; d++) { // For each query point
            max_sim = max(query_distance(data, query, i, d, r, n, q), max_sim);
        }
        canducate_user temp_user;
            temp_user.sn = i;
            temp_user.sim = query_distance(data, query, i, d, r, n, q);
            if (top_k.size<=k){
                top_k.push(temp_user);
            }
            else{
                if(top_k.top<temp_user.sim){
                    top_k.pop();
                    top_k.push(temp_user);
                }
            }
    }
}

void knn_distance(float *data[], priority_queue<canducate_user> &top_k, int n, int r, float *query[], int q, int k) {
    serial_distances(data, top_k, n, r, query, q, k);

}
