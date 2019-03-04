//
// Created by andyshen on 19-3-4.
//

#ifndef KNN_CUDAOPT_H
#define KNN_CUDAOPT_H

#include <unordered_map>
#include <vector>
#include "data.h"
#include "knn.h"
#include "sort.h"
#include "util.h"
#include "lsh.h"
#include "v_lsh.h"
#include "hashFunc.h"

using namespace std;
void gen_ExAudiences_cudabase(priority_queue<canducate_user> &top_k,
                              unordered_map<int,vector<int>> user_maps_seed,
                              vector<int> user_bucket_info,
                              priority_queue<uncertain_user> &indexUser,
                              vector<bucket_info> centroid_angle,
                              int n_bit,
                              int n_feats,
                              int k,
                              int n,
                              int q,
                              float * data,
                              float *queries);

void gen_ExAudiences_cuda(priority_queue<canducate_user> &top_k,
                          unordered_map<int,vector<int>> user_maps_pool,
                          unordered_map<int,vector<int>> user_maps_seed,
                          priority_queue<uncertain_user> &indexUser,
                          vector<bucket_info> centroid_angle,
                          int n_bit,
                          int n_feats,
                          int k,
                          int n,
                          int q,
                          float * data,
                          float *queries);

void gen_ExAudiences_cudaOpt(priority_queue<canducate_user> &top_k,
                             unordered_map<int,vector<int>> user_maps_pool,
                             unordered_map<int,vector<int>> user_maps_seed,
                             priority_queue<uncertain_user> &indexUser,
                             vector<bucket_info> centroid_angle,
                             int n_bit,
                             int n_feats,
                             int k,
                             int n,
                             int q,
                             float * data,
                             float *queries);
#endif //KNN_CUDAOPT_H
