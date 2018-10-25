//
// Created by Hang Shen on 2018/7/11.
//
#ifndef LSH_H
#define LSH_H
#include <unordered_map>
#include <vector>
#include<queue>
#include "util.h"
int signature_bit(float *data, float **planes, int index, int n_feats, int n_plane);
float** gen_signature_matrix(int n_feats, int n_plane);
void user_map(std::unordered_map<int, std::vector<int>> &user_maps,
              float *data,
              float **hash_func,
              int n_users,
              int n_feats,
              int n_plane);
void gen_ExAudiences(std::priority_queue<canducate_user> &top_k,
                     std::unordered_map<int, std::vector<int>> user_maps_pool,
                     std::unordered_map<int, std::vector<int>> user_maps_seed,
                     int n_bit,
                     int n_feats,
                     int k,
                     float * data,
                     float *queries);
float get_cosine_dis(int seed_index,
                     int pool_index,
                     int n_feats,
                     float *data,
                     float *queries);
#endif