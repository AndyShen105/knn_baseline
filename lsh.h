//
// Created by Hang Shen on 2018/7/11.
//
#ifndef LSH_H
#define LSH_H
#include <unordered_map>
#include <vector>
#include<queue>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <fstream>
#include <regex>
#include <vector>
#include "data.h"
#include "hashFunc.h"
#include "util.h"
canducate_user calculate_similarity(std::vector<int> seed, int pool_index, int n_feats, float * data, float *queriesint);
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
void pre_user_pool(std::vector<int> &user_bucket_info,
                   float *data,
                   float **hash_func,
                   int n_users,
                   int n_feats,
                   int n_plane);
void gen_ExAudiences_lsh_based(std::priority_queue<canducate_user> &top_k,
                                std::unordered_map<int,std::vector<int>> &user_seed_pool,
                                std::vector<int> user_bucket_info,
                                int n_user,
                                int n_feats,
                                int k,
                                float * data,
                                float *queries);
void pre_user_pool(std::vector<int> &user_bucket_info,
                   float *data,
                   float **hash_func,
                   int n_users,
                   int n_feats,
                   int n_plane);

#endif