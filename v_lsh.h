//
// Created by Hang Shen on 2018/7/11.
//
#ifndef V_LSH_H
#define V_LSH_H
#include <unordered_map>
#include "util.h"
#include "lsh.h"



void gen_ExAudiences_vlsh(std::priority_queue<canducate_user> &top_k,
                          std::unordered_map<int,std::vector<int>> user_maps_pool,
                          std::unordered_map<int,std::vector<int>> user_maps_seed,
                          std::vector<bucket_info> centroid_angle,
                          int n_bit,
                          int n_feats,
                          int k,
                          float * data,
                          float *queries);
void calculate_centroid_angle(std::vector<bucket_info> &centroid_angle, std::unordered_map<int, std::vector<int>> user_maps_seed, float *queries, int n_feats, int n_bit);
#endif