//
// Created by Hang Shen on 2018/7/11.
//
#ifndef LSH_H
#define LSH_H

int signature_bit(float *data, float **planes, int index, int n_feats, int n_plane);
float** gen_signature_matrix(int n_feats, int n_plane);
void save_hashFunc(float *sigMatrix, int n_feats, int n_plane);
void load_hashFunc(float *sigMatrix);
void user_map(float *data,float **hash_func, int n_users, int n_feats, int n_plane);
// void gen_ExAudiences(priority_queue<canducate_user> &top_k,
//                     unordered_map<int,vector<int>> user_maps_pool, 
//                     unordered_map<int,vector<int>> user_maps_seed, 
//                     int n_bit,
//                     int n_feats,
//                     int k,
//                     float * data,
//                     float *queries);
#endif