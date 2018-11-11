//
// Created by andyshen on 18-11-10.
//

#ifndef KNN_HASHFUNC_H
#define KNN_HASHFUNC_H
#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <iomanip>
#include <regex>
#define LONG_MAX      2147483647L
#include <unordered_map>
using namespace std;
int signature_bit(float *data, float **planes, int index, int n_feats, int n_plane);
float** gen_signature_matrix(int n_feats, int n_plane);
void user_map(unordered_map<int, vector<int>> &user_maps,
              float *data,
              float **hash_func,
              int n_users,
              int n_feats,
              int n_plane);
#endif //KNN_HASHFUNC_H
