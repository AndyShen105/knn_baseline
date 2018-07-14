//
// Created by Hang Shen on 2018/7/11.
//
#ifndef LSH_H
#define LSH_H

int signature_bit(float *data, float **planes, int index, int n_feats, int n_plane);
float** gen_signature_matrix(int n_feats, int n_plane);
void save_hashFunc(float *sigMatrix, int n_feats, int n_plane);
void load_hashFunc(float *sigMatrix);

#endif