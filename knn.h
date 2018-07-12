//
// Created by Hang Shen on 2018/3/15.
//
#ifndef KNN_H
#define KNN_H

#include "util.h"
#include <queue>
using namespace std;
void serial_distances(float *data[], priority_queue <canducate_user> &distances, int n, int r, float *query[], int q, int k);
void knn_distance(float **data, priority_queue <canducate_user> &distances, int n, int r, float **query, int i, int k);

#endif