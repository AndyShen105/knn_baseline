//
// Created by Hang Shen on 2018/3/15.
//
#ifndef KNN_H
#define KNN_H

void serial_distances(float *data[], float *distances[], int n, int r, int k, float *query[], int q);
void knn_distance(float **data, float **distances, int n, int r, int k, float **query, int i);

#endif