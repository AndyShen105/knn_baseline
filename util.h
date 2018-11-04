//
// Created by Hang Shen on 2018/3/15.
//
#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <iomanip>
#include <queue>

typedef struct canducate_user{
        int sn;
        float sim;
        friend bool operator < (canducate_user s1,canducate_user s2){
            return s1.sim > s2.sim;
        }
    }canducate_user;

typedef struct uncertain_user{
    int bucket_no;
    int index;
    float upperbound;
    friend bool operator < (uncertain_user s1, uncertain_user s2){
        return s1.upperbound > s2.upperbound;
    }
}uncertain_user;

struct bucket_info{
    int sn; // series number of bucket
    float* centroid; // centroid vector of bucket
    float theta_b;
    float centroid_sqrt;
    bucket_info(int size){
        centroid = new float[size]();
    }
    void centroid_copy(float* centroid_copy, int size){
        for (int i=0; i<size; i++){
            centroid[i] = centroid_copy[i];
        }
    }
};

float calutate_acc(std::vector<int>lsh, std::vector<int> v_lsh);
#endif