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

#endif