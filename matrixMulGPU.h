//
// Created by Hang Shen on 2018/3/17.
//

#ifndef BASELINE_MATRIXMULCPU_H
#define BASELINE_MATRIXMULCPU_H

#include <cuda_runtime.h>
int matrixMultiply(float **data, float **querys, float **distance, int block_size, dim3 &dimsA, dim3 &dimsB);
#endif //BASELINE_MATRIXMULCPU_H
