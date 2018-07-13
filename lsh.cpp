//
// Created by Hang Shen on 2018/7/12.
//

#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <fstream>
#include <data.h>
#include "lsh.h"

//set the random interval
std::uniform_real_distribution<double> dist(-1.0, 1.0);

int signature_bit(float *data, float **planes, int n_feats, int n_plane){
    /*
    LSH signature generation using random projection
    Returns the signature bits.
    */
    int sig = 0;
    for(int i=0; i<n_plane; i++){
        sig <<= 1;
        float dot_product = 0.0;
        for(int j=0; j<n_feats;j++){
            dot_product += data[j]*planes[i][j];
        }
        if(dot_product >= 0){
            sig |= 1;
        }
    }
    return sig;
}

float** gen_signature_matrix(int n_feats, int n_plane){

    //Initialize with non-deterministic seeds
    std::mt19937 rng;  
    float **pMatrix = new float*[n_plane];
    rng.seed(std::random_device{}()); 
    for(int i=0; i<n_plane; i++){
        pMatrix[i]=new float[n_feats];
        for(int j=0; j<n_feats;j++){
            pMatrix[i][j]= dist(rng);
        }
    }
    return pMatrix;
}

void save_hashFunc(float **sigMatrix, int n_feats, int n_plane){
    std::ofstream outFile;
    outFile.open("data/sigMatrix.txt");
    for(int i = 0; i < n_plane; i++){
        for(int j=0; j<n_feats-1;j++){
            outFile << sigMatrix[i][j]<<" ";
        }
        outFile << sigMatrix[i][n_feats-1]<<std::endl;
    }
    outFile.close();
}

void load_hashFunc(float **sigMatrix, int n_feats, int n_plane){
    delim_to_array(sigMatrix, "data/sigMatrix.txt", n_feats, n_plane, ' ');

}

int main(){
    float **p=gen_signature_matrix(4,4);
    
    for(int i=0; i<4; i++){
       
        for(int j=0; j<4;j++){
            std::cout<<p[i][j]<<" ";
        }
         std::cout<<std::endl;
    }
    save_hashFunc(p, 4, 4);
    float u[4]={0.1, -0.4, 0.5, 0.2};
    int sig=signature_bit(u, p, 4,4);
    std::cout<<sig;
    return 0;
}