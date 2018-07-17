//
// Created by Hang Shen on 2018/7/12.
//

#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <regex>
#include <queue>
#include <util.h>
#include "data.h"
#include "lsh.h"
using namespace std;

//set the random interval
std::uniform_real_distribution<double> dist(-1.0, 1.0);

int signature_bit(float *data, float **planes, int index, int n_feats, int n_plane){
    /*
    LSH signature generation using random projection
    Returns the signature bits.
    */
    int sig = 0;
    for(int i=0; i<n_plane; i++){
        sig <<= 1;
        float dot_product = 0.0;
        for(int j=0; j<n_feats;j++){
            dot_product += data[j+index]*planes[i][j];
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
        pMatrix[i] = new float[n_feats];
        for(int j=0; j<n_feats;j++){
            pMatrix[i][j] = dist(rng);
        }
    }
    return pMatrix;
}

void save_hashFunc(float **sigMatrix, int n_feats, int n_plane){
    std::ofstream outFile;
    outFile.open("data/sigMatrix.txt");
    for(int i = 0; i < n_plane; i++){
        for(int j=0; j<n_feats-1;j++){
            outFile << sigMatrix[i][j] <<" ";
        }
        outFile << sigMatrix[i][n_feats-1] << std::endl;
    }
    outFile.close();
}

void load_hashFunc(float **sigMatrix, int n_feats, int n_plane, char delimiter){
    string line;
    std::ifstream in("data/sigMatrix.txt"); 
    int i=0;
    while(getline(in, line)) { 
        std::stringstream lineStream(line);
        std::string item;
        int j=0;
        while (getline(lineStream, item, delimiter)) {
            sigMatrix[i][j]=stof(item);
            j++;
        }
        i++;
    }

}

void user_map(unordered_map<int, vector<int>> &user_maps, float *data,float **hash_func, int n_users, int n_feats, int n_plane){
    int index = 0;
    //unordered_map<int, vector<int>> user_maps;
    for (int i = 0; i < n_users; i++) { // For each user
        index = i*n_feats;
        int bucket_no = signature_bit(data, hash_func, index, n_feats, n_plane);
        user_maps[bucket_no].push_back(i);
    }
    
}

float get_cosine_dis(int seed_index,
                    int pool_index,
                    int n_feats,
                    float * data,
                    float *queries){
    float dis = 0.0f;
    float datanorm = 0.0f;
    float querynorm = 0.0f;
    for(int i=0; i<n_feats; i++){
        dis += data[seed_index*n_feats+i]*queries[pool_index*n_feats+i];
        datanorm = data[seed_index*n_feats+i]*data[seed_index*n_feats+i];
        querynorm = queries[pool_index*n_feats+i]*queries[pool_index*n_feats+i];
    }   
    return dis/(sqrt(datanorm) * sqrt(querynorm));
}

void gen_ExAudiences(priority_queue<canducate_user> &top_k,
                    unordered_map<int, 
                    vector<int>> user_maps_pool, 
                    unordered_map<int, 
                    vector<int>> user_maps_seed, 
                    int n_bit,
                    int n_feats,
                    int k,
                    float * data,
                    float *queries){
    int n_cycle = pow(2, n_bit);
    for(int i=0; i<n_cycle; i++){
        vector<int> seed = user_maps_seed[i];
        vector<int> pool = user_maps_pool[i];
        for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
            float max_sim = 0.0;
            for(vector<int>::const_iterator seed_index=seed.cbegin(); seed_index!=seed.cend(); seed_index++){
                float max_sim = max(get_cosine_dis(*seed_index, *pool_index, n_feats, data, queries), max_sim);
            }
            if(isnan(max_sim)){
            continue;
            }
            canducate_user temp_user;
            temp_user.sn = i;
            temp_user.sim = max_sim;
            if (top_k.size()<=k){
                top_k.push(temp_user);
            }
            else{
                if(top_k.top().sim<temp_user.sim){
                    top_k.pop();
                    top_k.push(temp_user);
                }
            }
        }
    }
        
}

#define CLOCKS_PER_SECOND 1000000.0
int main(){
    float * data,*queries;
    queries = (float *)malloc((int64_t)sizeof(float)*50*247753);
    csv_to_array(&queries, "/home/andy_shen/data/MovieLens/q.txt", 247753, 50);
    data = (float *)malloc((int64_t)sizeof(float)*50*33670);
    csv_to_array(&queries, "/home/andy_shen/data/MovieLens/p.txt", 33670, 50);
    cout<<"read data done"<<endl;
    clock_t time1 = clock();
    float** sig_maritx = gen_signature_matrix(50, 5);
    clock_t time2 = clock();
    unordered_map<int, vector<int>> user_maps_pool;
    user_map(user_maps_pool, queries, sig_maritx, 247753, 50, 5);
    unordered_map<int, vector<int>> user_maps_seed;
    user_map(user_maps_seed, queries, sig_maritx, 33670, 50, 5);
    // for(unordered_map<int, vector<int>>::iterator iter=user_maps.begin();iter!=user_maps.end();iter++)
    //     cout<<"key value is"<<iter->first<<" the mapped value is "<< iter->second.size()<<endl;
    clock_t time3 = clock();
    cout<<"init matrix time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;
    cout<<"init matrix time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;
    return 0;
}