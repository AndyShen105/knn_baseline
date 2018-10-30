//
// Created by Hang Shen on 2018/7/12.
//

#include <iostream>
#include <iomanip>
#include <math.h>
#include <random>
#include <fstream>
#include <regex>
#include "data.h"
#include "lsh.h"
using namespace std;

//set the random interval
std::normal_distribution  <double> dist(-1.0, 1.0);

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

void user_map(unordered_map<int, vector<int>> &user_maps, 
            float *data,
            float **hash_func, 
            int n_users, 
            int n_feats, 
            int n_plane){
    int index = 0;
    //unordered_map<int, vector<int>> user_maps;
    for (int i = 0; i < n_users; i++) { // For each user
        index = i*n_feats;
        int bucket_no = signature_bit(data, hash_func, index, n_feats, n_plane);
        user_maps[bucket_no].push_back(i);
    }
    
}

canducate_user calculate_similarity(vector<int> seed, int pool_index, int n_feats, float * data, float *queries){
    float max_sim = 0.0f;
    for(vector<int>::const_iterator seed_index=seed.cbegin(); seed_index!=seed.cend(); seed_index++){
        max_sim = max(get_cosine_dis(*seed_index, pool_index, n_feats, data, queries), max_sim);
    }
    canducate_user temp_user;
    temp_user.sn = pool_index;
    if(isnan(max_sim)){
        temp_user.sim = 0.0;

    }
    else{
        temp_user.sim = max_sim;
    }
    return temp_user;
}

float get_cosine_dis(int seed_index,
                     int pool_index,
                     int n_feats,
                     float *data,
                     float *queries){
    float dis = 0.0f;
    float datanorm = 0.0f;
    float querynorm = 0.0f;
    for(int i=0; i<n_feats; i++){
        dis += data[seed_index*n_feats+i]*queries[pool_index*n_feats+i];
        datanorm += data[seed_index*n_feats+i]*data[seed_index*n_feats+i];
        querynorm += queries[pool_index*n_feats+i]*queries[pool_index*n_feats+i];
    }
    if (datanorm == 0.0 || querynorm == 0.0)
        return -1000.0;
    return dis/(sqrt(datanorm) * sqrt(querynorm));
}

void gen_ExAudiences(priority_queue<canducate_user> &top_k,
                    unordered_map<int,vector<int>> user_maps_pool, 
                    unordered_map<int,vector<int>> user_maps_seed, 
                    int n_bit,
                    int n_feats,
                    int k,
                    float * data,
                    float *queries){
    int n_cycle = pow(2, n_bit);
    canducate_user temp_user;
    for(int i=0; i<n_cycle; i++){
        vector<int> seed = user_maps_seed[i];
        vector<int> pool = user_maps_pool[i];
        for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){

            temp_user = calculate_similarity(seed, *pool_index, n_feats, data, queries);
            if (top_k.size() == k && temp_user.sim > top_k.top().sim )
                top_k.pop();
            if (top_k.size() < k  )
                top_k.push(temp_user);
        }
        //cout<<"Bucket"<<i<<" finished "<<endl;
    }
        
}

#define CLOCKS_PER_SECOND 1000000.0
#if DEBUG
int main(){
    float *data,*queries;
    data = (float *)malloc((int64_t)sizeof(float)*50*247753);
    csv_to_array(&data, "/home/andyshen/data/MovieLens/q.txt", 247753, 50);
    queries = (float *)malloc((int64_t)sizeof(float)*50*33670);
    csv_to_array(&queries, "/home/andyshen/data/MovieLens/p.txt", 33670, 50);
    cout<<"read data done"<<endl;
    clock_t time1 = clock();
    float** sig_maritx = gen_signature_matrix(50, 5);
    clock_t time2 = clock();
    unordered_map<int, vector<int>> user_maps_pool;
    user_map(user_maps_pool, queries, sig_maritx, 247753, 50, 5);
    unordered_map<int, vector<int>> user_maps_seed;
    user_map(user_maps_seed, queries, sig_maritx, 33670, 50, 5);
    clock_t time3 = clock();
    cout<<"init sigmatrix time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;
    cout<<"init lshmatrix time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;
    priority_queue<canducate_user> top_k;
    gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, 5, 50, 1000, data, queries);
    clock_t time4 = clock();
    while(!top_k.empty()){
        cout<<"NO.:"<<top_k.top().sn<<" Sim:"<<top_k.top().sim<<endl;
        top_k.pop();
    }
    cout<<"query time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;
    return 0;
}
#endif