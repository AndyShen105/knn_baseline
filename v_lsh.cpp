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
#include "util.h"
#include "data.h"
#include "lsh.h"
#include<cmath>

using namespace std;

#define DEBUG 0
//set the random interval
std::uniform_real_distribution<double> dist(-1.0, 1.0);
//------------------
int all_count = 0;
int save_calu_times = 0;
//------------------
struct bucket_info{
    int sn; // series number of bucket
    float* centroid; // centroid vector of bucket
    float theta_b;
    float centroid_sqrt;
};

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
    return dis/(sqrt(datanorm) * sqrt(querynorm));
}

//calculate cetroid and angle of seed pool
void calculate_centroid_angle(vector<bucket_info> &centroid_angle, unordered_map<int,vector<int>> user_maps_seed, float *queries, int n_feats, int n_bit){
    int n_cycle = pow(2, n_bit);
    
    for(int i=0; i<n_cycle; i++){
        vector<int> seed = user_maps_seed[i];
        // if size of bucket < 50, skip
        if(seed.size()<50)
            continue;
        float *centroid = new float[n_feats]();
        for(vector<int>::const_iterator seed_index=seed.cbegin(); seed_index!=seed.cend(); seed_index++){
            for(int j=0; j<n_feats; j++){
                centroid[j] += queries[(*seed_index)*n_feats+j];
            }
        }
        // mean of centroid
        float centroid_sqrt=0.0;
        for(int j=0; j<n_feats; j++){
            centroid[j] = centroid[j]/seed.size();
//            #if DEBUG
//            cout<<centroid[j]<<" ";
//            #endif
            centroid_sqrt += pow(centroid[j], 2);
        }

        //get angle B
        float theta_b = 0.0;
        for(vector<int>::const_iterator seed_index=seed.cbegin(); seed_index!=seed.cend(); seed_index++){
            float user_sqrt = 0.0;
            float dot = 0.0;
            for(int j=0; j<n_feats; j++){
                user_sqrt += pow(queries[(*seed_index)*n_feats+j], 2);
                dot += queries[(*seed_index)*n_feats+j] * centroid[j];
            }
            float cos = dot/(sqrt(user_sqrt)*sqrt(centroid_sqrt));
            theta_b = max(acos(cos), theta_b);
        }
        bucket_info temp;
        temp.sn = i;
        temp.centroid = centroid;
        temp.theta_b = theta_b;
        temp.centroid_sqrt = centroid_sqrt;
        centroid_angle.push_back(temp);
        cout<<centroid_angle.size()<<endl;
        delete centroid;
    }
    #if DEBUG
    //print each theta b of bucket
    for(vector<bucket_info>::const_iterator seed_index=centroid_angle.cbegin(); seed_index!=centroid_angle.cend(); seed_index++){
        cout<<(*seed_index).sn<< "  "<<(*seed_index).theta_b<<endl;
    }
    #endif  
}

//calculate upperbound of each user in userpool.
void calculate_upperbound(vector<int> pool, 
                        float * data, 
                        float* centroid,
                        float centroid_sqrt,
                        float theta_b, 
                        float* upper_bound_list, 
                        int n_feats){
    int pool_sn = 0;
    for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
        float user_sqrt = 0.0;

        float dot = 0.0;
        for(int j=0; j<n_feats; j++){
            user_sqrt += pow(data[(*pool_index)*n_feats+j], 2);
            dot += data[(*pool_index)*n_feats+j] * centroid[j];
        }
        //cout<<"centroid_sqrt: "<<centroid_sqrt<<endl;
        float upper_bound=0.0;
        if(user_sqrt!=0) {
            float cos_theta = dot/(sqrt(user_sqrt)*sqrt(centroid_sqrt));
            float theta_ic = acos(cos_theta);
            if (theta_ic-theta_b>=0)
                upper_bound = cos(theta_ic-theta_b);
            else
                upper_bound = 1;
        }
        upper_bound_list[pool_sn++] = upper_bound;
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
        temp_user.sim = 0;

    }
    else{
        temp_user.sim = max_sim;
    }
    return temp_user;
}

void gen_ExAudiences(priority_queue<canducate_user> &top_k,
                    unordered_map<int,vector<int>> user_maps_pool, 
                    unordered_map<int,vector<int>> user_maps_seed, 
                    vector<bucket_info> centroid_angle,
                    int n_bit,
                    int n_feats,
                    int k,
                    float * data,
                    float *queries){
    int n_cycle = pow(2, n_bit);
    canducate_user temp_user;
    cout<< "k: "<<k<<endl;
    for(int i=0; i<n_cycle; i++){
        cout<<" ---------------------start---------------"<<endl;
        vector<int> seed = user_maps_seed[i];
        vector<int> pool = user_maps_pool[i];
        cout<<"pool size: "<<pool.size()<<endl;
        cout<<"seed size: "<<seed.size()<<endl;
        // if size<50, lsh, else v-lsh
        if (seed.size()<50){
            for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
                temp_user = calculate_similarity(seed, *pool_index, n_feats, data, queries);
                if (top_k.size()<k){
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
        else{
            float * centroid ;
            float theta_b ;
            int bucket_sn ;
            float centroid_sqrt;
            for(vector<bucket_info>::const_iterator index=centroid_angle.cbegin(); index!=centroid_angle.cend(); index++){
                if(i == (*index).sn){
                    centroid = (*index).centroid;
                    theta_b = (*index).theta_b;
                    centroid_sqrt = (*index).centroid_sqrt;
                    bucket_sn = (*index).sn;
                    break;
                }

            }
            float* upper_bound_list = new float(pool.size());
            cout<<"centroid_sqrt: "<<centroid_sqrt<<endl;
            calculate_upperbound(pool, data, centroid, centroid_sqrt, theta_b,  upper_bound_list,  n_feats);
            cout<<"caculate upper bound"<<endl;
//            for(int m=0; m<pool.size();m++)
//                cout<<upper_bound_list[m]<<endl;
            for(int j=0; j<pool.size(); j++){
                //cout<<"upper_bound_list: "<<upper_bound_list[j]<<endl;
                if(upper_bound_list[j]>top_k.top().sim){
                    temp_user.sn = pool[j];
                    temp_user.sim = upper_bound_list[j];
                    all_count++;
                    save_calu_times += seed.size();
                    
                }
                else{
                    temp_user = calculate_similarity(seed, pool[j], n_feats, data, queries);

                }

                if (j%1000==0)
                    cout<<"j: "<<j<<endl;
                if (top_k.size()>=k)
                    top_k.pop();
                top_k.push(temp_user);
                               
            }
            if(centroid==NULL)
                cout<<"invaild pointer"<<endl;
            cout<<"delete cebtroid"<<endl;
            //delete centroid;
            delete upper_bound_list;
        }
        cout<<"seed capacity: "<<seed.capacity()<<endl;
        cout<<"Bucket"<<i<<" finished "<<endl;
    }
        
}

#define CLOCKS_PER_SECOND 1000000.0

float calculate_angle(vector<float> a, vector<float> b){
    float a_n = 0.0;
    float b_n = 0.0;
    float dot = 0.0;
    for(vector<float>::const_iterator i=a.cbegin(); i!=a.cend(); i++){
         for(vector<float>::const_iterator j=b.cbegin(); j!=b.cend(); j++){
             dot+=(*i) * (*j);
             b_n+=(*j) * (*j);
        }
        a_n+=(*i) * (*i);
    }
    float cos = dot/(sqrt(a_n)*sqrt(b_n));
    float the = acos(0.0);
    std::cout<<cos<<endl;
    return the;
}
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
    for(int i=0; i<32; i++){
        cout<<user_maps_seed[i].size()<<endl;
    }
    clock_t time3 = clock();
    cout<<"init sigmatrix time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;
    cout<<"init lshmatrix time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;

    vector<bucket_info> centroid_angle;
    calculate_centroid_angle(centroid_angle, user_maps_seed, queries, 50, 5);


    priority_queue<canducate_user> top_k;


    gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, centroid_angle,  5, 50, 1000, data, queries);
    clock_t time4 = clock();
    while(!top_k.empty()){
        cout<<"NO.:"<<top_k.top().sn<<" Sim:"<<top_k.top().sim<<endl;
        top_k.pop();
    }
    cout<<"query time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;

    cout<<"all_count: "<<all_count<<endl;
    cout<<"save times: "<< save_calu_times<<endl;
    return 0;
    

}