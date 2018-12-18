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
#include "v_lsh.h"
#include "hashFunc.h"
#include<cmath>

using namespace std;

#define DEBUG 0
//------------------

//------------------

void calculate_cetroid_s1(float *centroid, float *queries, float &centroid_sqrt, int n_feats, vector<int> seed){
    /*calucte centroid with center of cosine*/
    int centroid_index = 0;
    float max_cos = 0.0f;
    for(vector<int>::const_iterator seed_index=seed.cbegin(); seed_index!=seed.cend(); seed_index++){
        float min_sim = 10.0f;
        for(vector<int>::const_iterator index=seed.cbegin(); index!=seed.cend(); index++){
            min_sim = min(get_cosine_dis(*seed_index, *index, n_feats, queries, queries), min_sim);
        }
        if(max_cos<min_sim){
            centroid_index = *seed_index;
            max_cos = min_sim;
        }
    }
    for(int i=0; i<n_feats; i++){
        centroid[i] = queries[centroid_index*n_feats+i];
        centroid_sqrt += pow(centroid[i], 2);
    }
}

void calculate_cetroid_s2(float *centroid, float *queries, float &centroid_sqrt, int n_feats, vector<int> seed){
    /*calucte centroid with 均值向量*/
    for(int j=0; j<n_feats; j++){
        for(vector<int>::const_iterator seed_index=seed.cbegin(); seed_index!=seed.cend(); seed_index++){
            centroid[j] += queries[(*seed_index)*n_feats+j];
        }
        centroid[j] = centroid[j]/seed.size();
        centroid_sqrt += pow(centroid[j], 2);
    }
}


void calculate_centroid_angle(vector<bucket_info> &centroid_angle, unordered_map<int,vector<int>> user_maps_seed, float *queries, int n_feats, int n_bit){
    //calculate cetroid and angle of seed pool
    int n_cycle = pow(2, n_bit);
    for(int i=0; i<n_cycle; i++){
        vector<int> seed = user_maps_seed[i];
        float *centroid = new float[n_feats]();
        float centroid_sqrt=0.0;
        calculate_cetroid_s2(centroid, queries, centroid_sqrt, n_feats, seed);
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
        bucket_info temp =  bucket_info(seed.size());
        temp.sn = i;
        temp.centroid_copy(centroid, seed.size());
        temp.theta_b = theta_b;
        temp.centroid_sqrt = centroid_sqrt;
        centroid_angle.push_back(temp);

        delete []centroid;
    }
    #if DEBUG
    //print each theta b of bucket
    for(vector<bucket_info>::const_iterator seed_index=centroid_angle.cbegin(); seed_index!=centroid_angle.cend(); seed_index++){
        cout<<(*seed_index).sn<< "  "<<(*seed_index).theta_b<<endl;
    }
    #endif  
}


void calculate_upperbound(vector<int> pool, 
                        float * data, 
                        float* centroid,
                        float centroid_sqrt,
                        float theta_b, 
                        float* upper_bound_list, 
                        int n_feats){
    //calculate upperbound of each user in userpool.
    int pool_sn = 0;
    for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
        float user_sqrt = 0.0;

        float dot = 0.0;
        for(int j=0; j<n_feats; j++){
            user_sqrt += pow(data[(*pool_index)*n_feats+j], 2);
            dot += data[(*pool_index)*n_feats+j] * centroid[j];
        }

        float upper_bound=-1000.0;
        if(user_sqrt!=0.0) {
            float cos_theta = dot/(sqrt(user_sqrt)*sqrt(centroid_sqrt));
            float theta_ic = acos(cos_theta);

            if (theta_ic-theta_b>0.0) {
                upper_bound = cos(theta_ic - theta_b);
            }
            else{
                upper_bound = 1.0;
            }

        }

        upper_bound_list[pool_sn++] = upper_bound;
    }
}

float calculate_upperbound_per_user(int index,
                          float * data,
                          float* centroid,
                          float centroid_sqrt,
                          float theta_b,
                          int n_feats){
    //calculate upperbound of user .
    float user_sqrt = 0.0;
    float dot = 0.0;
    for(int j=0; j<n_feats; j++){
        user_sqrt += pow(data[index*n_feats+j], 2);
        dot += data[index*n_feats+j] * centroid[j];
    }
    float upper_bound=-1000.0;
    if(user_sqrt!=0.0) {
        float cos_theta = dot/(sqrt(user_sqrt)*sqrt(centroid_sqrt));
        float theta_ic = acos(cos_theta);

        if (theta_ic-theta_b>0.0) {
            upper_bound = cos(theta_ic - theta_b);
        }
        else{
            upper_bound = 1.0;
        }

    }
    return upper_bound;
}

void gen_ExAudiences_vlsh(priority_queue<canducate_user> &top_k,
                    unordered_map<int,vector<int>> user_maps_pool,
                    unordered_map<int,vector<int>> user_maps_seed,
                    priority_queue<uncertain_user> &user_pool,
                    vector<bucket_info> centroid_angle,
                    int n_bit,
                    int n_feats,
                    int k,
                    float * data,
                    float *queries){

    int n_cycle = pow(2, n_bit);
    long long all_count = 0;
    long long save_calu_times = 0;
    long long sum_save_times = 0;
    canducate_user temp_user;
    uncertain_user user;

    for(int i=0; i<n_cycle; i++){

        if (user_maps_seed[i].empty()){
            continue;
        }
        vector<int> &seed = user_maps_seed[i];
        vector<int> &pool = user_maps_pool[i];

        // if size<10, lsh, else v-lsh
        if (seed.size()<10){
            for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
                sum_save_times+=seed.size();
                temp_user = calculate_similarity(seed, *pool_index, n_feats, data, queries);
                if (top_k.size() == k && temp_user.sim > top_k.top().sim )
                    top_k.pop();
                if (top_k.size() < k )
                    top_k.push(temp_user);
            }
        }
        else{
            for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
                float upper_bound = calculate_upperbound_per_user(*pool_index, data, centroid_angle[i].centroid,
                                                                  centroid_angle[i].centroid_sqrt,
                                                                  centroid_angle[i].theta_b, n_feats);
                if (upper_bound == 1.0) {
                    temp_user = calculate_similarity(seed, i, n_feats, data, queries);
                    if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
                        top_k.pop();
                    if (top_k.size() < k && temp_user.sim != -1000.0)
                        top_k.push(temp_user);
                } else {
                    user.bucket_no = i;
                    user.index = *pool_index;
                    user.upperbound = upper_bound;
                    user_pool.push(user);
                }

            }
        }
    }
    cout<<"size of candicate user with uppbound"<<user_pool.size()<<endl;
    while(!user_pool.empty()){
        user = user_pool.top();
        user_pool.pop();
        int bucket_no = user.bucket_no;
        int index = user.index;
        float upperbound = user.upperbound;
        if( upperbound >= top_k.top().sim ){
            vector<int> &seed = user_maps_seed[bucket_no];
            temp_user = calculate_similarity(seed, index, n_feats, data, queries);
            if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
                top_k.pop();
            if (top_k.size() < k && temp_user.sim != -1000.0)
                top_k.push(temp_user);
        }else{
            all_count++;
            save_calu_times += user_maps_seed[bucket_no].size();
        }
    }
    cout<<"all_count: "<<all_count<<endl;
    cout<<"save times: "<< save_calu_times<<endl;

}
/*
void gen_ExAudiences_vlsh_first(priority_queue<canducate_user> &top_k,
                          int n_user,
                          unordered_map<int,vector<int>> &user_maps_pool,
                          priority_queue<uncertain_user> &user_pool,
                          const  vector<bucket_info> &centroid_angle,
                          int n_feats,
                          int k,
                          float **hash_func,
                          float * data,
                          float *queries) {

    canducate_user temp_user;
    temp_user.sim = -1000.0;
    uncertain_user user;
    int bucket_no;
    int tempcount=0;
    for (int i = 0; i < n_user; i++) {
        if(i%10000==0)
            cout<<i<<endl;
        bucket_no = signature_bit(data, hash_func, i * n_feats, n_feats, 5);
        vector<int> &seed = user_maps_pool[5];

        if (seed.size() < 10) {

            temp_user = calculate_similarity(seed, i, n_feats, data, queries);
        } else {

            float upper_bound = calculate_upperbound_per_user(i, data, centroid_angle[bucket_no].centroid,
                                                              centroid_angle[bucket_no].centroid_sqrt,
                                                              centroid_angle[bucket_no].theta_b, n_feats);
            if (upper_bound == 1.0) {
                temp_user = calculate_similarity(seed, i, n_feats, data, queries);
            } else {
                user.bucket_no = bucket_no;
                user.index = i;
                user.upperbound = upper_bound;
                user_pool.push(user);

            }

        }
        if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
            top_k.pop();
        if (top_k.size() < k && temp_user.sim != -1000.0)
            top_k.push(temp_user);
    }
    cout<<"sum cosine calculation of vlsh: "<<tempcount<<endl;
}

 *
 */


#if DEBUG
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

    user_map(user_maps_seed, queries, sig_maritx, 33670,50, 5);
    for(int i=0; i<32; i++){
        cout<<user_maps_seed[i].size()<<endl;
    }
    clock_t time3 = clock();
    cout<<"init sigmatrix time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;
    cout<<"init lshmatrix time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;

    vector<bucket_info> centroid_angle;
    calculate_centroid_angle(centroid_angle, user_maps_seed, queries, 50, 5);

#if DEBUG
    for(vector<bucket_info>::const_iterator index=centroid_angle.cbegin(); index!=centroid_angle.cend(); index++){
        cout<<"centroid_sqrt: "<<(*index).centroid_sqrt<<endl;
        for(int i=0; i<50 ;i++)
            cout<< (*index).centroid[i]<<" ";
        cout<<endl;
    }
#endif

    priority_queue<canducate_user> top_k;


    gen_ExAudiences_vlsh(top_k, user_maps_pool, user_maps_seed, centroid_angle,  5, 50, 1000, data, queries);
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
#endif
