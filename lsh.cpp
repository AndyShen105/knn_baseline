//
// Created by Hang Shen on 2018/7/12.
//


#include "lsh.h"


using namespace std;
void save_topk( priority_queue<canducate_user> top_k, string file_path){
    std::ofstream outFile;
    outFile.open(file_path);
    while(!top_k.empty()){
        int temp = top_k.top().sn;
        outFile << temp <<std::endl;
    }
    outFile.close();
}
void load_baseline(std::vector<int> &baseline, string file_path, char delimiter){
    string line;
    std::ifstream in(file_path);
    int i=0;
    while(getline(in, line)) {
        std::stringstream lineStream(line);
        std::string item;
        while (getline(lineStream, item, delimiter)) {
            baseline.push_back(stoi(item));
        }
    }

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
    int ss_index = seed_index*n_feats;
    int sp_index = pool_index*n_feats;
    for(int i=0; i<n_feats; i++){
        int s_index = ss_index+i;
        int p_index = sp_index+i;
        dis += data[p_index]*queries[s_index];
        datanorm += data[p_index]*data[p_index];
        querynorm += queries[s_index]*queries[s_index];
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
    clock_t time1 = clock();
    long long sum_calculations_cosin = 0;
    for(int i=0; i<n_cycle; i++){
        vector<int> &seed = user_maps_seed[i];
        vector<int> &pool = user_maps_pool[i];
        int temp = seed.size()*pool.size();
        if(temp>0)
            sum_calculations_cosin += seed.size()*pool.size();
        for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
            temp_user = calculate_similarity(seed, *pool_index, n_feats, data, queries);

            if (top_k.size() == k && temp_user.sim > top_k.top().sim )
                top_k.pop();
            if (top_k.size() < k )
                top_k.push(temp_user);
        }

    }
    cout<<"lsh sum calculate cosine : "<<sum_calculations_cosin<<endl;
}

void pre_user_pool(vector<int> &user_bucket_info,
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
        user_bucket_info.push_back(bucket_no);
    }

}

void gen_ExAudiences_lsh_based(priority_queue<canducate_user> &top_k,
                               unordered_map<int,vector<int>> &user_seed_pool,
                               vector<int> user_bucket_info,
                               int n_user,
                               int n_feats,
                               int k,
                               float * data,
                               float *queries) {

    canducate_user temp_user;
    temp_user.sim = -1000.000;
    int bucket_no;
    long long tempcount=0;
    for (int i = 0; i < n_user; i++) {
        if(i%10000==0)
            cout<<i<<endl;
        bucket_no = user_bucket_info[i];
        vector<int> &seed = user_seed_pool[bucket_no];
        temp_user = calculate_similarity(seed, i, n_feats, data, queries);

        if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
            top_k.pop();
        if (top_k.size() < k && temp_user.sim != -1000.0)
            top_k.push(temp_user);
    }
    cout<<"sum cosine calculation of vlsh: "<<tempcount<<endl;
}
#define CLOCKS_PER_SECOND 1000000.0
//
//int main(){
//    std::vector<int> result;
//    load_baseline(result, "data/sigMatrix.txt", ' ');
//    cout<<result[0]<<endl;
//    for(int i=0; i<result.size();i++){
//        cout<<result[i]<<endl;
//    }
//    cout<<"sdad"<<endl;
//    return 0;
//}
