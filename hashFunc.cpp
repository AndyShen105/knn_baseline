//
// Created by andyshen on 18-11-10.
//

#include "hashFunc.h"


//set the random interval
//uniform_real_distribution<double>
//std::normal_distribution  <double> dist(-1.0, 1.0);
std::uniform_real_distribution <double> dist(-1.0, 1.0);

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

float** gen_best_local_signature_matrix(int n_feats, int n_plane, int n_sig,int n_query, float *data, float *queries){

    vector<float**> all_sig;
    int index = 0;
    long min_value = 2147483647L;


    for(int i=0; i<n_sig; i++){
        unordered_map<int, vector<int>> user_maps_seed;
        float** sig_maritx = gen_signature_matrix(n_feats, n_plane);
        user_map(user_maps_seed, queries, sig_maritx, n_query, n_feats, n_plane);
        long sum_calcu=0;
        for (int i=0; i<pow(2, n_plane); i++){
            sum_calcu += pow(user_maps_seed[i].size(), 2);
        }
        //cout<<sum_calcu<<endl;
        if(sum_calcu < min_value){
            min_value = sum_calcu;
            index = i;
        }
        all_sig.push_back(sig_maritx);
    }

    for(int i=0; i<n_sig; i++){
        if (i==index)
            continue;
        free(all_sig[i]); all_sig[i] = nullptr;
    }

    return all_sig[index];

}


//
////
//int main(){
//#define CLOCKS_PER_SECOND 1000000.0
//#include "data.h"
//    float *queries, *data;
//    data = (float *)malloc((int64_t)sizeof(float)*50*552339);
//    csv_to_array(&data, "/home/andyshen/data/Yelp/q.txt", 552339, 50);
//    queries = (float *)malloc((int64_t)sizeof(float)*50*77079);
//    csv_to_array(&queries, "/home/andyshen/data/Yelp/p.txt", 77079, 50);
////    data = (float *)malloc((int64_t)sizeof(float)*50*247753);
////    csv_to_array(&data, "/home/andyshen/data/MovieLens/q.txt", 247753, 50);
////    queries = (float *)malloc((int64_t)sizeof(float)*50*33670);
////    csv_to_array(&queries, "/home/andyshen/data/MovieLens/p.txt", 33670, 50);
//    cout<<"read data done"<<endl;
//    clock_t time1 = clock();
//    int arr[15]={1,5,10,20,30,40,50,60,70,80,90,100,120,140,160};
////    for (int n=0;n<015;n++){
////        for(int i=0; i<10; i++) {
////            float **temp = gen_best_local_signature_matrix(50, 5, arr[n], 77079, data, queries);
////            clock_t time2 = clock();
////            //cout << "time: " << (time2 - time1) / CLOCKS_PER_SECOND << endl;
////            unordered_map<int, vector<int>> user_maps_seed;
////            unordered_map<int, vector<int>> user_maps_pool;
////            user_map(user_maps_pool, data, temp, 552339, 50, 5);
////            user_map(user_maps_seed, queries, temp, 77079, 50, 5);
////            long sum = 0;
////            for (int i = 0; i < 32; i++) {
////                sum += user_maps_seed[i].size() * user_maps_pool[i].size();
////                // cout<<user_maps_seed[i].size()<<" "<<user_maps_pool[i].size()<<endl;
////            }
////            cout  << sum<< "," ;
////        }
////        cout<<endl;
////    }
//
//    long all_sum=0;
//    long minv=LONG_MAX;
//    long maxv=0;
//    for(int i=0; i<50; i++){
//        float** sig_maritx = gen_signature_matrix(50,5);
//        unordered_map<int, vector<int>> user_maps_seeda;
//        unordered_map<int, vector<int>> user_maps_poola;
//        user_map(user_maps_poola, data, sig_maritx,552339, 50, 5);
//        user_map(user_maps_seeda, queries, sig_maritx, 77079, 50, 5);
//        long temsum=0;
//        for(int i=0; i<32; i++)
//        {temsum+=user_maps_seeda[i].size()*user_maps_poola[i].size(); }
//        cout<<"random "<<i<<":"<<temsum<<endl;
//        minv=min(minv,temsum);
//        maxv=max(maxv,temsum);
//        all_sum+=temsum;
//
//        free(sig_maritx);sig_maritx= nullptr;
//    }
//    cout<<"avr: "<<all_sum/50<<endl;
//    cout<<"min: "<<minv<<endl;
//    cout<<"max: "<<maxv<<endl;
//    return 0;
//}