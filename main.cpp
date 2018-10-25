#include <iostream>
#include<vector>
#include<queue>
#include <time.h>
#include "data.h"
#include "main.h"
#include "knn.h"
#include "sort.h"
#include "util.h"
#include "lsh.h"
#include "v_lsh.h"


int main(int argc, const char * argv[]) {

    /* Parameters */
    int n, r, q, k, s, flag;
    float *data, *queries;
    n = atoi( argv[1] ); // Size of dataset (rows)
    r = atoi( argv[2] ); // Number of dimensions (columns)
    k = atoi( argv[3] ); // Nearest neighbors
    q = atoi( argv[4] ); // Number of queries
    flag = atoi( argv[5] );
    string inputFile = argv[6];
    string queryFile = argv[7];

    cout<<"Size of dataset: "<<n<<endl;
    cout<<"Number of dimensions : "<<r<<endl;
    cout<<"Top-k: "<<k<<endl;
    cout<<"Number of Seeds: "<<q<<endl;
    cout<<"State: "<<flag<<endl;
    cout<<"InputFile: "<< inputFile <<endl;
    cout<<"QueryFile: "<< queryFile <<endl;

    data = (float *)malloc((int64_t)sizeof(float)*n*r); check_alloc(data);
    queries = (float *)malloc((int64_t)sizeof(float)*q*r); check_alloc(queries);

    clock_t timeStart = clock();
    
    csv_to_array(&queries, queryFile, q, r);
    csv_to_array(&data, inputFile, n, r);
 
    clock_t timeEnd = clock();
    double duration = (timeEnd-timeStart)/CLOCKS_PER_SECOND ;
    cout<<"Duration of reading data: "<< duration <<endl;
    cout<<"finsh reading data"<<endl;

    std::cout << "Data size:  " << (sizeof(float)*n*r)/(1024*1024) << "MB (" << n << " x " << r << ")" << std::endl;
    std::cout << "Query size: " << (sizeof(float)*q*r)/(1024*1024) << "MB (" << q << " x " << r << ")" << std::endl;
//    std::cout << "Dist. size: " << (sizeof(float)*n*q)/(1024*1024) << "MB (" << n << " x " << q << ")" << std::endl;
//    std::cout << "Index size: " << (sizeof(float)*k*q)/(1024*1024) << "MB (" << k << " x " << q << ")" << std::endl;
//    std::cout << "Total:      " << ((sizeof(float)*n*r) + (sizeof(float)*q*r) + (sizeof(float)*n*q) + (sizeof(float)*k*q))/(1024*1024) << "MB" << std::endl;



    float** sig_maritx = gen_signature_matrix(50, 5);
    priority_queue<canducate_user> top_k;
    switch  (flag){
        case 0:
        {
            clock_t timeStart2 = clock();
            knn_distance(&data, top_k, n, r, &queries, q, k);
            clock_t timeEnd2 = clock();
            double processTime = (timeEnd2-timeStart2)/CLOCKS_PER_SECOND ;
            cout<<"Duration of process time: "<< processTime <<endl;
            break;
        }
        case 1:
        {
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, queries, sig_maritx, n, r, 5);
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx, q, r, 5);
            gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, 5, 50, k, data, queries);

            break;
        }
        case 2:
        {
            clock_t time1 = clock();
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, queries, sig_maritx, n, r, 5);
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx, q, r, 5);
            clock_t time2 = clock();
            cout<<"pre-process time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;

            gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, 5, 50, k, data, queries);
            clock_t time3 = clock();
            cout<<"lsh query time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;
            vector<int> lsh_re;
            while(!top_k.empty()){
                lsh_re.push_back(top_k.top().sn);
                top_k.pop();
            }
            //v-lsh
            vector<int > v_lsh_re;
            vector<bucket_info> centroid_angle;
            calculate_centroid_angle(centroid_angle, user_maps_seed, queries, r, 5);
            gen_ExAudiences_vlsh(top_k, user_maps_pool, user_maps_seed, centroid_angle,  5, r, k, data, queries);
            clock_t time4 = clock();
            cout<<"v_lsh query time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;
            while(!top_k.empty()){
                v_lsh_re.push_back(top_k.top().sn);
                top_k.pop();
            }
            cout<<"accuary: "<<calutate_acc(lsh_re, v_lsh_re)<<endl;

            break;
        }

    }
    while(!top_k.empty()){
        cout<<"NO.:"<<top_k.top().sn<<" Sim:"<<top_k.top().sim<<endl;
        top_k.pop();
    }
    free(data); data = nullptr;
    free(queries); queries = nullptr;

    return EXIT_SUCCESS;
}