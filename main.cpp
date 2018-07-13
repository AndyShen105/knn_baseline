#include <iostream>
#include<vector>
#include<queue>
#include <time.h>
#include "data.h"
#include "main.h"
#include "knn.h"
#include "sort.h"
#include "util.h"



int main(int argc, const char * argv[]) {

    /* Parameters */
    int n, r, q, k, s;
    float *data, *queries;
    n = atoi( argv[1] ); // Size of dataset (rows)
    r = atoi( argv[2] ); // Number of dimensions (columns)
    k = atoi( argv[3] ); // Nearest neighbors
    q = atoi( argv[4] ); // Number of queries

    string inputFile = argv[5];
    string queryFile = argv[6];
    cout<<"inputFile: "<< inputFile <<endl;
    cout<<"queryFile: "<< queryFile <<endl;

    data = (float *)malloc((int64_t)sizeof(float)*n*r); check_alloc(data);
    
    queries = (float *)malloc((int64_t)sizeof(float)*q*r); check_alloc(queries);
    float* thresholds = (float *)malloc(sizeof(int)*q); check_alloc(thresholds);

    clock_t timeStart = clock();
    
    csv_to_array(&queries, queryFile, q, r);
    csv_to_array(&data, inputFile, n, r);
 
    clock_t timeEnd = clock();
    double duration = (timeEnd-timeStart)/CLOCKS_PER_SECOND ;
    cout<<"Duration of reading data: "<< duration <<endl;
    cout<<"finsh reading data"<<endl;

    std::cout << "Data size:  " << (sizeof(float)*n*r)/(1024*1024) << "MB (" << n << " x " << r << ")" << std::endl;
    std::cout << "Query size: " << (sizeof(float)*q*r)/(1024*1024) << "MB (" << q << " x " << r << ")" << std::endl;
    std::cout << "Dist. size: " << (sizeof(float)*n*q)/(1024*1024) << "MB (" << n << " x " << q << ")" << std::endl;
    std::cout << "Index size: " << (sizeof(float)*k*q)/(1024*1024) << "MB (" << k << " x " << q << ")" << std::endl;
    std::cout << "Total:      " << ((sizeof(float)*n*r) + (sizeof(float)*q*r) + (sizeof(float)*n*q) + (sizeof(float)*k*q))/(1024*1024) << "MB" << std::endl;

    priority_queue<canducate_user> top_k;
    clock_t timeStart2 = clock();
    /*
    data: audience pool
    query: audience seeds
    distance: similarty
    n: size of audience pool
    r: dimensions of audience
    q: numbers of queries
    */
    knn_distance(&data, top_k, n, r, &queries, q, k);
    clock_t timeEnd2 = clock();
    double processTime = (timeEnd2-timeStart2)/CLOCKS_PER_SECOND ;
    cout<<"Duration of process time: "<< processTime <<endl;
    while(!top_k.empty()){
        cout<<"NO.:"<<top_k.top().sn<<" Sim:"<<top_k.top().sim<<endl;
        top_k.pop();
    }
    free(data); data = nullptr;
    free(queries); queries = nullptr;

    return EXIT_SUCCESS;
}