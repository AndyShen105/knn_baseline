#include <iostream>
#include<vector>
#include <time.h>
#include "data.h"
#include "main.h"
#include "knn.h"
#include "sort.h"

using namespace std;


int main(int argc, const char * argv[]) {

    /* Parameters */
    int n, r, q, k, s;
    float *data, *queries;
    n = atoi( argv[1] ); // Size of dataset (rows
    r = atoi( argv[2] ); // Number of dimensions (columns)
    k = atoi( argv[3] ); // Nearest neighbors
    q = atoi( argv[4] ); // Number of queries

    string inputFile = argv[6];
    string queryFile = argv[7];

    data = (float *)malloc((int64_t)sizeof(float)*n*r); check_alloc(data);
    queries = (float *)malloc((int64_t)sizeof(float)*q*r); check_alloc(queries);
    float* thresholds = (float *)malloc(sizeof(int)*q); check_alloc(thresholds);

    clock_t timeStart = clock();
    csv_to_array(&data, inputFile, n, r);
    csv_to_array(&queries, queryFile, q, r);
    clock_t timeEnd = clock();
    double duration = (timeEnd-timeStart)/CLOCKS_PER_SECOND ;
    cout<<"Duration of reading data: "<< duration <<endl;
    cout<<"finsh reading data"<<endl;

    std::cout << "Data size:  " << (sizeof(float)*n*r)/(1024*1024) << "MB (" << n << " x " << r << ")" << std::endl;
    std::cout << "Query size: " << (sizeof(float)*q*r)/(1024*1024) << "MB (" << q << " x " << r << ")" << std::endl;
    std::cout << "Dist. size: " << (sizeof(float)*n*q)/(1024*1024) << "MB (" << n << " x " << q << ")" << std::endl;
    std::cout << "Index size: " << (sizeof(float)*k*q)/(1024*1024) << "MB (" << k << " x " << q << ")" << std::endl;
    std::cout << "Total:      " << ((sizeof(float)*n*r) + (sizeof(float)*q*r) + (sizeof(float)*n*q) + (sizeof(float)*k*q))/(1024*1024) << "MB" << std::endl;

    float* bdistances = (float *)malloc((int64_t)sizeof(float)*n*q); check_alloc(bdistances);
    int *knn = (int *)malloc(sizeof(int)*k*q); check_alloc(knn);

    clock_t timeStart2 = clock();
    knn_distance(&data, &bdistances, n, r, k, &queries, q);
    clock_t timeEnd2 = clock();
    double processTime = (timeEnd2-timeStart2)/CLOCKS_PER_SECOND ;
    cout<<"Duration of process time: "<< processTime <<endl;

    clock_t timeStart3 = clock();
    mksort(&bdistances, &thresholds, n, q, k);
    clock_t timeEnd3 = clock();
    double sortTime = (timeEnd3-timeStart3)/CLOCKS_PER_SECOND ;
    cout<<"Duration of sort time: "<< sortTime <<endl;

    free(bdistances); bdistances = nullptr;
    free(data); data = nullptr;
    free(queries); queries = nullptr;

    return EXIT_SUCCESS;
}