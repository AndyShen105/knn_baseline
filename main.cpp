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
    int n, r, q, k;
    float *data, *queries;
    n = 10; // Size of dataset (rows
    r = 50; // Number of dimensions (columns)
    k = 5; // Nearest neighbors
    q = 5; // Number of queries

    string inputFile = "/Users/hangshen/Desktop/data_ming/recommdation/MFRetrieval-master/Netflix/q2.txt";
    string queryFile = "/Users/hangshen/Desktop/data_ming/recommdation/MFRetrieval-master/Netflix/p2.txt";

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

    clock_t timeStart2 = clock();
    mksort(&bdistances, &thresholds, n, q, k);
    clock_t timeEnd2 = clock();
    double processTime = (timeEnd2-timeStart2)/CLOCKS_PER_SECOND*1000.0;
    cout<<"Duration of process time: "<< processTime <<endl;

    clock_t timeStart3 = clock();
    mksort(&bdistances, &thresholds, n, q, k);
    clock_t timeEnd3 = clock();
    double sortTime = (timeEnd2-timeStart2)/CLOCKS_PER_SECOND*1000.0 ;
    cout<<"Duration of sort time: "<< sortTime <<endl;
    return 0;
}