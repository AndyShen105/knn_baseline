//creat by andyshen on 3.3.2019
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <helper_cuda.h>
#include "cudaopt.h"


#define CLOCKS_PER_SECOND 1000000.0
#define DEBUG 0
__global__ void matrixMulCosine(float* d_data, float * d_query, float * d_result, int *uIndex, int *qIndex, int q, int n, int n_feats);
__global__ void calcuSimCuda(float* d_data, float * d_query, float * d_result, int *uIndex, int *qIndex, int q, int n_feats);
__global__ void preProcess(float* d_data, float * d_query, float * d_Udot, float * d_Sdot, int n, int q, int n_feats);
__global__ void matrixMulCosineOpt(float* d_data, float * d_query, float * d_result, int *uIndex, int *qIndex, int q, int n,
                                   float* d_Udot, float* d_Sdot,int n_feats);

int* vector2int(vector<int> v){
    int* pV;
    int size = v.size();
    pV = (int *)malloc(sizeof(int)*size);
    for(int i=0; i<size; i++){
        pV[i] = v[i];
    }
    return pV;
}


void gen_ExAudiences_cudabase(priority_queue<canducate_user> &top_k,
                              unordered_map<int,vector<int>> user_maps_seed,
                              vector<int> user_bucket_info,
                              priority_queue<uncertain_user> &indexUser,
                              vector<bucket_info> centroid_angle,
                              int n_bit,
                              int n_feats,
                              int k,
                              int n,
                              int q,
                              float * data,
                              float *queries) {

    long long all_count = 0;
    long long save_calu_times = 0;
    long long sum_save_times = 0;
    long query_Bytes = sizeof(float) * q * n_feats;
    long data_Bytes = sizeof(float) * n * n_feats;

    // 申请数据device内存
    float *d_data, *d_query;
    cudaMalloc((void **) &d_query, query_Bytes);
    cudaMalloc((void **) &d_data, data_Bytes);
    // copy数据到device
    cudaMemcpy((void *) d_query, (void *) queries, query_Bytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void *) d_data, (void *) data, data_Bytes, cudaMemcpyHostToDevice);

    canducate_user temp_user;
    uncertain_user indexuser;

    for (int i = 0; i < n; i++)
    {
        int bucket_no = user_bucket_info[i];
        vector<int> &seed = user_maps_seed[bucket_no];
        if (seed.size()==0) {
            continue;
        }
        // if size<10, lsh, else v-lsh
        if (seed.size() < 1000) {
            sum_save_times += seed.size();
            temp_user = calculate_similarity(seed, i, n_feats, data, queries);
            if (top_k.size() == k && temp_user.sim > top_k.top().sim)
                top_k.pop();
            if (top_k.size() < k)
                top_k.push(temp_user);
        } else {
            int *uIndex, *qIndex;
            uIndex = (int *) malloc(sizeof(int));
            uIndex[0] = i;

            float upper_bound = calculate_upperbound_per_user(i, data, centroid_angle[bucket_no].centroid,
                                                                  centroid_angle[bucket_no].centroid_sqrt,
                                                                  centroid_angle[bucket_no].theta_b, n_feats);
            //如果上界不为1更新indexuser，跳过循环
            if (upper_bound != 1.0) {
                indexuser.bucket_no = bucket_no;
                indexuser.index = i;
                indexuser.upperbound = upper_bound;
                indexUser.push(indexuser);
            }else{
                //申请索引内存
                float *result, *d_result;
                int *d_qIndex, *d_uIndex;
                long result_Bytes = sizeof(float) * seed.size();
                long qIndex_Bytes = sizeof(int) * seed.size();
                long uIndex_Bytes = sizeof(int);
                result = (float *) malloc(result_Bytes);
                cudaMalloc((void **) &d_result, result_Bytes);
                cudaMalloc((void **) &d_qIndex, qIndex_Bytes);
                cudaMalloc((void **) &d_uIndex, uIndex_Bytes);
                //copy索引数据到device
                qIndex = vector2int(seed);
                cudaMemcpy((void *) d_qIndex, (void *) qIndex, qIndex_Bytes, cudaMemcpyHostToDevice);
                cudaMemcpy((void *) d_uIndex, (void *) uIndex, uIndex_Bytes, cudaMemcpyHostToDevice);

                dim3 blockSize(1024);
                dim3 gridSize((seed.size() + blockSize.x - 1) / blockSize.x);

                calcuSimCuda << < gridSize, blockSize >> >
                                            (d_data, d_query, d_result, d_uIndex, d_qIndex, seed.size(), n_feats);

                cudaMemcpy((void *) result, (void *) d_result, result_Bytes, cudaMemcpyDeviceToHost);
                cudaDeviceSynchronize();
                float tempSim = 0.0;
                for (int j = 0; j < seed.size(); j++) {
                    tempSim = max(tempSim, result[j]);
                }
                temp_user.sn = i;
                temp_user.sim = tempSim;
                if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
                    top_k.pop();
                if (top_k.size() < k && temp_user.sim != -1000.0)
                    top_k.push(temp_user);
                //释放申请的内存
                cudaFree(d_result);
                cudaFree(d_qIndex);
                cudaFree(d_uIndex);
                free(result);
                free(qIndex);
                free(uIndex);
            }

        }

    }
    cudaFree(d_data);
    cudaFree(d_query);
    cout << "size of candicate user with uppbound" << indexUser.size() << endl;
    while (!indexUser.empty()) { indexuser = indexUser.top();
        indexUser.pop();
        int bucket_no = indexuser.bucket_no;
        int index = indexuser.index;
        float upperbound = indexuser.upperbound;
        if (upperbound >= top_k.top().sim) {
            vector<int> &seed = user_maps_seed[bucket_no];
            temp_user = calculate_similarity(seed, index, n_feats, data, queries);
            if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
                top_k.pop();
            if (top_k.size() < k && temp_user.sim != -1000.0)
                top_k.push(temp_user);
        } else {
            all_count++;
            save_calu_times += user_maps_seed[bucket_no].size();
        }
    }
    cout << "all_count: " << all_count << endl;
    cout << "save times: " << save_calu_times << endl;
}

void gen_ExAudiences_cuda(priority_queue<canducate_user> &top_k,
                          unordered_map<int,vector<int>> user_maps_pool,
                          unordered_map<int,vector<int>> user_maps_seed,
                          priority_queue<uncertain_user> &indexUser,
                          vector<bucket_info> centroid_angle,
                          int n_bit,
                          int n_feats,
                          int k,
                          int n,
                          int q,
                          float * data,
                          float *queries){

    int n_cycle = pow(2, n_bit);
    long long all_count = 0;
    long long save_calu_times = 0;
    long long sum_save_times = 0;
    long query_Bytes = sizeof(float)*q*n_feats;
    long data_Bytes = sizeof(float)*n*n_feats;


    // 申请数据device内存
    float  *d_data, *d_query;
    cudaMalloc((void**)&d_query, query_Bytes);
    cudaMalloc((void**)&d_data, data_Bytes);
    // copy数据到device
    cudaMemcpy((void*)d_query, (void*) queries, query_Bytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)d_data, (void*) data, data_Bytes, cudaMemcpyHostToDevice);

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

            int *uIndex, *qIndex;
            uIndex = (int *)malloc(sizeof(int)*pool.size());
            int uCount = 0;
            for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
                float upper_bound = calculate_upperbound_per_user(*pool_index, data, centroid_angle[i].centroid,
                                                                  centroid_angle[i].centroid_sqrt,
                                                                  centroid_angle[i].theta_b, n_feats);
                if (upper_bound == 1.0) {
                    uIndex[uCount] = *pool_index;
                    uCount++;
                } else {
                    user.bucket_no = i;
                    user.index = *pool_index;
                    user.upperbound = upper_bound;
                    indexUser.push(user);
                }
            }
            //申请索引内存
            float *result, *d_result;
            int *d_qIndex, *d_uIndex;
            long result_Bytes = sizeof(float)*seed.size()*uCount;
            long qIndex_Bytes = sizeof(int)*seed.size();
            long uIndex_Bytes = sizeof(int)*pool.size();
            result = (float *)malloc(result_Bytes);
            cudaMalloc((void**)&d_result, result_Bytes);
            cudaMalloc((void**)&d_qIndex, qIndex_Bytes);
            cudaMalloc((void**)&d_uIndex, uIndex_Bytes);

            //copy索引数据到device
            qIndex = vector2int(seed);
            cudaMemcpy((void*)d_qIndex, (void*) qIndex, qIndex_Bytes, cudaMemcpyHostToDevice);
            cudaMemcpy((void*)d_uIndex, (void*) uIndex, uIndex_Bytes, cudaMemcpyHostToDevice);

            dim3 blockSize(1024);
            dim3 gridSize((uCount*seed.size() + blockSize.x - 1) / blockSize.x);

            matrixMulCosine<< < gridSize, blockSize >> >(d_data, d_query, d_result, d_uIndex, d_qIndex, seed.size(), uCount, n_feats);
            cudaDeviceSynchronize();
            cudaMemcpy((void*)result, (void*)d_result, result_Bytes, cudaMemcpyDeviceToHost);
            for(int i=0; i<uCount; i++){
                float tempSim = 0.0;
                for(int j=0; j<seed.size(); j++){
                    tempSim = max(tempSim, result[i*seed.size()+j]);
                }
                temp_user.sn = uIndex[i];
                temp_user.sim = tempSim;
                if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
                    top_k.pop();
                if (top_k.size() < k && temp_user.sim != -1000.0)
                    top_k.push(temp_user);

            }
            //释放申请的内存
            cudaFree(d_result);
            cudaFree(d_qIndex);
            cudaFree(d_uIndex);
            free(result);
            free(qIndex);
            free(uIndex);


        }
    }
    cudaFree(d_data);
    cudaFree(d_query);
    cout<<"size of candicate user with uppbound"<<indexUser.size()<<endl;
    while(!indexUser.empty()){
        user = indexUser.top();
        indexUser.pop();
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




void gen_ExAudiences_cudaOpt(priority_queue<canducate_user> &top_k,
                          unordered_map<int,vector<int>> user_maps_pool,
                          unordered_map<int,vector<int>> user_maps_seed,
                          priority_queue<uncertain_user> &indexUser,
                          vector<bucket_info> centroid_angle,
                          int n_bit,
                          int n_feats,
                          int k,
                          int n,
                          int q,
                          float * data,
                          float *queries){

    int n_cycle = pow(2, n_bit);
    long long all_count = 0;
    long long save_calu_times = 0;
    long long sum_save_times = 0;
    long query_Bytes = sizeof(float)*q*n_feats;
    long data_Bytes = sizeof(float)*n*n_feats;

    // 申请数据device内存
    float  *d_data, *d_query, *d_Udot, *d_Sdot;
    cudaMalloc((void**)&d_query, query_Bytes);
    cudaMalloc((void**)&d_data, data_Bytes);
    cudaMalloc((void**)&d_Udot, sizeof(float)*n);
    cudaMalloc((void**)&d_Sdot, sizeof(float)*q);
    // copy数据到device
    cudaMemcpy((void*)d_query, (void*) queries, query_Bytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)d_data, (void*) data, data_Bytes, cudaMemcpyHostToDevice);

    dim3 preProblockSize(1024);
    dim3 preProgridSize((n+q+preProblockSize.x - 1) / preProblockSize.x);

    preProcess<< < preProgridSize, preProblockSize >> >(d_data, d_query, d_Udot, d_Sdot, n, q, n_feats);
    cudaDeviceSynchronize();

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

            int *uIndex, *qIndex;
            uIndex = (int *)malloc(sizeof(int)*pool.size());
            int uCount = 0;
            for(vector<int>::const_iterator pool_index=pool.cbegin(); pool_index!=pool.cend(); pool_index++){
                float upper_bound = calculate_upperbound_per_user(*pool_index, data, centroid_angle[i].centroid,
                                                                  centroid_angle[i].centroid_sqrt,
                                                                  centroid_angle[i].theta_b, n_feats);
                if (upper_bound == 1.0) {
                    uIndex[uCount] = *pool_index;
                    uCount++;
                } else {
                    user.bucket_no = i;
                    user.index = *pool_index;
                    user.upperbound = upper_bound;
                    indexUser.push(user);
                }
            }
            //申请索引内存
            float *result, *d_result;
            int *d_qIndex, *d_uIndex;
            long result_Bytes = sizeof(float)*seed.size()*uCount;
            long qIndex_Bytes = sizeof(int)*seed.size();
            long uIndex_Bytes = sizeof(int)*pool.size();
            result = (float *)malloc(result_Bytes);
            cudaMalloc((void**)&d_result, result_Bytes);
            cudaMalloc((void**)&d_qIndex, qIndex_Bytes);
            cudaMalloc((void**)&d_uIndex, uIndex_Bytes);

            //printf("%s\n", cudaGetErrorString(cudaStatus));

            //copy索引数据到device
            qIndex = vector2int(seed);
            cudaError_t cudaStatus1=cudaMemcpy((void*)d_qIndex, (void*) qIndex, qIndex_Bytes, cudaMemcpyHostToDevice);
            cudaError_t cudaStatus2=cudaMemcpy((void*)d_uIndex, (void*) uIndex, uIndex_Bytes, cudaMemcpyHostToDevice);
            //cout<<cudaStatus1<<" "<<cudaStatus2<<endl;
            dim3 blockSize(1024);
            dim3 gridSize((uCount*seed.size() + blockSize.x - 1) / blockSize.x);
            matrixMulCosineOpt<< < gridSize, blockSize >> >(d_data, d_query,d_result, d_uIndex, d_qIndex, seed.size(), uCount, d_Udot, d_Sdot, n_feats);
            cudaDeviceSynchronize();

            cudaMemcpy((void*)result, (void*)d_result, result_Bytes, cudaMemcpyDeviceToHost);
            for(int i=0; i<uCount; i++){
                float tempSim = 0.0;
                for(int j=0; j<seed.size(); j++){
                    tempSim = max(tempSim, result[i*seed.size()+j]);
                }
                temp_user.sn = uIndex[i];
                temp_user.sim = tempSim;
                if (top_k.size() == k && temp_user.sim > top_k.top().sim && temp_user.sim != -1000.0)
                    top_k.pop();
                if (top_k.size() < k && temp_user.sim != -1000.0)
                    top_k.push(temp_user);

            }
            //释放申请的内存
            cudaFree(d_result);
            cudaFree(d_qIndex);
            cudaFree(d_uIndex);
            free(result);
            free(qIndex);
            free(uIndex);
        }
    }
    cudaFree(d_Udot);
    cudaFree(d_Sdot);
    cudaFree(d_data);
    cudaFree(d_query);
    cout<<"size of candicate user with uppbound"<<indexUser.size()<<endl;
    while(!indexUser.empty()){
        user = indexUser.top();
        indexUser.pop();
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


void gen_ExAudiences (int n_feats, int q, int n, float * data, float *queries){

    long query_Bytes = sizeof(float)*q*n_feats;
    long data_Bytes = sizeof(float)*n*n_feats;
    long result_Bytes = sizeof(float)*q*n;

    // 申请device内存
    float *result, *d_data, *d_query, *d_result;
    result = (float *)malloc(result_Bytes);
    cudaMalloc((void**)&d_query, query_Bytes);
    cudaMalloc((void**)&d_data, data_Bytes);
    cudaMalloc((void**)&d_result, result_Bytes);

    clock_t time1 = clock();
    cudaMemcpy((void*)d_query, (void*) queries, query_Bytes, cudaMemcpyHostToDevice);
    cudaMemcpy((void*)d_data, (void*) data, data_Bytes, cudaMemcpyHostToDevice);
    clock_t time2 = clock();

    dim3 blockSize(1024);
    cout<<n*q<<endl;
    dim3 gridSize((n*q + blockSize.x - 1) / blockSize.x);

    cout<<blockSize.x<<endl;
    cout<<gridSize.x<<endl;

    cout<<"tranfer data time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;

    //matrixMulCosine << < gridSize, blockSize >> >( d_data, d_query, d_result, q, n, n_feats);
    cudaMemcpy((void*)result, (void*)d_result, result_Bytes, cudaMemcpyDeviceToHost);
    clock_t time3 = clock();
    cout<<"matrixmul time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;
    for(int i=0; i<50; i++){
        cout<<result[i]<<endl;
    }
    cudaFree(d_data);
    cudaFree(d_query);
    cudaFree(d_result);
    free(result);
    // 释放host内存

}

__device__ float dotProduct(float* x, float * y, int x_index, int y_index, int n_feats){
    float temp = 0.0;
    x_index = x_index*n_feats;
    y_index = y_index*n_feats;
    for(int i=0; i<n_feats; i++){
        temp += x[x_index+i]*y[y_index+i];
    }
    return temp;
}

__global__ void calcuSimCuda(float* d_data, float * d_query, float * d_result, int *uIndex, int *qIndex, int q, int n_feats)
{

    // 获取索引
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int n_index = *uIndex;
    int q_index = qIndex[index];
    //printf("index: %d n_index: %d q_index: %d \r\n", index, uIndex[index/q], qIndex[index%q]);
    //关闭多余的线程
    if(index>=q)
        return;

    //初始化结果为0
    d_result[index]=0.0;
    float x,y,dot,temp;
    x = dotProduct(d_query, d_query, q_index, q_index, n_feats);
    y = dotProduct(d_data, d_data, n_index, n_index, n_feats);
    dot = dotProduct(d_data, d_query, n_index, q_index, n_feats);
    temp = dot/(sqrt(x)*sqrt(y));
/*
    printf("n_index: %d \r\n",n_index/50);
    printf("d_index: %d \r\n",q_index/50);
    if(index/q==0)
        printf("index: %d  temp: %f \r\n",index/q, temp);
 */
    d_result[index] = max(d_result[index], temp);

}

__global__ void matrixMulCosine(float* d_data, float * d_query, float * d_result, int *uIndex, int *qIndex, int q, int n, int n_feats)
{

    // 获取索引
    printf("n_index: %d \r\n",uIndex[0]);
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int n_index = uIndex[index/q];
    int q_index = qIndex[index%q];
    //printf("index: %d n_index: %d q_index: %d \r\n", index, uIndex[index/q], qIndex[index%q]);
    //关闭多余的线程
    if(index>=n*q)
        return;

    //初始化结果为0
    d_result[index]=0.0;
    float x,y,dot,temp;
    x = dotProduct(d_query, d_query, q_index, q_index, n_feats);
    y = dotProduct(d_data, d_data, n_index, n_index, n_feats);
    dot = dotProduct(d_data, d_query, n_index, q_index, n_feats);
    temp = dot/(sqrt(x)*sqrt(y));
/*
    printf("n_index: %d \r\n",n_index/50);
    printf("d_index: %d \r\n",q_index/50);
    if(index/q==0)
        printf("index: %d  temp: %f \r\n",index/q, temp);
 */
    d_result[index] = max(d_result[index], temp);

}

/*----------------------------------------------------------------*/
__global__ void preProcess(float* d_data, float * d_query, float * d_Udot, float * d_Sdot, int n, int q, int n_feats){
    int index;
    index = threadIdx.x + blockIdx.x * blockDim.x;
    if(index>n+q)
        return;

    float x;
    if(index<q){
        x = dotProduct(d_query, d_query, index, index, n_feats);
        d_Sdot[index]=sqrt(x);
    }else{
        index = index-q;
        x = dotProduct(d_data, d_data, index, index, n_feats);
        d_Udot[index]=sqrt(x);


    }
}

__global__ void matrixMulCosineOpt(float* d_data, float * d_query, float * d_result, int *d_uIndex, int *d_qIndex, int q, int n,
                                   float* d_Udot, float* d_Sdot,int n_feats)
{

    // 获取索引
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int n_index = d_uIndex[index/q];
    int q_index = d_qIndex[index%q];

    //关闭多余的线程
    if(index>=n*q)
        return;

    //初始化结果为0
    d_result[index]=0.0;
    float dot,temp;

    dot = dotProduct(d_data, d_query, n_index, q_index, n_feats);
    temp = dot/(d_Sdot[q_index]*d_Udot[n_index]);

    //printf("dot: %f ds: %f du: %f temp: %f\r\n",dot, d_Sdot[q_index], d_Udot[n_index], temp);
    d_result[index] = max(d_result[index], temp);

}
/*----------------------------------------------------------------*/
#if 0
int main(){
    float *data,*queries;
    data = (float *)malloc((int64_t)sizeof(float)*50*247753);
    csv_to_array(&data, "/home/andyshen/data/MovieLens/q.txt", 247753, 50);
    queries = (float *)malloc((int64_t)sizeof(float)*50*33670);
    csv_to_array(&queries, "/home/andyshen/data/MovieLens/p.txt", 33670, 50);
    cout<<"read data done"<<endl;

    clock_t time1 = clock();
    //float** sig_maritx = gen_signature_matrix(50, 5);
    float**  sig_maritx = gen_best_local_signature_matrix(50, 5, 50 , 33670, data, queries);
    vector<int> user_bucket_info;
    pre_user_pool(user_bucket_info, data, sig_maritx, 247753,50, 5);
    clock_t time2 = clock();
    unordered_map<int, vector<int>> user_maps_pool;
    user_map(user_maps_pool, queries, sig_maritx, 247753, 50, 5);
    unordered_map<int, vector<int>> user_maps_seed;

    user_map(user_maps_seed, queries, sig_maritx, 33670,50, 5);
    for(int i=0; i<32; i++){
        cout<<user_maps_seed[i].size()<<"  "<< user_maps_pool[i].size()<<endl;
    }
    clock_t time3 = clock();
    cout<<"init sigmatrix time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;
    cout<<"init lshmatrix time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;

    vector<bucket_info> centroid_angle;
    calculate_centroid_angle(centroid_angle, user_maps_seed, queries, 50, 5);

    priority_queue<canducate_user> top_k;
    priority_queue<uncertain_user> indexUser;
    //gen_ExAudiences_cuda(top_k, user_maps_pool, user_maps_seed, indexUser, centroid_angle,  5, 50, 1000, 247753, 33670, data, queries);
    gen_ExAudiences_cudabase(top_k, user_maps_seed, user_bucket_info,  indexUser, centroid_angle,  5, 50, 1000, 247753, 33670, data, queries);
    clock_t time4 = clock();
    cout<<"query time cuda: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;
    //gen_ExAudiences_cudaOpt(top_k, user_maps_pool, user_maps_seed, indexUser, centroid_angle,  5, 50, 1000, 247753, 33670, data, queries);
    priority_queue<uncertain_user> user_pool;
    gen_ExAudiences_vlsh(top_k, user_maps_pool, user_maps_seed, user_pool, centroid_angle, 5, 50, 1000, data, queries);
    clock_t time5 = clock();
    cout<<"query time cudaopt: "<<(time5-time4)/CLOCKS_PER_SECOND<<endl;
//    while(!top_k.empty()){
//        cout<<"NO.:"<<top_k.top().sn<<" Sim:"<<top_k.top().sim<<endl;
//        top_k.pop();
//    }

    return 0;


}
#endif