
#include "main.h"


int main(int argc, const char * argv[]) {

    /* Parameters */
    int n, r, q, k, n_samples, flag;
    float *data, *queries;
    n = atoi( argv[1] ); // Size of dataset (rows)
    r = atoi( argv[2] ); // Number of dimensions (columns)
    k = atoi( argv[3] ); // Nearest neighbors
    q = atoi( argv[4] ); // Number of queries
    flag = atoi( argv[5] ); //query by which mode
    n_samples = atoi( argv[6] ); //query by which mode
    string inputFile = argv[7];
    string queryFile = argv[8];

    cout<<"Size of dataset: "<<n<<endl;
    cout<<"Number of dimensions : "<<r<<endl;
    cout<<"Top-k: "<<k<<endl;
    cout<<"Number of Seeds: "<<q<<endl;
    cout<<"State: "<<flag<<endl;
    cout<<"n_sample: "<<n_samples<<endl;
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


    //gen sig maritx
    int n_bit = 5;
    priority_queue<canducate_user> top_k;

    switch  (flag){
        case 0:
        {   //brute force
            clock_t timeStart2 = clock();
            knn_distance(&data, top_k, n, r, &queries, q, k);
            clock_t timeEnd2 = clock();
            double processTime = (timeEnd2-timeStart2)/CLOCKS_PER_SECOND ;
            cout<<"Duration of process time: "<< processTime <<endl;
            break;
        }
        case 1:
        {   //lsh
            float** sig_maritx = gen_signature_matrix(50, n_bit);
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, data, sig_maritx, n, r, n_bit);
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx, q, r, n_bit);

            clock_t time1 = clock();
            gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, n_bit, r, k, data, queries);
            clock_t time2 = clock();
            cout<<"lsh query time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;

            break;
        }
        case 2:
        {   //lsh+suiji
            float**  sig_maritx_best = gen_best_local_signature_matrix(50, 5, n_samples , q, data, queries);
            clock_t time1 = clock();
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx_best, q, r, n_bit);
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, data, sig_maritx_best, n, r, n_bit);
            clock_t time2 = clock();


            vector<bucket_info> centroid_angle;
            priority_queue<uncertain_user> user_pool;
            clock_t time3 = clock();
            cout<<"pre time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;
            gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, n_bit, r, k, data, queries);
            clock_t time4 = clock();
            cout<<"lsh+suiji query time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;
            break;
        }
        case 3:
        {   //vlsh
            float**  sig_maritx_best = gen_best_local_signature_matrix(50, 5, n_samples, q, data, queries);
            clock_t time1 = clock();
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx_best, q, r, n_bit);
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, data, sig_maritx_best, n, r, n_bit);
            clock_t time2 = clock();
            //v-lsh
            vector<int > v_lsh_re;
            vector<bucket_info> centroid_angle;
            priority_queue<uncertain_user> user_pool;
            calculate_centroid_angle(centroid_angle, user_maps_seed, queries, r, n_bit);
            clock_t time3 = clock();
            cout<<"pre time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;
            gen_ExAudiences_vlsh(top_k, user_maps_pool, user_maps_seed, user_pool, centroid_angle, n_bit, r, k, data, queries);
            clock_t time4 = clock();
            cout<<"v_lsh query time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;
            break;
        }

        case 4:
        {   //vlsh and lsh
            clock_t time1 = clock();
            float**  sig_maritx_best = gen_best_local_signature_matrix(50, 5, n_samples , q, data, queries);
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx_best, q, r, n_bit);
            clock_t time2 = clock();
            cout<<"pre-process time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;

            //init user and see map which used by lsh and vlsh
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, data, sig_maritx_best, n, r, n_bit);


            clock_t time3 = clock();
            //query with lsh
            gen_ExAudiences(top_k, user_maps_pool, user_maps_seed, n_bit, r, k, data, queries);
            clock_t time4 = clock();
            cout<<"lsh query time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;
            vector<int> lsh_re;
            while(!top_k.empty()){
                lsh_re.push_back(top_k.top().sn);
                top_k.pop();
            }

            //v-lsh
            vector<int > v_lsh_re;
            vector<bucket_info> centroid_angle;
            // init centroid info
            clock_t time5 = clock();
            calculate_centroid_angle(centroid_angle, user_maps_seed, queries, r, n_bit);
            clock_t time6 = clock();
            cout<<"calculate centroid angle time: "<<(time6-time5)/CLOCKS_PER_SECOND<<endl;

            //query by vlsh

            priority_queue<uncertain_user> user_pool;
            //ProfilerStart("vlsh.prof");
            gen_ExAudiences_vlsh(top_k, user_maps_pool, user_maps_seed, user_pool, centroid_angle, n_bit, r, k, data, queries);
            //ProfilerStop();
            clock_t time7 = clock();
            cout<<"v_lsh query time: "<<(time7-time6)/CLOCKS_PER_SECOND<<endl;
            while(!top_k.empty()){
                v_lsh_re.push_back(top_k.top().sn);
                top_k.pop();
            }

            cout<<"size of vlsh: "<<v_lsh_re.size()<<endl;
            cout<<"size of  lsh: "<<lsh_re.size()<<endl;
            //calculate accary between lsh and vlsh
            cout<<"accuary: "<<calutate_acc(lsh_re, v_lsh_re)<<endl;

            break;
        }
        case 5:
        {
            clock_t time1 = clock();
            float**  sig_maritx_best = gen_best_local_signature_matrix(r, n_bit, n_samples , q, data, queries);
            unordered_map<int, vector<int>> user_maps_seed;
            user_map(user_maps_seed, queries, sig_maritx_best, q, r, n_bit);
            vector<int> user_bucket_info;
            pre_user_pool(user_bucket_info, data, sig_maritx_best, n, r, n_bit);
            clock_t time2 = clock();
            cout<<"pre-process time: "<<(time2-time1)/CLOCKS_PER_SECOND<<endl;

            //init user and see map which used by lsh and vlsh
            unordered_map<int, vector<int>> user_maps_pool;
            user_map(user_maps_pool, data, sig_maritx_best, n, r, n_bit);
            clock_t time3 = clock();
            cout<<"calculate user pool time: "<<(time3-time2)/CLOCKS_PER_SECOND<<endl;

            vector<bucket_info> centroid_angle;
            calculate_centroid_angle(centroid_angle, user_maps_seed, queries, r, n_bit);
            clock_t time4 = clock();
            cout<<"calculate centroid angle time: "<<(time4-time3)/CLOCKS_PER_SECOND<<endl;

            priority_queue<canducate_user> top_k;
            priority_queue<uncertain_user> indexUser;

            gen_ExAudiences_vlsh(top_k, user_maps_pool, user_maps_seed, indexUser, centroid_angle, n_bit, r, k, data, queries);
            clock_t time5 = clock();
            cout<<"cpu query time: "<<(time5-time4)/CLOCKS_PER_SECOND<<endl;

            //gen_ExAudiences_cuda(top_k, user_maps_pool, user_maps_seed, indexUser, centroid_angle,  5, 50, 1000, 247753, 33670, data, queries);
            gen_ExAudiences_cudabase(top_k, user_maps_seed, user_bucket_info,  indexUser, centroid_angle, n_bit, r, k, n, q, data, queries);
            clock_t time6 = clock();
            cout<<"query time cuda: "<<(time6-time5)/CLOCKS_PER_SECOND<<endl;

            gen_ExAudiences_cudaOpt(top_k, user_maps_pool, user_maps_seed, indexUser, centroid_angle, n_bit, r, k, n, q, data, queries);
            priority_queue<uncertain_user> user_pool;
            clock_t time7 = clock();
            cout<<"query time cudaopt: "<<(time7-time6)/CLOCKS_PER_SECOND<<endl;
            break;
        }

    }
    free(data); data = nullptr;
    free(queries); queries = nullptr;

    return EXIT_SUCCESS;
}