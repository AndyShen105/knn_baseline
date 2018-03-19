//
// Created by Hang Shen on 2018/3/15.
//


/* Return squared distance */
inline float query_distance(float *data[], float *query[], int i, int j, int r, int n, int q) {
    float sqDistance = 0.0f;
    for (int d = 0; d < r; d++) {
        double temp=(*data)[i*r+d] * (*query)[d*q+j];
        sqDistance+=temp;
    }
    return sqDistance;
}

/* Distances array size = n*q */
void serial_distances(float *data[], float *distances[], int n, int r, int k, float *query[], int q) {
    for (int i = 0; i < n; i++) { // For each object
        for (int d = 0; d < q; d++) { // For each query point
            (*distances)[i*q+d] = query_distance(data, query, i, d, r, n, q);
        }
    }
}

void knn_distance(float *data[], float *distances[], int n, int r, int k, float *query[], int q) {
    serial_distances(data, distances, n, r, k, query, q);

}
