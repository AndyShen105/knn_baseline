#include <iostream>
#include <iomanip>
#include "main.h"
#include "math.h"
/* Heap sorting */

/* Restricted log2: >0 values only */
inline unsigned int ulog2(unsigned int value) {
    unsigned int result = 0;
    while (value) {
        value >>= 1;
        result++;
    }
    return result - 1;
}

inline unsigned int heap_parent(unsigned int node) {
    return (node - 1) / 2;
}

inline unsigned int heap_leftchild(unsigned int node) {
    return 2 * node + 1;
}

inline unsigned int heap_rightchild(unsigned int node) {
    return 2 * node + 2;
}

inline bool heap_isminlevel(unsigned int node) {
    return ulog2(node + 1) % 2 == 0;
}

inline bool heap_ismaxlevel(unsigned int node) {
    return !heap_isminlevel(node);
}

//return a
inline unsigned int heap_maxidx(float *heapData[], unsigned int hsize) {
    //std::cout<<(*heapData)[1]<<std::endl;
    //std::cout<<(*heapData)[2]<<std::endl;
    //std::cout<<heapData[0]<<std::endl;
    if (hsize == 1) return 0;
    else if (hsize == 2) return 1;
    else return ((*heapData)[1] > (*heapData)[2] ? 1 : 2);
}

inline void heap_swap(float *heapData[], int *heapIndex[], int a, int b) {
    float tmpVal = (*heapData)[a];
    int tmpIdx = (*heapIndex)[a];
    (*heapData)[a] = (*heapData)[b];
    (*heapIndex)[a] = (*heapIndex)[b];
    (*heapData)[b] = tmpVal;
    (*heapIndex)[b] = tmpIdx;
}

void heap_tricklemin(float *heapData[], int *heapIndex[], int start, int end) {
    int root = start;
    int minChildIdx, minGChildIdx;
    float minChildVal, minGChildVal;

    trickle:
    if (heap_leftchild(root) < end) { // root has at least 1 child
        /* Check children */
        minChildIdx = heap_leftchild(root);
        minChildVal = (*heapData)[minChildIdx];

        if (heap_rightchild(root) < end && minChildVal > (*heapData)[heap_rightchild(root)]) {
            minChildIdx = heap_rightchild(root);
            minChildVal = (*heapData)[minChildIdx];
        }

        /* Check grandchildren */
        if (heap_leftchild(heap_leftchild(root)) < end) {
            minGChildIdx = heap_leftchild(heap_leftchild(root));
            minGChildVal = (*heapData)[minGChildIdx];

            for (int i = heap_leftchild(heap_leftchild(root)) + 1; (i < end) && (i < heap_leftchild(heap_leftchild(root)) + 4); i++) {
                if (minGChildVal > (*heapData)[i]) {
                    minGChildIdx = i;
                    minGChildVal = (*heapData)[i];
                }
            }
        }
        if (heap_leftchild(heap_leftchild(root)) >= end || minChildVal <= minGChildVal) {
            if (minChildVal < (*heapData)[root]) { // Swap parent and child
                heap_swap(heapData, heapIndex, root, minChildIdx);
            }
        }
        else { // Manage grandchildren
            if (minGChildVal < (*heapData)[root]) {
                heap_swap(heapData, heapIndex, root, minGChildIdx);
                if (minGChildVal > (*heapData)[heap_parent(minGChildIdx)]) {
                    heap_swap(heapData, heapIndex, heap_parent(minGChildIdx), minGChildIdx);
                }
                goto trickle;
            }
        }
    }
}

void heap_tricklemax(float *heapData[], int *heapIndex[], int start, int end) {
    int root = start;
    int maxChildIdx, maxGChildIdx;
    float maxChildVal, maxGChildVal;

    trickle:
    if (heap_leftchild(root) < end) { // root has at least 1 child
        /* Check children */
        maxChildIdx = heap_leftchild(root);
        maxChildVal = (*heapData)[maxChildIdx];

        if (heap_rightchild(root) < end && maxChildVal < (*heapData)[heap_rightchild(root)]) {
            maxChildIdx = heap_rightchild(root);
            maxChildVal = (*heapData)[maxChildIdx];
        }

        /* Check grandchildren */
        if (heap_leftchild(heap_leftchild(root)) < end) {
            maxGChildIdx = heap_leftchild(heap_leftchild(root));
            maxGChildVal = (*heapData)[maxGChildIdx];

            for (int i = heap_leftchild(heap_leftchild(root)) + 1; (i < end) && (i < heap_leftchild(heap_leftchild(root)) + 4); i++) {
                if (maxGChildVal < (*heapData)[i]) {
                    maxGChildIdx = i;
                    maxGChildVal = (*heapData)[i];
                }
            }
        }
        if (heap_leftchild(heap_leftchild(root)) >= end || maxChildVal >= maxGChildVal) {
            if (maxChildVal > (*heapData)[root]) { // Swap parent and child
                heap_swap(heapData, heapIndex, root, maxChildIdx);
            }
        }
        else { // Manage grandchildren
            if (maxGChildVal > (*heapData)[root]) {
                heap_swap(heapData, heapIndex, root, maxGChildIdx);
                if (maxGChildVal < (*heapData)[heap_parent(maxGChildIdx)]) {
                    heap_swap(heapData, heapIndex, heap_parent(maxGChildIdx), maxGChildIdx);
                }
                goto trickle;
            }
        }
    }
}

/* Turn an array into a min-max heap in height-linear time */
void heapify_minmax(float *heapData[], int *heapIndex[], int hsize) {
    int start = (int)floor((float)(hsize - 2) / 2);
    DBG(std::cout << "n of level of heap:"<< start << std::endl;);
    for (int i = start; i >= 0; i--) {
        if (heap_isminlevel(i)) {
            heap_tricklemin(heapData, heapIndex, i, hsize);
        }
        else {
            heap_tricklemax(heapData, heapIndex, i, hsize);
        }
    }
}

/* Insert into fixed size heap by overwriting largest element */
void heap_insert(float *heapData[], int *heapIndex[], float newPt, int newIdx, int hsize) {
    /* Overwrite largest element if applicable */
    int maxIdx = heap_maxidx(heapData, hsize), maxChildIdx;
    float maxChildVal;
    if ((*heapData)[maxIdx] > newPt) {
        (*heapData)[maxIdx] = newPt;
        (*heapIndex)[maxIdx] = newIdx;

        /* Compare with root and maintain ordering */
        if (newPt < (*heapData)[0]) { // New smallest point
            heap_swap(heapData, heapIndex, 0, maxIdx);
            heap_tricklemin(heapData, heapIndex, 0, hsize);
        }
        heap_tricklemax(heapData, heapIndex, maxIdx, hsize);

        /* Compare with children of sibling */
        maxChildIdx = (maxIdx == 1 ? 5 : 3);
        if (maxChildIdx < hsize && (maxChildIdx+1) < hsize) {
            maxChildVal = (*heapData)[maxChildIdx];

            if (maxChildVal < (*heapData)[maxChildIdx + 1]) {
                maxChildIdx = maxChildIdx + 1;
                maxChildVal = (*heapData)[maxChildIdx];
            }

            if (maxChildVal > (*heapData)[maxIdx]) {
                heap_swap(heapData, heapIndex, maxIdx, maxChildIdx);
            }
        }
    }
}

void heap_delete(float *heapData[], int *heapIndex[], int hsize) {
    /* Replace root with last element */
    (*heapData)[0] = (*heapData)[hsize - 1];
    (*heapIndex)[0] = (*heapIndex)[hsize - 1];

    heap_tricklemin(heapData, heapIndex, 0, hsize);
}

/* Partial heapsort with separate heap */
void serial_minmaxheapsort_index(float *data[], int *index[], int n, int q, int k) {
    /* Take first k elements of data and heapify */
    float *heapData = (float *)malloc(sizeof(float)*k); // Cache
    int *heapIndex = (int *)malloc(sizeof(int)*k);

    for (int d = 0; d < q; d++) { // Each query
        for (int i = 0; i < k; i++) {
            heapData[i] = (*data)[i*q+d];
            heapIndex[i] = i;
        }

        heapify_minmax(&heapData, &heapIndex, k);

        /* Iteratively insert remaining heap elements, starting at k */
        for (int j = k; j < n; j++) {
            heap_insert(&heapData, &heapIndex, (*data)[j*q+d], j, k);
        }

        /* Output */
        int heapIdx = 0;
        for (int i = 0; i < k; i++) {
            (*index)[i*q+d] = heapIndex[0];
            heap_delete(&heapData, &heapIndex, k-i);
        }
    }
}

/* Return max only */
void serial_minmaxheapsort_max(float *data[], float *maxval[], int n, int q, int k) {
    /* Take first k elements of data and heapify */
    float *heapData = (float *)malloc(sizeof(float)*k); // Cache
    int *heapIndex = (int *)malloc(sizeof(int)*k); //record the index of top-k
    for (int d = 0; d < q; d++) { // Each query
        for (int i = 0; i < k; i++) {
            heapData[i] = (*data)[i*q+d];
            heapIndex[i] = i;
        }

        heapify_minmax(&heapData, &heapIndex, k);


        /* Iteratively insert remaining heap elements, starting at k */
        for (int j = k; j < n; j++) {
            heap_insert(&heapData, &heapIndex, (*data)[j*q+d], j, k);

        }

        /* Output the max value of min-value list*/
        (*maxval)[d] = heapData[heap_maxidx(&heapData, k)];
        DBG(std::cout<<"key value"<<maxval[0][0]<<std::endl;);
    }
}

/* Return only maximum index (for threshold) */
void mksort(float *data[], float *maxval[], int n, int q, int k) {
        DBG(std::cout << "Sort [minmax heapsort] with CPU." << std::endl;);
        serial_minmaxheapsort_max(data, maxval, n, q, k);

}
