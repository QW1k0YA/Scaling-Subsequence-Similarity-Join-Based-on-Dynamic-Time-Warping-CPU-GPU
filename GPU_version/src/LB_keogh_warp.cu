
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
__device__ FLOAT lb_keogh_warp(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;

    for (size_t i = 0; i < seqlen; i++) {
        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
    }

    return dist;
}

/**
 * Device function for computing the Keogh lower bound in a strided manner.
 * This function calculates the LB_Keogh lower bound by having each thread process
 * elements at regular intervals (strides) to enable parallel computation across
 * threads in a warp/block.
 *
 * @param q: Query subsequence
 * @param U: Upper envelope of reference subsequence
 * @param L: Lower envelope of reference subsequence
 * @param seqlen: Length of the sequences
 * @return FLOAT: Computed lower bound distance
 */
__device__ FLOAT lb_keogh_stride(const FLOAT *q, const FLOAT *U, const FLOAT *L, long long seqlen) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;
    int tid = threadIdx.x;
    int step = blockDim.x;

    for (size_t i = tid; i < seqlen; i+=step) {
        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
    }

    return dist;
}

__device__ FLOAT lb_keogh_warp_with_nomalise(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen,const FLOAT* mu,const FLOAT* invsig) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;

    for (size_t i = 0; i < seqlen; i++) {
        FLOAT d = 0;
        temp_u = (U[i] - mu[i]) * invsig[i];
        temp_l = (L[i] - mu[i]) * invsig[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
    }

    return dist;
}

/**
 * Device function for computing the normalized Keogh lower bound in a strided manner.
 * This function calculates the LB_Keogh lower bound with normalization applied,
 * where each value is normalized using the mean (mu) and inverse standard deviation (invsig).
 *
 * @param q: Query subsequence
 * @param U: Upper envelope of reference subsequence
 * @param L: Lower envelope of reference subsequence
 * @param seqlen: Length of the sequences
 * @param mu: Mean value for normalization
 * @param invsig: Inverse of standard deviation for normalization
 * @return FLOAT: Computed normalized lower bound distance
 */
__device__ FLOAT lb_keogh_stride_with_nomalise(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen,
                                               const FLOAT mu, const FLOAT invsig) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;
    int tid = threadIdx.x;
    int step = blockDim.x;

    for (size_t i = tid; i < seqlen; i+=step) {
        FLOAT d = 0;
        temp_u = (U[i] - mu) * invsig;
        temp_l = (L[i] - mu) * invsig;
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
    }

    return dist;
}