
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