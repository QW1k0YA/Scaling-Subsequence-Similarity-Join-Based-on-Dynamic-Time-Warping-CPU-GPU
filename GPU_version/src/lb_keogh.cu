
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
bool LB_KK(const vector<FLOAT> &q, const FLOAT *U, const FLOAT *L, long long seqlen, FLOAT threshold_2,
           vector<FLOAT> &cb, FLOAT special_shared_vector[], FLOAT &lbk) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;
    
    dist=special_shared_vector[0]+special_shared_vector[1];

    FLOAT temp_dist = 0;
    for (size_t i = 0; i < 3 ; i++) {

        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        temp_dist += d;
    }

    for (size_t i = seqlen - 3; i < seqlen ; i++) {
        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        temp_dist += d;
    }

    dist = MAX(temp_dist,dist);

    cb[0]=special_shared_vector[0];
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=special_shared_vector[1];
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    for (size_t i = 3; i < seqlen - 3 ; i++) {
        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }

        dist += d;
        if(threshold_2 < dist)
        {
            return 0;
        }
        cb[i] = d;
    }
    lbk = dist;
    return 1;
}

__device__ bool LB_KK(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen, FLOAT threshold_2,
                      FLOAT *cb, FLOAT special_shared_vector[], FLOAT &lbk) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;
    
    dist=special_shared_vector[0]+special_shared_vector[1];

    FLOAT temp_dist = 0;
    for (size_t i = 0; i < 3 ; i++) {

        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        temp_dist += d;
    }

    for (size_t i = seqlen - 3; i < seqlen ; i++) {
        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        temp_dist += d;
    }

    dist = MAX(temp_dist,dist);

    cb[0]=special_shared_vector[0];
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=special_shared_vector[1];
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    for (size_t i = 3; i < seqlen - 3 ; i++) {
        FLOAT d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }

        dist += d;
        if(threshold_2 < dist)
        {
            return 0;
        }
        cb[i] = d;
    }
    lbk = dist;
    return 1;
}

/**
 * Implements the first level of the Keogh-Koutra lower bound (LB_Keogh) optimization.
 * This function computes a fast lower bound estimation for DTW distance to prune
 * unnecessary full DTW computations. It checks early termination conditions based
 * on the first few and last few points of the sequences.
 *
 * @param q: Query subsequence
 * @param t: Target subsequence
 * @param U: Upper envelope of the target subsequence
 * @param L: Lower envelope of the target subsequence
 * @param seqlen: Length of the sequences
 * @param threshold_2: Squared threshold for early termination
 * @param cb: Cost buffer to store intermediate distances
 * @return bool: True if the lower bound is below threshold (further computation needed), False otherwise
 */
bool LB_KK_FIRST(const vector<FLOAT> &q, const vector<FLOAT> &t, const vector<FLOAT> &U, const vector<FLOAT> &L,
                 long long seqlen, FLOAT threshold_2, vector<FLOAT> &cb) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;

    int m = seqlen;
    FLOAT d = 0;
    FLOAT x0 = t[0] ;
    FLOAT y0 = t[(m - 1 )] ;
    FLOAT dleft=DIST(x0, q[0]);
    FLOAT dright=DIST(y0, q[m - 1]);

    FLOAT x1 = (t[( 1)] );
    d = min(DIST(x1, q[0]), DIST(x0, q[1]));
    d = min(d, DIST(x1, q[1]));
    dleft+=d;

    FLOAT y1 = (t[(m - 2 )]);
    d = min(DIST(y1, q[m - 1]), DIST(y0, q[m - 2]) );
    d = min(d, DIST(y1, q[m - 2]));
    dright+=d;

    if (dleft+dright>=threshold_2){
        return false;
    }
    else{
        FLOAT x2 = (t[(2)]);
        d = min(DIST(x0, q[2]), DIST(x1, q[2]));
        d = min(d, DIST(x2, q[2]));
        d = min(d, DIST(x2, q[1]));
        d = min(d, DIST(x2, q[0]));
        dleft += d;

        FLOAT y2 = (t[(m - 3 )]);
        d = min(DIST(y0, q[m - 3]), DIST(y1, q[m - 3]));
        d = min(d, DIST(y2, q[m - 3]));
        d = min(d, DIST(y2, q[m - 2]));
        d = min(d, DIST(y2, q[m - 1]));
        dright += d;

        if (dleft+dright> threshold_2){
            return false;
        }
    }

    cb[0]=dleft;
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=dright;
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    dist = dleft + dright;
    for (size_t i = 3; i < seqlen - 3; i++) {
        d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
        if(threshold_2 < dist)
        {
            return false;
        }
        cb[i] = d;
    }
    return true;
}

bool LB_KK_FIRST(const FLOAT* q, const FLOAT* t, const FLOAT* U, const FLOAT* L,
                 long long seqlen, FLOAT threshold_2, vector<FLOAT> &cb) {

    FLOAT  dist = 0.0;
    FLOAT temp_u;
    FLOAT temp_l;

    int m = seqlen;
    FLOAT d = 0;
    FLOAT x0 = t[0] ;
    FLOAT y0 = t[(m - 1 )] ;
    FLOAT dleft=DIST(x0, q[0]);
    FLOAT dright=DIST(y0, q[m - 1]);

    FLOAT x1 = (t[( 1)] );
    d = min(DIST(x1, q[0]), DIST(x0, q[1]));
    d = min(d, DIST(x1, q[1]));
    dleft+=d;

    FLOAT y1 = (t[(m - 2 )]);
    d = min(DIST(y1, q[m - 1]), DIST(y0, q[m - 2]) );
    d = min(d, DIST(y1, q[m - 2]));
    dright+=d;

    if (dleft+dright>=threshold_2){
        return false;
    }
    else{
        FLOAT x2 = (t[(2)]);
        d = min(DIST(x0, q[2]), DIST(x1, q[2]));
        d = min(d, DIST(x2, q[2]));
        d = min(d, DIST(x2, q[1]));
        d = min(d, DIST(x2, q[0]));
        dleft += d;

        FLOAT y2 = (t[(m - 3 )]);
        d = min(DIST(y0, q[m - 3]), DIST(y1, q[m - 3]));
        d = min(d, DIST(y2, q[m - 3]));
        d = min(d, DIST(y2, q[m - 2]));
        d = min(d, DIST(y2, q[m - 1]));
        dright += d;

        if (dleft+dright> threshold_2){
            return false;
        }
    }

    cb[0]=dleft;
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=dright;
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    dist = dleft + dright;
    for (size_t i = 3; i < seqlen - 3; i++) {
        d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
        if(threshold_2 < dist)
        {
            return false;
        }
        cb[i] = d;
    }
    return true;
}