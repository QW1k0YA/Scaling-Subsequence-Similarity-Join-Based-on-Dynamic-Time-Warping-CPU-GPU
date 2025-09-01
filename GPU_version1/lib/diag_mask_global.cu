
#include "matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "GPU_parameters.h"
using namespace std;
#define BIAS 3

__device__ void
diag_mask_global(const FLOAT *ts, int subseqlen, int diagID, bool *lb_vector, const FLOAT *mu, const FLOAT *UTS,
                 const FLOAT *LTS, const FLOAT *MASK, const FLOAT *TS2, const FLOAT *sumMASK,
                 const FLOAT *invsig, const FLOAT *sumU_sumL, const FLOAT *dr_bwdU_plus_dr_bwdL,
                 const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_bwd, const FLOAT *dc_fwd,
                 const FLOAT *dr_bwdMASK, const FLOAT *dr_fwdMASK, const FLOAT *dc_bwdTS2,
                 const FLOAT *dc_fwdTS2, const FLOAT *DUL2, FLOAT *lb_vector_new, FLOAT  bsf, FLOAT **subs,
                 FLOAT **special_shared_vector, FLOAT &cnt, const FLOAT *DUL,
                 const FLOAT *norm_U_plus_norm_L_global, int start_pos, int end_pos)
{
    FLOAT  M_NORM,t,dist,DUL_value;

    int diag = diagID;
    FLOAT  cov_U_plus_cov_L_fir = (elementWiseMultiply_p_plus_sum(ts + start_pos + diag - 1 + 3,
                                                                  UTS + start_pos + 3,LTS + start_pos + 3,subseqlen - 6));
    FLOAT  dot_TS_M_fir = (elementWiseMultiply_p_sum(ts + start_pos  + diag - 1 + 3,
                                                     MASK + start_pos + 3,subseqlen - 6));
    FLOAT  dot_TS2_M_fir = (elementWiseMultiply_p_sum(TS2 + start_pos + diag - 1 + 3,
                                                      MASK + start_pos + 3,subseqlen - 6));
    FLOAT  cov_U_plus_cov_L_sec = (elementWiseMultiply_p_plus_sum(ts + start_pos + 3 , UTS+ start_pos  + diag - 1 + 3,
                                                                  LTS + start_pos  + diag - 1 + 3,subseqlen - 6));
    FLOAT  dot_TS_M_sec = (elementWiseMultiply_p_sum(ts + start_pos + 3,
                                                     MASK + start_pos  + diag - 1 + 3,subseqlen - 6));
    FLOAT  dot_TS2_M_sec = (elementWiseMultiply_p_sum(TS2 + start_pos + 3,
                                                      MASK + start_pos + diag - 1 + 3,subseqlen - 6));

    int row;
    int col;
    FLOAT  lb;
    if(!lb_vector[0])
    {
        row = start_pos;
        col = start_pos + diag - 1;

        M_NORM = dot_TS2_M_fir - 2 * mu[col] * dot_TS_M_fir + sumMASK[row] * mu[col] * mu[col];
        M_NORM = M_NORM * invsig[col] * invsig[col];

        t = (cov_U_plus_cov_L_fir - mu[col] * sumU_sumL[row] +
             2 * mu[row] * (sumMASK[row] * mu[col] - dot_TS_M_fir)) * (invsig[col] * invsig[row]);
        dist = 2 * M_NORM - 2 * t;

        if(dist < DUL2[row]){
            lb = 0;
        }
        else{
            lb = 0.5*(sqrt(2*dist - DUL2[row]) - DUL[row]);
        }
        lb_vector_new [0] = lb;
        if(lb > bsf){
            lb_vector[0] = true;
            cnt++;
        }
        else
        {
            col = start_pos;
            row = start_pos + diag - 1;

            M_NORM = dot_TS2_M_sec - 2 * mu[col] * dot_TS_M_sec + sumMASK[row] * mu[col] * mu[col];
            M_NORM = M_NORM * invsig[col] * invsig[col];

            t = (cov_U_plus_cov_L_sec - mu[col] * sumU_sumL[row] +
                 2 * mu[row] * (sumMASK[row] * mu[col] - dot_TS_M_sec)) * (invsig[col] * invsig[row]);
            dist = 2 * M_NORM - 2 * t;

            if(dist < DUL2[row]){
                lb = 0;
            }
            else {
                lb = 0.5 * (sqrt(2 * dist - DUL2[row]) - DUL[row]);
            }
            lb_vector_new [0] = lb;
            if(lb > bsf)
            {
                lb_vector[0] = true;
                cnt++;
            }
        }
    }

    int lb_pos = 0;
    for (int low_index = start_pos + 1; low_index < end_pos; low_index++) {
        int high_index = diag + low_index - 1;
        lb_pos++;

        cov_U_plus_cov_L_fir = cov_U_plus_cov_L_fir - dr_bwdU_plus_dr_bwdL[low_index ] * dc_bwd[high_index ] +
                               (dr_fwdU_plus_dr_fwdL[low_index ]) * dc_fwd[high_index];
        dot_TS_M_fir = dot_TS_M_fir - dr_bwdMASK[low_index] * dc_bwd[high_index] + dr_fwdMASK[low_index] * dc_fwd[high_index];
        dot_TS2_M_fir = dot_TS2_M_fir - dr_bwdMASK[low_index] * dc_bwdTS2[high_index] + dr_fwdMASK[low_index] * dc_fwdTS2[high_index];

        cov_U_plus_cov_L_sec = cov_U_plus_cov_L_sec - dr_bwdU_plus_dr_bwdL[high_index] * dc_bwd[low_index] +
                               (dr_fwdU_plus_dr_fwdL[high_index]) * dc_fwd[low_index];
        dot_TS_M_sec = dot_TS_M_sec - dr_bwdMASK[high_index] * dc_bwd[low_index] + dr_fwdMASK[high_index] * dc_fwd[low_index];
        dot_TS2_M_sec = dot_TS2_M_sec - dr_bwdMASK[high_index] * dc_bwdTS2[low_index] + dr_fwdMASK[high_index] * dc_fwdTS2[low_index];

        if(lb_vector[lb_pos])
        {
            continue;
        }

        M_NORM = dot_TS2_M_fir - 2 * mu[high_index] * dot_TS_M_fir + sumMASK[low_index] * mu[high_index] * mu[high_index];
        M_NORM = M_NORM * invsig[high_index] * invsig[high_index];
        t = (cov_U_plus_cov_L_fir - mu[high_index] * sumU_sumL[low_index] +
             2 * mu[low_index] * (sumMASK[low_index] * mu[high_index] - dot_TS_M_fir)) * (invsig[high_index] * invsig[low_index]);
        dist = 2 * M_NORM - 2 * t  + norm_U_plus_norm_L_global[low_index];

        if(dist < DUL2[low_index])
        {
            lb = 0;
        }
        else
        {

            lb = 0.5*(sqrt(2 * dist - DUL2[low_index]) - DUL[low_index]) ;

            lb_vector_new [lb_pos] = lb;
            if(threadIdx.x == 0)

            if(lb > bsf)
            {

                lb_vector[lb_pos] = true;

                cnt++;
                continue;
            }
        }

        M_NORM = dot_TS2_M_sec - 2 * mu[low_index] * dot_TS_M_sec + sumMASK[high_index] * mu[low_index] * mu[low_index];
        M_NORM = M_NORM * invsig[low_index] * invsig[low_index];
        t = (cov_U_plus_cov_L_sec - mu[low_index] * sumU_sumL[high_index] +
             2 * mu[high_index] * (sumMASK[high_index] * mu[low_index] - dot_TS_M_sec)) * (invsig[low_index] * invsig[high_index]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L_global[high_index] ;

        if(dist < DUL2[high_index])
        {
            lb = 0;
        }
        else {
            lb = 0.5 * (sqrt(2 * dist - DUL2[high_index]) - DUL[high_index]);
            lb_vector_new [lb_pos] = max(lb, lb_vector_new [lb_pos]);

            if(lb > bsf)
            {
                lb_vector[lb_pos] = true;
                cnt++;

            }
        }

    }

    lb_pos = -1;
    auto temp_1 = end_pos;
    for (col = start_pos; col < temp_1; col++)
    {
        lb_pos++;
        row = diag + col - 1;

        if(lb_vector[lb_pos])
        {
            continue;
        }

        const FLOAT *t_ = subs[row];
        const FLOAT *q = subs[col];

        int m = subseqlen;
        FLOAT d;
        FLOAT threshold = bsf;
        FLOAT threshold2= threshold*threshold;

        FLOAT x0 = t_[0] ;
        FLOAT y0 = t_[(m - 1 )] ;

        const FLOAT dleft_orgin=DIST(x0, q[0]);
        const FLOAT dright_orgin=DIST(y0, q[m - 1]);
        FLOAT dleft = dleft_orgin;
        FLOAT dright = dright_orgin;

        FLOAT x1 = (t_[( 1)] );
        const FLOAT d_left_weak = min(DIST(x1, q[0]), DIST(x1, q[1]));
        d = min(d_left_weak, DIST(x0, q[1]));
        dleft+=d;

        FLOAT y1 = (t_[(m - 2 )]);
        const FLOAT d_right_weak = min(DIST(y1, q[m - 1]),  DIST(y1, q[m - 2]));
        d = min(d_right_weak,DIST(y0, q[m - 2]));
        dright+=d;

        if (dleft+dright + lb_vector_new [lb_pos] * lb_vector_new [lb_pos]>=threshold2){
            cnt++;
            lb_vector[lb_pos] = true;
            continue;
        }
        else{

            d = MIN(DIST(x1,q[0]) + DIST(t_[2],q[0]),d_left_weak + DIST(t_[2],q[1]));
            d = MIN(d,DIST(q[1],t_[1]) + DIST(q[2],t_[2]));
            d = MIN(d,dleft + DIST(q[2],t_[1]));
            d = MIN(d,DIST(q[1],t_[0]) + DIST(q[2],t_[0]));
            dleft = d + dleft_orgin;

            d = MIN(DIST(t_[m-2],q[m-1]) + DIST(t_[m-3],q[m-1]),d_right_weak + DIST(t_[m-3],q[m-2]));
            d = MIN(d,DIST(q[m-2],t_[m-2]) + DIST(q[m-3],t_[m-3]));
            d = MIN(d,dright + DIST(q[m-3],t_[m-2]));
            d = MIN(d,DIST(q[m-2],t_[m-1]) + DIST(q[m-3],t_[m-1]));
            dright = d + dright_orgin;

            if (dleft+dright  + lb_vector_new [lb_pos]* lb_vector_new [lb_pos] >=threshold2){
                lb_vector[lb_pos] = true;
                cnt++;
                continue;
            }

        }

    }

}

