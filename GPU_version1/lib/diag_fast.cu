
#include "iostream"
#include "ctime"
#include "vector"
#include "cmath"
#include "matrix.cuh"
#include "elseoperation.cuh"
#include "GPU_parameters.h"
#include "algorithm"
#include "chrono"
using namespace std;

__device__ void diag_fast_segment(const FLOAT *ts, int subseqlen, int diag_ID, const FLOAT *UTS, const FLOAT *LTS,
                                  const FLOAT *mu, const FLOAT *sumU_sumL, const FLOAT *invsig,
                                  const FLOAT *norm_U_plus_norm_L_trans, bool *lb_vector,
                                  const FLOAT *dr_bwdU_plus_dr_bwdL, const FLOAT *dc_bwd,
                                  const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_fwd, FLOAT &cnt,
                                  int start_pos, int end_pos, FLOAT bsf,
                                  const FLOAT *DUL, const FLOAT *DUL2) {

    int diag = diag_ID;
    FLOAT  lb;
    FLOAT  cov_U_plus_cov_L_fir = elementWiseMultiply_p_plus_sum(ts + start_pos + diag - 1,
                                                                UTS+ start_pos, LTS+ start_pos, subseqlen);
    FLOAT  cov_U_plus_cov_L_sec = (elementWiseMultiply_p_plus_sum(ts+ start_pos , UTS + start_pos + diag - 1,
                                                                  LTS + start_pos +  diag - 1, subseqlen));
    FLOAT fir_basic = cov_U_plus_cov_L_fir;
    FLOAT sec_basic = cov_U_plus_cov_L_sec;

    int row = start_pos;
    int col = start_pos + diag - 1;
    FLOAT  tt = (cov_U_plus_cov_L_fir - mu[col] * sumU_sumL[row]) * invsig[col];
    FLOAT  local_del = norm_U_plus_norm_L_trans[row] - 2 * tt;

    lb = 0.5*(sqrt(2*((local_del*invsig[row] + 2*subseqlen*mu[row]*mu[row]*invsig[row]*invsig[row]) + 2*subseqlen) - DUL2[row]) - DUL[row]);

    if (lb > bsf) {
        cnt++;
        lb_vector [row] = true;
    }
    else
    {
        col = start_pos;
        row = start_pos + diag - 1;
        tt = (cov_U_plus_cov_L_sec - mu[col] * sumU_sumL[row]) * invsig[col];
        local_del = norm_U_plus_norm_L_trans[row] - 2 * tt;
        lb = 0.5*(sqrt(2*((local_del + 2*subseqlen*mu[row]*mu[row]*invsig[row])*invsig[row] + 2*subseqlen) - DUL2[row]) - DUL[row]);

        if (lb>bsf) {
            lb_vector [col] = true;
            cnt++;
        }
    }

    for (int low_index = start_pos + 1; low_index < end_pos; low_index++) {

        int high_index = diag + low_index - 1;

        cov_U_plus_cov_L_fir = cov_U_plus_cov_L_fir -
                               dr_bwdU_plus_dr_bwdL[low_index] * dc_bwd[high_index] + dr_fwdU_plus_dr_fwdL[low_index] * dc_fwd[high_index];

        cov_U_plus_cov_L_sec = cov_U_plus_cov_L_sec -
                               dr_bwdU_plus_dr_bwdL[high_index] * dc_bwd[low_index] + dr_fwdU_plus_dr_fwdL[high_index] * dc_fwd[low_index];

        tt = (cov_U_plus_cov_L_fir - mu[high_index] * sumU_sumL[low_index]) * invsig[high_index];
        local_del = norm_U_plus_norm_L_trans[low_index] - 2 * tt;

        lb = 0.5*(sqrt(2*((local_del + 2*subseqlen*mu[low_index]*mu[low_index]*invsig[low_index])*invsig[low_index] + 2*subseqlen) - DUL2[low_index]) - DUL[low_index]);

        if (lb> bsf) {
            cnt++;
            lb_vector [low_index] = true;
            continue;
        }

        tt = (cov_U_plus_cov_L_sec - mu[low_index] * sumU_sumL[high_index]) * invsig[low_index];
        local_del = norm_U_plus_norm_L_trans[high_index] - 2 * tt;

        lb = 0.5*(sqrt(2*((local_del + 2*subseqlen*mu[high_index]*mu[high_index]*invsig[high_index])*invsig[high_index] + 2*subseqlen) - DUL2[high_index]) - DUL[high_index]);

        if (lb > bsf) {
            lb_vector [low_index] = true;
            cnt++;
        }

    }

}

