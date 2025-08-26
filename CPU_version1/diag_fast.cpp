
#include "iostream"
#include "ctime"
#include "vector"
#include "cmath"
#include "alldef/matrix.h"
#include "alldef/elseoperation.h"
#include "algorithm"
#include "chrono"
using namespace std;

void
diag_fast(const vector<DOUBLE> &ts, int subseqlen, int diag_ID, const vector<DOUBLE> &UTS, const vector<DOUBLE> &LTS,
          const vector<DOUBLE> &mu, const vector<DOUBLE> &sumU_sumL, const vector<DOUBLE> &invsig,
          const vector<DOUBLE> &norm_U_plus_norm_L_trans, const vector<DOUBLE> &del, vector<bool> &lb_vector,
          const vector<DOUBLE> &dr_bwdU_plus_dr_bwdL, const vector<DOUBLE> &dc_bwd,
          const vector<DOUBLE> &dr_fwdU_plus_dr_fwdL, const vector<DOUBLE> &dc_fwd, double &cnt,
          const vector<DOUBLE> &DUL, const vector<DOUBLE> &DUL2, double bsf) {

    int subcount = ts.size() - subseqlen + 1;
    int diag = diag_ID;
    DOUBLE lb;

    DOUBLE cov_U_plus_cov_L_fir = elementWiseMultiply_p_plus_sum(ts.data() + diag - 1,
                                                          UTS.data(), LTS.data(), subseqlen);
    DOUBLE cov_U_plus_cov_L_sec = (elementWiseMultiply_p_plus_sum(ts.data() , UTS.data() + diag - 1,
                                                                  LTS.data()+ diag - 1, subseqlen));
    DOUBLE fir_basic = cov_U_plus_cov_L_fir;
    DOUBLE sec_basic = cov_U_plus_cov_L_sec;
    int row = 0;
    int col = diag - 1;
    DOUBLE tt = (cov_U_plus_cov_L_fir - mu[col] * sumU_sumL[row]) * invsig[col];
    DOUBLE local_del = norm_U_plus_norm_L_trans[row] - 2 * tt;

    lb = 0.5*(sqrt(2*((local_del*invsig[row] + 2*subseqlen*mu[row]*mu[row]*invsig[row]*invsig[row]) + 2*subseqlen) - DUL2[row]) - DUL[row]);

    if (lb > bsf) {
        cnt++;
        lb_vector [row] = true;
    }
    else
    {
        col = 0;
        row = diag - 1;
        tt = (cov_U_plus_cov_L_sec - mu[col] * sumU_sumL[row]) * invsig[col];
        local_del = norm_U_plus_norm_L_trans[row] - 2 * tt;
        lb = 0.5*(sqrt(2*((local_del + 2*subseqlen*mu[row]*mu[row]*invsig[row])*invsig[row] + 2*subseqlen) - DUL2[row]) - DUL[row]);

        if (lb>bsf) {
            lb_vector [col] = true;
            cnt++;
        }
    }

    for (int low_index = 1; low_index < subcount - diag + 1; low_index++) {

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