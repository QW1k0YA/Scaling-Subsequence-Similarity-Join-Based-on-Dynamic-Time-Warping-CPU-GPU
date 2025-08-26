
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define BIAS 3

void
diag_mask_global(const vector<DOUBLE> &ts, int subseqlen, int diagID, vector<bool> &lb_vector, const vector<DOUBLE> &mu,
                 const vector<DOUBLE> &UTS, const vector<DOUBLE> &LTS, const vector<DOUBLE> &MASK,
                 const vector<DOUBLE> &TS2, const vector<DOUBLE> &sumMASK, const vector<DOUBLE> &invsig,
                 const vector<DOUBLE> &sumU_sumL, const vector<DOUBLE> &dr_bwdU_plus_dr_bwdL,
                 const vector<DOUBLE> &dr_fwdU_plus_dr_fwdL, const vector<DOUBLE> &dc_bwd, const vector<DOUBLE> &dc_fwd,
                 const vector<DOUBLE> &dr_bwdMASK, const vector<DOUBLE> &dr_fwdMASK, const vector<DOUBLE> &dc_bwdTS2,
                 const vector<DOUBLE> &dc_fwdTS2, const vector<DOUBLE> &DUL2, vector<DOUBLE> &lb_vector_new, DOUBLE bsf,
                 const vector<vector<DOUBLE>> &subs, double **special_shared_vector, double &cnt,
                 const vector<DOUBLE> &norm_U_plus_norm_L_global, vector<DOUBLE> &DUL)
{

    size_t len = ts.size();
    size_t subcount = len - subseqlen + 1;

    DOUBLE M_NORM,t,dist,DUL_value;

    int diag = diagID;

    DOUBLE cov_U_plus_cov_L_fir = (elementWiseMultiply_p_plus_sum(ts.data() + diag - 1 + 3,
                                                              UTS.data() + 3,LTS.data() + 3,subseqlen - 6));
    DOUBLE dot_TS_M_fir = (elementWiseMultiply_p_sum(ts.data() + diag - 1 + 3,
                                                 MASK.data() + 3,subseqlen - 6));
    DOUBLE dot_TS2_M_fir = (elementWiseMultiply_p_sum(TS2.data() + diag - 1 + 3,
                                                  MASK.data() + 3,subseqlen - 6));
    DOUBLE cov_U_plus_cov_L_sec = (elementWiseMultiply_p_plus_sum(ts.data() + 3 , UTS.data() + diag - 1 + 3,
                                                                  LTS.data()+ diag - 1 + 3,subseqlen - 6));
    DOUBLE dot_TS_M_sec = (elementWiseMultiply_p_sum(ts.data() + 3,
                                                     MASK.data() + diag - 1 + 3,subseqlen - 6));
    DOUBLE dot_TS2_M_sec = (elementWiseMultiply_p_sum(TS2.data() + 3,
                                                      MASK.data() + diag - 1 + 3,subseqlen - 6));

    int row;
    int col;
    DOUBLE lb;
    if(!lb_vector[0])
    {
        row = 0;
        col = diag - 1;

        M_NORM = dot_TS2_M_fir - 2 * mu[col] * dot_TS_M_fir + sumMASK[row] * mu[col] * mu[col];
        M_NORM = M_NORM * invsig[col] * invsig[col];

        t = (cov_U_plus_cov_L_fir - mu[col] * sumU_sumL[row] +
             2 * mu[row] * (sumMASK[row] * mu[col] - dot_TS_M_fir)) * (invsig[col] * invsig[row]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L_global[row];

        if(dist < DUL2[row]){
            lb = 0;
        }
        else{
            lb = 0.5*(sqrt(2*dist - DUL2[row]) - DUL[row]);
        }
        lb_vector_new [row] = lb;
        if(lb > bsf){
            lb_vector[row] = true;
            cnt++;
        }
        else
        {
            col = 0;
            row = diag - 1;

            M_NORM = dot_TS2_M_sec - 2 * mu[col] * dot_TS_M_sec + sumMASK[row] * mu[col] * mu[col];
            M_NORM = M_NORM * invsig[col] * invsig[col];

            t = (cov_U_plus_cov_L_sec - mu[col] * sumU_sumL[row] +
                 2 * mu[row] * (sumMASK[row] * mu[col] - dot_TS_M_sec)) * (invsig[col] * invsig[row]);
            dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L_global[col];

            if(dist < DUL2[row]){
                lb = 0;
            }
            else {
                lb = 0.5 * (sqrt(2 * dist - DUL2[row]) - DUL[row]);
            }
            lb_vector_new [col] = lb;
            if(lb > bsf)
            {
                lb_vector[col] = true;
                cnt++;
            }
        }
    }

    for (int low_index = 1; low_index < subcount - diag + 1; low_index++) {

        int high_index = diag + low_index - 1;

        cov_U_plus_cov_L_fir = cov_U_plus_cov_L_fir - dr_bwdU_plus_dr_bwdL[low_index ] * dc_bwd[high_index ] +
                               (dr_fwdU_plus_dr_fwdL[low_index ]) * dc_fwd[high_index];
        dot_TS_M_fir = dot_TS_M_fir - dr_bwdMASK[low_index] * dc_bwd[high_index] + dr_fwdMASK[low_index] * dc_fwd[high_index];
        dot_TS2_M_fir = dot_TS2_M_fir - dr_bwdMASK[low_index] * dc_bwdTS2[high_index] + dr_fwdMASK[low_index] * dc_fwdTS2[high_index];

        cov_U_plus_cov_L_sec = cov_U_plus_cov_L_sec - dr_bwdU_plus_dr_bwdL[high_index] * dc_bwd[low_index] +
                               (dr_fwdU_plus_dr_fwdL[high_index]) * dc_fwd[low_index];
        dot_TS_M_sec = dot_TS_M_sec - dr_bwdMASK[high_index] * dc_bwd[low_index] + dr_fwdMASK[high_index] * dc_fwd[low_index];
        dot_TS2_M_sec = dot_TS2_M_sec - dr_bwdMASK[high_index] * dc_bwdTS2[low_index] + dr_fwdMASK[high_index] * dc_fwdTS2[low_index];

        if(lb_vector[low_index])
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
            lb_vector_new [low_index] = lb;

            if(lb > bsf)
            {
                lb_vector[low_index] = true;
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
            lb_vector_new [low_index] = max(lb, lb_vector_new [low_index]);
            if(lb > bsf)
            {
                lb_vector[low_index] = true;
                cnt++;
            }
        }

    }

    auto temp_1 = subcount - diag + 1;
    for (col = 0; col < temp_1; col++)
    {

        row = diag + col - 1;

        if(lb_vector[col])
        {
            continue;
        }
        
        const vector<DOUBLE> &t_ = subs[row];
        const vector<DOUBLE> &q = subs[col];

        int m = subseqlen;
        double d;
        double threshold = bsf;
        double threshold2= threshold*threshold;

        double x0 = t_[0] ;
        double y0 = t_[(m - 1 )] ;

        const double dleft_orgin=DIST(x0, q[0]);
        const double dright_orgin=DIST(y0, q[m - 1]);
        double dleft = dleft_orgin;
        double dright = dright_orgin;

        double x1 = (t_[( 1)] );
        const double d_left_weak = min(DIST(x1, q[0]), DIST(x1, q[1]));
        d = min(d_left_weak, DIST(x0, q[1]));
        dleft+=d;

        double y1 = (t_[(m - 2 )]);
        const double d_right_weak = min(DIST(y1, q[m - 1]),  DIST(y1, q[m - 2]));
        d = min(d_right_weak,DIST(y0, q[m - 2]));
        dright+=d;

        if (dleft+dright + lb_vector_new [col] * lb_vector_new [col]>=threshold2){
            cnt++;
            lb_vector[col] = true;
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

            if (dleft+dright  + lb_vector_new [col]* lb_vector_new [col] >=threshold2){
                lb_vector[col] = true;
                cnt++;
                continue;
            }
            special_shared_vector[col][0]=dleft; special_shared_vector[col][1]=dright;
        }

    }

}

