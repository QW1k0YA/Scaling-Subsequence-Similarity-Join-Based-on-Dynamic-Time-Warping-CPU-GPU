
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allunder/underDIAGV3.h"
#include "alldef/allstruct.h"

using namespace std;
void dist_mpx_DIAGV3(const vector<DOUBLE>& ts,int minlag,int subseqlen,
                                 int warpmax,DOUBLE bsf,const vector<vector<bool>>& lb_profile_previous,RETURN_FAST& result)
{
    int debug_sum = 0;
    auto local_fun = std::chrono::high_resolution_clock::now();
    size_t subcount = ts.size() - subseqlen + 1;

    vector<vector<bool>> lb_profile_dummy = lb_profile_previous;
    subcount = static_cast<int>(ts.size()) - subseqlen + 1;

    vector<DOUBLE> mu = movmean(ts, 0, subseqlen - 1, 1);
    vector<DOUBLE> sig = movstd(ts, 0, subseqlen - 1, 1);

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> MP(subcount, INFINITY);
    vector<DOUBLE> MPI(subcount, INFINITY);

    vector<vector<DOUBLE>> subs(subseqlen, vector<DOUBLE>(subcount, 0.0));

    for (int i = 1; i <= subcount; i++) {
        for (int j = i; j <= i + subseqlen - 1; j++) {
            subs[j - i][i - 1] = (ts[j - 1] - mu[i - 1]) / sig[i - 1];
        }
    }
    int Len = static_cast<int>(ts.size());

    RETURN_COM fromcom;

    for (int pos = 1; pos <= Len; pos++) {
        for (int j = 1; j <= fromcom.len_of_table; j++) {
            fromcom.count_table[pos - 1][j - 1] =
                    fromcom.count_table[pos - 1][j - 1] / fromcom.count_table[pos - 1][fromcom.len_of_table - 1];
        }
    }

    vector<DOUBLE> TS2 = elementWiseMultiply(ts, ts);
    
    const DOUBLE *dc_bwd_less = ts.data();
    const DOUBLE *dc_fwd = ts.data() + subseqlen - 1;
    const DOUBLE *dc_bwdTS2_less = TS2.data();
    const DOUBLE *dc_fwdTS2 = TS2.data() + subseqlen - 1;

    vector<DOUBLE> norm_U_plus_norm_L(subcount, 0);
    vector<DOUBLE> del(subcount, 0);

    int effective_max, effective_min;
    int pos;
    int posLL;
    int posUU;
    DOUBLE ths_diff;

    vector<DOUBLE>  sumU_sumL, sumU2_sumL2, sumMASK, dr_bwdMASK;
    double *dr_fwdU_plus_dr_fwdL;
    double  *dr_fwdMASK;
    vector<DOUBLE> raw_DIFF_UL, raw_DIFF_UL2, DUL2_raw, DUL2;
    DOUBLE DUL, cov_U_plus_cov_L, dot_TS_M, dot_TS2_M, M_NORM, t, dist;
    int diag;
    vector<DOUBLE> normal_DIFF;

    for (int diag_ID = subseqlen + 1; diag_ID <= subcount - 1; diag_ID++) {
        vector<DOUBLE> temp_normal_DIFF(Len,NAN);
        normal_DIFF = temp_normal_DIFF;
        effective_max = Len - (diag_ID - 1);

        for (int i = 1; i <= effective_max; i++) {
            pos = i + (diag_ID - 1);
            posLL = fromcom.pos_LL[i - 1];
            posUU = fromcom.pos_UU[i - 1];
            normal_DIFF[i - 1] = fromcom.count_table[pos - 1][posUU - 1] - fromcom.count_table[pos - 1][posLL - 1];
        }

        ths_diff = (meanWithoutNaN(normal_DIFF) + 0.5) / 2;

        vector<DOUBLE> MASK_local;
        for (DOUBLE value: normal_DIFF) {
            if (value <= ths_diff) {
                MASK_local.push_back(1);
            } else {
                MASK_local.push_back(0);
            }
        }

        effective_max = Len - (diag_ID-1);
        for(int i = 1;i <= effective_max;i++)
        {
            pos = i + (diag_ID-1);
            if(fromcom.attached_mask[pos-1])
            {

            }
        }

        vector<DOUBLE> UTS_local = elementWiseMultiply(fromcom.UTS, MASK_local);
        vector<DOUBLE> LTS_local = elementWiseMultiply(fromcom.LTS, MASK_local);
        vector<DOUBLE> UTS_plus_LTS_local = plusvector(UTS_local, LTS_local);

        const DOUBLE *dr_bwdU_plus_dr_bwdL_less = UTS_plus_LTS_local.data();
        dr_fwdU_plus_dr_fwdL = UTS_plus_LTS_local.data() + subseqlen - 1;

        sumU_sumL = movsum(UTS_plus_LTS_local, subseqlen - 1);
        sumU2_sumL2 = movsum(
                plusvector(elementWiseMultiply(UTS_local, UTS_local), elementWiseMultiply(LTS_local, LTS_local)),
                subseqlen - 1);

        sumMASK = movsum(MASK_local, subseqlen - 1);
        dr_bwdMASK = addElementToFront(extr_vfromv(MASK_local, 1, subcount - 1), 0);
        dr_fwdMASK = MASK_local.data() +  subseqlen - 1;

        vector<DOUBLE> temp_norm_U_plus_norm_L(subcount,0.0);
        norm_U_plus_norm_L = temp_norm_U_plus_norm_L;

        for (int row = 1; row <= subcount; row++) {
            norm_U_plus_norm_L[row - 1] = (sumU2_sumL2[row - 1] - 2 * sumU_sumL[row - 1] * mu[row - 1] +
                                           2 * sumMASK[row - 1] * mu[row - 1] * mu[row - 1]) * pow((invsig[row - 1]),2);
        }

        raw_DIFF_UL = substractvector(UTS_local, LTS_local);
        raw_DIFF_UL2 = elementWiseMultiply(raw_DIFF_UL, raw_DIFF_UL);
        DUL2_raw = movsum(raw_DIFF_UL2, subseqlen - 1);
        DUL2 = elementWiseDivison_vv(DUL2_raw, elementWiseMultiply(sig,sig));
        DOUBLE del_ths;
        for (int row = 1; row <= effective_max; row++) {
            DUL = sqrt(MAX(DUL2[row - 1], 0.0));
            del_ths = ((2*bsf + DUL)*(2*bsf + DUL) + DUL2[row-1])*0.5;
            del[row - 1] = del_ths;
        }

        diag = diag_ID;

        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() + diag - 1,UTS_local.data(),LTS_local.data(),subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data()+ diag - 1,
                                                    MASK_local.data(),subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                                     MASK_local.data(),subseqlen));
        int row = 1;
        int col = diag;

        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig[col - 1] * invsig[col - 1];
        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

        if (dist > del[row - 1]) {
            lb_profile_dummy[row - 1][col - 1] = true;
        } else if (dist < MP[row - 1]) {
            MP[row - 1] = dist;
            MPI[row - 1] = col;
        }

        for (row = 2; row <= subcount - diag + 1; row++) {
            col = diag + row - 1;

            if(row == 1 || col == 1)
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L + dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
            }
            else
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL_less[row - 2] * dc_bwd_less[col - 2] +
                                   dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
            }

            if(col == 1)
            {
                dot_TS_M = dot_TS_M + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }
            else
            {
                dot_TS_M = dot_TS_M - dr_bwdMASK[row - 1] * dc_bwd_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M - dr_bwdMASK[row - 1] * dc_bwdTS2_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }

            M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
            M_NORM = M_NORM * invsig[col - 1]* invsig[col - 1];
            t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
                 2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
            dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

            if (dist > del[row - 1]) {
                lb_profile_dummy[row - 1][col - 1] = true;
            } else if (dist < MP[row - 1]) {
                MP[row - 1] = dist;
                MPI[row - 1] = col;
            }
        }

        vector<DOUBLE> normal_DIFF_temp(Len, NAN);
        normal_DIFF = normal_DIFF_temp;

        effective_min = 1 + (diag_ID - 1);

        for (int i = effective_min; i <= Len; i++) {
            pos = i - (diag_ID - 1);
            posLL = fromcom.pos_LL[i - 1];
            posUU = fromcom.pos_UU[i - 1];
            normal_DIFF[i - 1] = fromcom.count_table[pos - 1][posUU - 1] - fromcom.count_table[pos - 1][posLL - 1];
        }

        ths_diff = (meanWithoutNaN(normal_DIFF) + 0.5) / 2;

        vector<DOUBLE> MASK_local_temp;
        MASK_local = MASK_local_temp;
        for (DOUBLE value: normal_DIFF) {
            if (value <= ths_diff) {
                MASK_local.push_back(1);
            } else {
                MASK_local.push_back(0);
            }
        }

        UTS_local = elementWiseMultiply(fromcom.UTS, MASK_local);
        LTS_local = elementWiseMultiply(fromcom.LTS, MASK_local);

        UTS_plus_LTS_local = plusvector(UTS_local, LTS_local);
        dr_bwdU_plus_dr_bwdL_less = UTS_plus_LTS_local.data();
        dr_fwdU_plus_dr_fwdL = UTS_plus_LTS_local.data() + subseqlen - 1;

        sumU_sumL = movsum(UTS_plus_LTS_local, subseqlen - 1);
        sumU2_sumL2 = movsum(
                plusvector(elementWiseMultiply(UTS_local, UTS_local), elementWiseMultiply(LTS_local, LTS_local)),
                subseqlen - 1);

        sumMASK = movsum(MASK_local, subseqlen - 1);
        dr_bwdMASK = addElementToFront(extr_vfromv(MASK_local, 1, subcount - 1), 0);
        dr_fwdMASK = MASK_local.data() +  subseqlen - 1;

        for (row = effective_min; row <= subcount; row++) {
            norm_U_plus_norm_L[row - 1] = (sumU2_sumL2[row - 1] - 2 * sumU_sumL[row - 1] * mu[row - 1] +
                                           2 * sumMASK[row - 1] * mu[row - 1] * mu[row - 1]) * (invsig[row - 1])* (invsig[row - 1]);
        }

        raw_DIFF_UL = substractvector(UTS_local, LTS_local);
        raw_DIFF_UL2 = elementWiseMultiply(raw_DIFF_UL, raw_DIFF_UL);

        DUL2_raw = movsum(raw_DIFF_UL2, subseqlen - 1);
        DUL2 = elementWiseMultiply(DUL2_raw, elementWiseMultiply(invsig ,invsig));

        for (row = effective_min; row <= subcount; row++) {
            DUL = sqrt(MAX(DUL2[row - 1], 0.0));
            del[row - 1] = ((2 * bsf + DUL) * (2 * bsf + DUL) + DUL2[row - 1]) * 0.5;
        }

        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() ,UTS_local.data()+ diag - 1,
                                                                 LTS_local.data()+ diag - 1,subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data(),
                                                    MASK_local.data() + diag - 1,subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                                     MASK_local.data() + diag - 1,subseqlen));
        col = 1;
        row = diag;

        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig[col - 1]* invsig[col - 1];
        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

        if (dist > del[row - 1]) {
            lb_profile_dummy[col - 1][row - 1] = true;

        } else if (dist < MP[row - 1]) {
            MP[row - 1] = dist;
            MPI[row - 1] = col;
        }

        for (col = 2; col <= subcount - diag + 1; col++) {
            row = diag + col - 1;

            M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
            M_NORM = M_NORM * invsig[col - 1]* invsig[col - 1];
            t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
                 2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
            dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];
            if(row == 1 || col == 1)
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L + dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
            }
            else
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL_less[row - 2] * dc_bwd_less[col - 2] +
                                   dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
            }

            if(col == 1)
            {
                dot_TS_M = dot_TS_M + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }
            else
            {
                dot_TS_M = dot_TS_M - dr_bwdMASK[row - 1] * dc_bwd_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M - dr_bwdMASK[row - 1] * dc_bwdTS2_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }
            if (dist > del[row - 1]) {
                lb_profile_dummy[col - 1][row - 1] = true;
            } else if (dist < MP[row - 1]) {
                MP[row - 1] = dist;
                MPI[row - 1] = col;
            }
        }

    }

    size_t temp_2;
    for (size_t row = 1; row <= subcount; row++) {
        temp_2 = MIN(row + minlag, subcount);
        for (size_t col = row; col <= temp_2; col++) {
            lb_profile_dummy[row - 1][col - 1] = true;
        }
    }

    lb_profile_dummy = matrixOr(lb_profile_dummy, transposeMatrix_bool(lb_profile_dummy));

    DOUBLE SSS2 = sum_bool_Matrix(lb_profile_dummy);
    DOUBLE P__ratio = SSS2 / (subcount * subcount);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - local_fun);

    printf("Local_MASK V3 MPI purns :%5.6f , cost: %5.3fs \n", P__ratio, DOUBLE(duration.count())/100000);

    result.lb_profile_dummy = lb_profile_dummy;
    result.MPI = MPI;
    result.MP = MP;
    result.P__ratio = P__ratio;
}
