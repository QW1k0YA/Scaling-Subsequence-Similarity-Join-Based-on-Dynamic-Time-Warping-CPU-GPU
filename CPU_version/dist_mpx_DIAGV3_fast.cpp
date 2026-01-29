
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allunder/underfast.h"
#include "alldef/allstruct.h"

using namespace std;

void dist_mpx_DIAGV3_fast(const vector<DOUBLE>& ts,size_t minlag,size_t subseqlen,
                          size_t warpmax,DOUBLE bsf,const vector<vector<bool>>& lb_profile_previous,    RETURN_FAST& result)
{
    int debug_sum = 0;
    auto local_fun = std::chrono::high_resolution_clock::now();
    auto subcount = ts.size() - subseqlen + 1;

    vector<vector<bool>> lb_profile_dummy = lb_profile_previous;
    subcount = ts.size() - subseqlen + 1;

    vector<DOUBLE> mu = movmean(ts, 0, subseqlen - 1, 1);
    vector<DOUBLE> sig = movstd(ts, 0, subseqlen - 1, 1);

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<vector<DOUBLE>> sig_row2 = convertTo2DRowOrder(elementWiseMultiply(sig, sig));
    vector<vector<DOUBLE>> invsig_m2 = elementWiseDivision(1.0, sig_row2);
    vector<DOUBLE> invsig2 = extractVector(invsig_m2, 0, 1);

    vector<DOUBLE> MP(subcount, INFINITY);
    vector<DOUBLE> MPI(subcount, INFINITY);

    vector<vector<DOUBLE>> subs(subseqlen, vector<DOUBLE>(subcount, 0.0));

    size_t temp_1;
    for (size_t i = 1; i <= subcount; i++) {
        temp_1 =  i + subseqlen - 1;
        for (size_t j = i; j <= temp_1; j++) {
            subs[j - i][i - 1] = (ts[j - 1] - mu[i - 1]) / sig[i - 1];
        }
    }
    size_t Len = ts.size();

    RETURN_COM fromcom;
    compute_shared_dataV3(ts, subseqlen, warpmax, mu, sig, subs,fromcom);

    for (size_t pos = 1; pos <= Len; pos++) {
        for (size_t j = 1; j <= fromcom.len_of_table; j++) {
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

    size_t effective_max_row, effective_min_row;
    size_t pos;
    size_t posLL;
    size_t posUU;
    DOUBLE ths_diff;

    vector<DOUBLE>  UTS_local,LTS_local,UTS_plus_LTS_local,sumU_sumL, sumU2_sumL2, sumMASK;
    double *dr_fwdU_plus_dr_fwdL;
    double *dr_bwdU_plus_dr_bwdL_less;
    double  *dr_fwdMASK, *dr_bwdMASK_less;

    DOUBLE DUL, cov_U_plus_cov_L, dot_TS_M, dot_TS2_M, M_NORM, t, dist;
    size_t diag;
    vector<DOUBLE> normal_DIFF;

    for (size_t diag_ID = subseqlen + 1; diag_ID <= subcount - 1; diag_ID++) {
        vector<DOUBLE> temp_normal_DIFF(Len, NAN);
        normal_DIFF = temp_normal_DIFF;
        effective_max_row = Len - (diag_ID - 1);

        for (size_t i = 1; i <= effective_max_row; i++) {
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

        elementWiseMultiply(fromcom.UTS, MASK_local, UTS_local);
        elementWiseMultiply(fromcom.LTS, MASK_local, LTS_local);

        plusvector(UTS_local, LTS_local, UTS_plus_LTS_local);

        dr_bwdU_plus_dr_bwdL_less = UTS_plus_LTS_local.data();
        dr_fwdU_plus_dr_fwdL = UTS_plus_LTS_local.data() + subseqlen - 1;

        sumU_sumL = movsum(UTS_plus_LTS_local, subseqlen - 1);
        sumU2_sumL2 = movsum(
                plusvector(elementWiseMultiply(UTS_local, UTS_local), elementWiseMultiply(LTS_local, LTS_local)),
                subseqlen - 1);

        sumMASK = movsum(MASK_local, subseqlen - 1);
        dr_bwdMASK_less = MASK_local.data();
        dr_fwdMASK = MASK_local.data() +  subseqlen - 1;

        for (size_t row = 1; row <= effective_max_row; row++) {
            norm_U_plus_norm_L[row - 1] = (sumU2_sumL2[row - 1] - 2 * sumU_sumL[row - 1] * mu[row - 1] +
                                           2 * sumMASK[row - 1] * mu[row - 1] * mu[row - 1]) * (invsig2[row - 1]);
        }

        vector<DOUBLE> raw_DIFF_UL = substractvector(UTS_local, LTS_local);
        vector<DOUBLE>raw_DIFF_UL2 = elementWiseMultiply(raw_DIFF_UL, raw_DIFF_UL);
        vector<DOUBLE>DUL2_raw = movsum(raw_DIFF_UL2, subseqlen - 1);
        vector<DOUBLE>DUL2 = elementWiseMultiply(DUL2_raw, invsig2);

        for (size_t row = 1; row <= effective_max_row; row++) {
            DUL = sqrt(DUL2[row - 1]);
            del[row - 1] = ((2 * bsf + DUL) * (2 * bsf + DUL) + DUL2[row - 1]) * 0.5;
        }

        diag = diag_ID;

        cov_U_plus_cov_L = (elementWiseMultiply_p_sum(ts.data() + diag - 1,
                                                            plusvector_p(UTS_local.data(),
                                                                         LTS_local.data(),subseqlen),subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data()+ diag - 1,
                                                    MASK_local.data(),subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                                     MASK_local.data(),subseqlen));

        size_t row = 1;
        size_t col = diag;

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
                dot_TS_M = dot_TS_M + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }
            else
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL_less[row - 2] * dc_bwd_less[col - 2] +
                                   dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
                dot_TS_M = dot_TS_M - dr_bwdMASK_less[row - 2] * dc_bwd_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M - dr_bwdMASK_less[row - 2] * dc_bwdTS2_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }

            M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
            M_NORM = M_NORM * invsig2[col - 1];
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

        effective_min_row = 1 + (diag_ID - 1);

        for (size_t i = effective_min_row; i <= Len; i++) {
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

        elementWiseMultiply(fromcom.UTS, MASK_local, UTS_local);
        elementWiseMultiply(fromcom.LTS, MASK_local, LTS_local);

        plusvector(UTS_local, LTS_local, UTS_plus_LTS_local);

        dr_bwdU_plus_dr_bwdL_less = UTS_plus_LTS_local.data();
        dr_fwdU_plus_dr_fwdL = UTS_plus_LTS_local.data() + subseqlen - 1;

        sumU_sumL = movsum(UTS_plus_LTS_local, subseqlen - 1);
        sumU2_sumL2 = movsum(
                plusvector(elementWiseMultiply(UTS_local, UTS_local), elementWiseMultiply(LTS_local, LTS_local)),
                subseqlen - 1);

        sumMASK = movsum(MASK_local, subseqlen - 1);
        dr_bwdMASK_less = MASK_local.data();
        dr_fwdMASK = MASK_local.data() +  subseqlen - 1;
        for (row = effective_min_row; row <= subcount; row++) {
            norm_U_plus_norm_L[row - 1] = (sumU2_sumL2[row - 1] - 2 * sumU_sumL[row - 1] * mu[row - 1] +
                                           2 * sumMASK[row - 1] * mu[row - 1] * mu[row - 1]) * (invsig2[row - 1]);
        }

        raw_DIFF_UL = substractvector(UTS_local, LTS_local);
        raw_DIFF_UL2 = elementWiseMultiply(raw_DIFF_UL, raw_DIFF_UL);

        DUL2_raw = movsum(raw_DIFF_UL2, subseqlen - 1);
        DUL2 = elementWiseMultiply(DUL2_raw, invsig2);

        for (row = effective_min_row; row <= subcount; row++) {
            DUL = sqrt(DUL2[row - 1]);
            del[row - 1] = ((2 * bsf + DUL) * (2 * bsf + DUL) + DUL2[row - 1]) * 0.5;
        }

        cov_U_plus_cov_L = (elementWiseMultiply_p_sum(ts.data(),
                                                            plusvector_p(UTS_local.data() + diag - 1,
                                                                         LTS_local.data() + diag - 1,subseqlen),subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data(),
                                                    MASK_local.data() + diag - 1,subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                                     MASK_local.data() + diag - 1,subseqlen));

        col = 1;
        row = diag;

        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig2[col - 1];
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

            if(row == 1 || col == 1)
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L + dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
                dot_TS_M = dot_TS_M + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }
            else
            {
                cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL_less[row - 2] * dc_bwd_less[col - 2] +
                                   dr_fwdU_plus_dr_fwdL[row - 1] * dc_fwd[col - 1];
                dot_TS_M = dot_TS_M - dr_bwdMASK_less[row - 2] * dc_bwd_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
                dot_TS2_M = dot_TS2_M - dr_bwdMASK_less[row - 2] * dc_bwdTS2_less[col - 2] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            }

            M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
            M_NORM = M_NORM * invsig2[col - 1];
            t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
                 2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
            dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

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

    printf("Local_MASK V3 fast version MPI purns :%5.6f  cost: %5.3f s \n", P__ratio, DOUBLE(duration.count()/1000000));

    result.lb_profile_dummy = lb_profile_dummy;
    result.MPI = MPI;
    result.MP = MP;
    result.P__ratio = P__ratio;

}

