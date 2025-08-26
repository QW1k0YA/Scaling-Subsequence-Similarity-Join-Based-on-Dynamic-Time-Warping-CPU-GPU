
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;

void dist_LB_in_one_Line_O1_V2_onlylb(const vector<DOUBLE>& ts, int minlag,
                                      int subseqlen,  DOUBLE bsf,int diagID,vector<vector<bool>>& lb_profile_dummy,
                                      const vector<vector<DOUBLE>> & count_table,const vector<DOUBLE>& pos_UU,const vector<DOUBLE>& pos_LL,
                                      const vector<DOUBLE>& mu,const vector<DOUBLE>& sig,int len_of_table,const vector<DOUBLE>& UTS_,const vector<DOUBLE>& LTS_)
{

    size_t len = ts.size();
    size_t subcount = len - subseqlen + 1;
    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    if(diagID < 0)
    {
        cerr << "diag_ID must >= 0 in distLB";
        exit(0);
    } else if (diagID <= minlag)
    {
        cerr << "ONLY IF diag_ID  > minlag , IT HAS MEANINGS  in distLB";
        exit(0);
    } else if(diagID > subcount)
    {
        cerr << "diag_ID CANNOT > subcount  in distLB" ;
        exit(0);
    }

    vector<DOUBLE> normal_DIFF(len,0.0);
    size_t posLL,posUU;
    size_t pos;
    DOUBLE diff2,diff_ratio;
    for(size_t i = 1;i <= len;i++)
    {
        pos = i + diagID -1;
        if(pos > len)
        {
            continue;
        }
        posLL = pos_LL[i-1];
        posUU = pos_UU[i-1];
        diff2 = count_table[pos-1][posUU-1] - count_table[pos-1][posLL-1];
        diff_ratio = diff2/count_table[pos-1][len_of_table-1];
        normal_DIFF[i-1] = diff_ratio;
    }
    DOUBLE avg_diff = meanWithoutNaN(normal_DIFF);
    DOUBLE ths_diff = (avg_diff + 0.5)/2;
    vector<DOUBLE> MASK(normal_DIFF.size());
    int i = 0;

    for(auto value:normal_DIFF)
    {
        if(value <= ths_diff)
        {
            MASK[i] = 1;
        }
        else
        {
            MASK[i] = 0;
        }
        i++;
    }

    vector<DOUBLE> MP(subcount,INFINITY);
    vector<DOUBLE> MPI(subcount,INFINITY);

    vector<DOUBLE> UTS = elementWiseMultiply(UTS_,MASK);
    vector<DOUBLE> LTS = elementWiseMultiply(LTS_,MASK);

    DOUBLE *dr_bwdU_plus_dr_bwdL_less, *dr_fwdU_plus_dr_fwdL;
    const double* dc_bwd_less = ts.data();
    
    const double *dc_fwd = ts.data() + subseqlen - 1;

    dr_bwdU_plus_dr_bwdL_less = plusvector_p(UTS.data(), LTS.data(),subcount-1);
    
    dr_fwdU_plus_dr_fwdL = plusvector_p(UTS.data() + subseqlen - 1, LTS.data() + subseqlen - 1,UTS.size() - subseqlen + 1);

    vector<DOUBLE>sumU = movsum(UTS, subseqlen - 1);
    vector<DOUBLE>sumL = movsum(LTS, subseqlen - 1);
    vector<DOUBLE>sumU2 = movsum(elementWiseMultiply(UTS, UTS), subseqlen - 1);
    vector<DOUBLE>sumL2 = movsum(elementWiseMultiply(LTS, LTS), subseqlen - 1);
    vector<DOUBLE>sumU2_sumL2 = plusvector(sumU2, sumL2);
    vector<DOUBLE>sumU_sumL = plusvector(sumU, sumL);

    vector<DOUBLE> sumMASK = movsum(MASK, subseqlen - 1);
    vector<DOUBLE> norm_U_plus_norm_L(subcount,0.0);

    for(size_t row = 1;row <= subcount;row++)
    {
        norm_U_plus_norm_L[row-1] = (sumU2_sumL2[row-1] - 2*sumU_sumL[row-1]*mu[row-1]
                                     + 2*sumMASK[row-1]*mu[row-1]*mu[row-1])*invsig[row-1]*invsig[row-1];
    }

    vector<DOUBLE> del(subcount,0.0);
    vector<DOUBLE> raw_DIFF_UL = substractvector(UTS,LTS);
    vector<DOUBLE> raw_DIFF_UL2 = elementWiseMultiply(raw_DIFF_UL,raw_DIFF_UL);
    vector<DOUBLE> DUL2_raw = movsum(raw_DIFF_UL2, subseqlen - 1);
    vector<DOUBLE> DUL2 = elementWiseDivison_vv(DUL2_raw, elementWiseMultiply(sig,sig));
    DOUBLE DUL,del_ths;

    for(size_t row = 1;row <= subcount;row ++)
    {
        DUL = sqrt(MAX(DUL2[row - 1], 0.0));
        del_ths = ((2*bsf+DUL)*(2*bsf+DUL)+ DUL2[row-1])*0.5;
        del[row-1] = del_ths;
    }
    vector<DOUBLE> TS2;
    elementWiseMultiply(ts, ts, TS2);
    DOUBLE *dr_bwdMASK_less = MASK.data();
    DOUBLE *dr_fwdMASK = MASK.data() + subseqlen - 1;
    DOUBLE *dc_bwdTS2_less = TS2.data();
    DOUBLE *dc_fwdTS2 = TS2.data() + subseqlen - 1;

    int diag = diagID;

    DOUBLE cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() + diag - 1,
                                                              UTS.data(),LTS.data(),subseqlen));
    DOUBLE dot_TS_M = (elementWiseMultiply_p_sum(ts.data()+ diag - 1,
                                                 MASK.data(),subseqlen));
    DOUBLE dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                                  MASK.data(),subseqlen));

    int row = 1;
    int col = diag;

    DOUBLE M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
    M_NORM = M_NORM * invsig[col - 1] * invsig[col - 1];

    DOUBLE t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
                2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
    DOUBLE dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

    if(dist > del[row-1])
    {
        lb_profile_dummy [row-1][col-1] = true;
    }
    else if(dist < MP[row-1])
    {
        MP[row-1] = dist;
        MPI[row-1] = col;
    }

    auto temp_1 = subcount - diag + 1;
    for (row = 2; row <= temp_1; row++) {
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
        M_NORM = M_NORM * invsig[col - 1]* invsig[col - 1];
        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

        if (dist > del[row - 1]) {
            lb_profile_dummy [row - 1][col - 1] = true;
        } else if (dist < MP[row - 1]) {
            MP[row - 1] = dist;
            MPI[row - 1] = col;
        }
    }

    cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() ,UTS.data()+ diag - 1,
                                                       LTS.data()+ diag - 1,subseqlen));
    dot_TS_M = (elementWiseMultiply_p_sum(ts.data(),
                                          MASK.data() + diag - 1,subseqlen));
    dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                           MASK.data() + diag - 1,subseqlen));
    col = 1;
    row = diag;

    M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
    M_NORM = M_NORM * invsig[col - 1]*invsig[col-1];
    t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
         2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
    dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

    if (dist > del[row - 1]) {
        lb_profile_dummy [col - 1][row - 1] = true;

    } else if (dist < MP[row - 1]) {
        MP[row - 1] = dist;
        MPI[row - 1] = col;
    }
    temp_1 = subcount - diag + 1;
    for (col = 2; col <= temp_1; col++) {
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
        M_NORM = M_NORM * invsig[col - 1]*invsig[col-1];
        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

        if (dist > del[row - 1]) {
            lb_profile_dummy [col - 1][row - 1] = true;
        } else if (dist < MP[row - 1]) {
            MP[row - 1] = dist;
            MPI[row - 1] = col;
        }
    }

}

void dist_LB_in_one_Line_O1_V2(const vector<DOUBLE>& ts, int minlag,
                                               int subseqlen,  DOUBLE bsf,int diagID,const vector<vector<bool>>& lb_profile_dummy,
                                               const vector<vector<DOUBLE>> & count_table,const vector<DOUBLE>& pos_UU,const vector<DOUBLE>& pos_LL,
                                               const vector<DOUBLE>& mu,const vector<DOUBLE>& sig,int len_of_table,const vector<DOUBLE>& UTS_,const vector<DOUBLE>& LTS_,LB_RETURN& result)
{
    result.lb_profile_dummy = lb_profile_dummy;
    size_t len = ts.size();
    size_t subcount = len - subseqlen + 1;
    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    if(diagID < 0)
    {
        cerr << "diag_ID must >= 0 in distLB";
        exit(0);
    } else if (diagID <= minlag)
    {
        cerr << "ONLY IF diag_ID  > minlag , IT HAS MEANINGS  in distLB";
        exit(0);
    } else if(diagID > subcount)
    {
        cerr << "diag_ID CANNOT > subcount  in distLB" ;
        exit(0);
    }

    vector<DOUBLE> normal_DIFF(len,0.0);
    size_t posLL,posUU;
    size_t pos;
    DOUBLE diff2,diff_ratio;
    for(size_t i = 1;i <= len;i++)
    {
        pos = i + diagID -1;
        if(pos > len)
        {
            continue;
        }
        posLL = pos_LL[i-1];
        posUU = pos_UU[i-1];
        diff2 = count_table[pos-1][posUU-1] - count_table[pos-1][posLL-1];
        diff_ratio = diff2/count_table[pos-1][len_of_table-1];
        normal_DIFF[i-1] = diff_ratio;
    }
    DOUBLE avg_diff = meanWithoutNaN(normal_DIFF);
    DOUBLE ths_diff = (avg_diff + 0.5)/2;
    vector<DOUBLE> MASK(normal_DIFF.size());
    int i = 0;

    for(auto value:normal_DIFF)
    {
        if(value <= ths_diff)
        {
            MASK[i] = 1;
        }
        else
        {
            MASK[i] = 0;
        }
        i++;
    }

    vector<DOUBLE> MP(subcount,INFINITY);
    vector<DOUBLE> MPI(subcount,INFINITY);

    vector<DOUBLE> UTS = elementWiseMultiply(UTS_,MASK);
    vector<DOUBLE> LTS = elementWiseMultiply(LTS_,MASK);

    vector<DOUBLE>dc_bwd = addElementToFront(extr_vfromv(ts, 1, subcount - 1), 0.0);
    vector<DOUBLE>dc_fwd = extr_vfromv(ts, subseqlen, ts.size());

    DOUBLE *dr_bwdU_plus_dr_bwdL_less, *dr_fwdU_plus_dr_fwdL;

    dr_bwdU_plus_dr_bwdL_less = plusvector_p(UTS.data(), LTS.data(),subcount-1);
    
    dr_fwdU_plus_dr_fwdL = plusvector_p(UTS.data() + subseqlen - 1, LTS.data() + subseqlen - 1,UTS.size() - subseqlen + 1);

    vector<DOUBLE> sumU = movsum(UTS, subseqlen - 1);
    vector<DOUBLE> sumL = movsum(LTS, subseqlen - 1);
    vector<DOUBLE> sumU2 = movsum(elementWiseMultiply(UTS, UTS), subseqlen - 1);
    vector<DOUBLE> sumL2 = movsum(elementWiseMultiply(LTS, LTS), subseqlen - 1);
    vector<DOUBLE> sumU2_sumL2 = plusvector(sumU2, sumL2);
    vector<DOUBLE> sumU_sumL = plusvector(sumU, sumL);

    vector<DOUBLE> sumMASK = movsum(MASK, subseqlen - 1);
    vector<DOUBLE> norm_U_plus_norm_L(subcount,0.0);

    for(size_t row = 1;row <= subcount;row++)
    {
        norm_U_plus_norm_L[row-1] = (sumU2_sumL2[row-1] - 2*sumU_sumL[row-1]*mu[row-1]
                + 2*sumMASK[row-1]*mu[row-1]*mu[row-1])*invsig[row-1]*invsig[row-1];
    }

    vector<DOUBLE> del(subcount,0.0);
    vector<DOUBLE> raw_DIFF_UL = substractvector(UTS,LTS);
    vector<DOUBLE> raw_DIFF_UL2 = elementWiseMultiply(raw_DIFF_UL,raw_DIFF_UL);
    vector<DOUBLE> DUL2_raw = movsum(raw_DIFF_UL2, subseqlen - 1);
    vector<DOUBLE> DUL2 = elementWiseDivison_vv(DUL2_raw, elementWiseMultiply(sig,sig));
    DOUBLE DUL,del_ths;

    for(size_t row = 1;row <= subcount;row ++)
    {
        DUL = sqrt(MAX(DUL2[row - 1], 0.0));
        del_ths = ((2*bsf+DUL)*(2*bsf+DUL)+ DUL2[row-1])*0.5;
        del[row-1] = del_ths;
    }

    vector<DOUBLE> TS2 = elementWiseMultiply(ts,ts);
    vector<DOUBLE> dr_bwdMASK = addElementToFront(extr_vfromv(MASK, 1, subcount - 1), 0.0);
    vector<DOUBLE> dr_fwdMASK = extr_vfromv(MASK, subseqlen, MASK.size());

    vector<DOUBLE> dc_bwdTS2 = addElementToFront(extr_vfromv(TS2, 1, subcount - 1), 0.0);
    vector<DOUBLE> dc_fwdTS2 = extr_vfromv(TS2, subseqlen, TS2.size());

    int diag = diagID;

    DOUBLE cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() + diag - 1,
                                                             UTS.data(),LTS.data(),subseqlen));
    DOUBLE dot_TS_M = (elementWiseMultiply_p_sum(ts.data()+ diag - 1,
                                                MASK.data(),subseqlen));
    DOUBLE dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                                 MASK.data(),subseqlen));

    int row = 1;
    int col = diag;

    DOUBLE M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
    M_NORM = M_NORM * invsig[col - 1] * invsig[col - 1];

    DOUBLE t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
         2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
    DOUBLE dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

    if(dist > del[row-1])
    {
        result.lb_profile_dummy [row-1][col-1] = true;
    }
    else if(dist < MP[row-1])
    {
        MP[row-1] = dist;
        MPI[row-1] = col;
    }

    auto temp_1 = subcount - diag + 1;
    for (row = 2; row <= temp_1; row++) {
        col = diag + row - 1;
        if(row == 1)
        {
            cov_U_plus_cov_L = cov_U_plus_cov_L + (dr_fwdU_plus_dr_fwdL[row - 1]) * dc_fwd[col - 1];
        }
        else
        {
            cov_U_plus_cov_L = cov_U_plus_cov_L - (dr_bwdU_plus_dr_bwdL_less[row - 2] * dc_bwd[col - 1] +
                                                   (dr_fwdU_plus_dr_fwdL[row - 1]) * dc_fwd[col - 1]);
        }

        dot_TS_M = dot_TS_M - dr_bwdMASK[row - 1] * dc_bwd[col - 1] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
        dot_TS2_M = dot_TS2_M - dr_bwdMASK[row - 1] * dc_bwdTS2[col - 1] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig[col - 1];
        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

        if (dist > del[row - 1]) {
            result.lb_profile_dummy [row - 1][col - 1] = true;
        } else if (dist < MP[row - 1]) {
            MP[row - 1] = dist;
            MPI[row - 1] = col;
        }
    }

    cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() ,UTS.data()+ diag - 1,
                                                             LTS.data()+ diag - 1,subseqlen));
    dot_TS_M = (elementWiseMultiply_p_sum(ts.data(),
                                                MASK.data() + diag - 1,subseqlen));
    dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                                 MASK.data() + diag - 1,subseqlen));
    col = 1;
    row = diag;

    M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
    M_NORM = M_NORM * invsig[col - 1]*invsig[col-1];
    t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
         2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
    dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

    if (dist > del[row - 1]) {
        result.lb_profile_dummy [col - 1][row - 1] = true;

    } else if (dist < MP[row - 1]) {
        MP[row - 1] = dist;
        MPI[row - 1] = col;
    }
    temp_1 = subcount - diag + 1;
    for (col = 2; col <= temp_1; col++) {
        row = diag + col - 1;
        if(row == 1)
        {
            cov_U_plus_cov_L = cov_U_plus_cov_L +
                               (dr_fwdU_plus_dr_fwdL[row - 1]) * dc_fwd[col - 1];
        }
        else
        {
            cov_U_plus_cov_L = cov_U_plus_cov_L - (dr_bwdU_plus_dr_bwdL_less[row - 2]) * dc_bwd[col - 1] +
                               (dr_fwdU_plus_dr_fwdL[row - 1]) * dc_fwd[col - 1];
        }

        dot_TS_M = dot_TS_M - dr_bwdMASK[row-1]* dc_bwd[col-1] + dr_fwdMASK[row-1]* dc_fwd[col-1];

        dot_TS2_M = dot_TS2_M - dr_bwdMASK[row-1]* dc_bwdTS2[col-1] + dr_fwdMASK[row-1]* dc_fwdTS2[col-1];

        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig[col - 1]*invsig[col-1];
        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row - 1];

        if (dist > del[row - 1]) {
            result.lb_profile_dummy [col - 1][row - 1] = true;
        } else if (dist < MP[row - 1]) {
            MP[row - 1] = dist;
            MPI[row - 1] = col;
        }
    }

    result.MP = MP;
    result.MPI = MPI;
    result.Pratio = NAN;

}

