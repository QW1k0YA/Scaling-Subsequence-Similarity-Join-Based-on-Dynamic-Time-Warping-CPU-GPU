
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
using namespace std;

void dist_mpx_v2O1_selectedV6_(const vector<DOUBLE>& a, int minlag, int subseqLen, int warpmax, DOUBLE bsf,
                               double adjust_factor, RETURN_V4& result)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    size_t subcount = a.size() - subseqLen + 1;

    vector<vector<bool>> del_test_profile(subcount,vector<bool>(subcount,false));
    vector<vector<bool>> lb_profile_dummy(subcount,vector<bool>(subcount,false));
    subcount = a.size() - subseqLen + 1;

    vector<DOUBLE> mu,sig;
    mu = movmean(a, 0, subseqLen - 1, 1);
    sig = movstd(a, 0, subseqLen - 1, 1);

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);

    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);

    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> MP(subcount,INFINITY);
    vector<DOUBLE> MPI(subcount,INFINITY);
    vector<DOUBLE> row_min_DEL(subcount,INFINITY);

    vector<DOUBLE> UTS_global, LTS_global;
    if (warpmax >= 0) {
        UTS_global = movmax(a, warpmax, warpmax);
        LTS_global = movmin(a, warpmax, warpmax);
    }

    vector<DOUBLE> cnt(a.size(), 0.0);

    vector<vector<DOUBLE>> subs(subseqLen, vector<DOUBLE>(subcount, 0.0));
    for (int i = 1; i <= subcount; i++) {
        for (int j = i; j <= i + subseqLen - 1; j++) {
            subs[j - i][i - 1] = (a[j - 1] - mu[i - 1]) / sig[i - 1];
        }
    }

    vector<DOUBLE> ttt(vector<DOUBLE>(static_cast<int>(a.size()) - 2 * subseqLen, subseqLen));
    vector<DOUBLE> aa(subseqLen);
    for(int i = 0; i < subseqLen; i++)
    {
        aa[i] = i+1;
    }

    vector<DOUBLE> bb = linspace(subseqLen, 1, subseqLen);

    vector<DOUBLE> mod = mergeVectors(aa,ttt,bb);

    int len_of_table = 600;
    vector<vector<DOUBLE>> count_table(subseqLen, vector<DOUBLE>(len_of_table, 0.0));

    int i;
    DOUBLE value;
    int pp;
    for(i = 1; i <= subseqLen; i++)
    {
        for(int j = 1;j <= subs[0].size();j++) {
            value = subs[i - 1][j - 1];
            pp = static_cast<int>(floor(value * 100 + len_of_table / 2));
            if (pp <= 0) {
                pp = 1;
            }
            if (pp > len_of_table) {
                pp = len_of_table;
            }
            count_table[i-1][pp-1] = count_table[i-1][pp-1] + 1;
        }
    }
    DOUBLE sum1;
    for(i = 1; i <= subseqLen; i++)
    {
        sum1 = 0;
        for(int j = 1;j <=len_of_table;j++)
        {
            sum1 += count_table[i-1][j-1];
            count_table[i-1][j-1] = sum1;
        }
    }

    vector<DOUBLE> normal_DIFF(a.size(), 0.0);
    DOUBLE UU,LL,posLL,posUU,diff2;
    for(i = 1;i <= subcount;i++)
    {
        for(int j =1; j <= subseqLen; j++)
        {
            LL = (LTS_global[i + j - 2] - mu[i - 1]) * invsig[i - 1];
            UU = (UTS_global[i + j - 2] - mu[i - 1]) * invsig[i - 1];

            posLL=floor(LL*100+ len_of_table/2);
            posUU=floor(UU*100+ len_of_table/2);

            if(posLL <= 0 )
            {
                posLL = 1;
            }
            if(posLL > len_of_table)
            {
                posLL = len_of_table;
            }
            if(posUU<=0)
            {
                posUU = 1;
            }
            if(posUU>len_of_table)
            {
                posUU = len_of_table;
            }

            diff2 = count_table[j-1][posUU-1] - count_table[j-1][posLL-1];
            normal_DIFF[i+j-2] += diff2 / a.size();
            cnt[i+j-2] += 1;
        }
    }

    normal_DIFF = elementWiseDivison_vv(normal_DIFF,mod);
    vector<DOUBLE> DIFF = normal_DIFF;

    DOUBLE avg_diff,std_diff,ths_diff;
    avg_diff = Mean(DIFF);
    std_diff = stddev(DIFF);
    ths_diff = avg_diff + adjust_factor*std_diff;

    DOUBLE qqq,www;
    if(adjust_factor == -100)
    {
        ths_diff = 0.5;
        qqq = avg_diff;
        www = std_diff;
    }

    i = 0;
    vector<DOUBLE> MASK(DIFF.size(),0.0);
    for(auto value:DIFF)
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

    UTS_global = elementWiseMultiply(UTS_global, MASK);
    LTS_global = elementWiseMultiply(LTS_global, MASK);

    vector<DOUBLE> dr_bwdU = addElementToFront(extr_vfromv(UTS_global, 1, subcount - 1), 0.0);
    vector<DOUBLE> dr_bwdL = addElementToFront(extr_vfromv(LTS_global, 1, subcount - 1), 0.0);

    vector<DOUBLE> dr_fwdU = extr_vfromv(UTS_global, subseqLen, UTS_global.size());
    vector<DOUBLE> dr_fwdL = extr_vfromv(LTS_global, subseqLen, LTS_global.size());

    vector<DOUBLE>dc_bwd = addElementToFront(extr_vfromv(a, 1, subcount - 1), 0.0);
    vector<DOUBLE>dc_fwd = extr_vfromv(a, subseqLen, a.size());

    vector<DOUBLE> dr_bwdU_plus_dr_bwdL, dr_fwdU_plus_dr_fwdL;
    dr_bwdU_plus_dr_bwdL = plusvector(dr_bwdU, dr_bwdL);
    dr_fwdU_plus_dr_fwdL = plusvector(dr_fwdU, dr_fwdL);

    vector<DOUBLE> sumU, sumL, sumU2, sumL2, sumU2_sumL2, sumU_sumL;
    sumU = movsum(UTS_global, subseqLen - 1);
    sumL = movsum(LTS_global, subseqLen - 1);
    sumU2 = movsum(elementWiseMultiply(UTS_global, UTS_global), subseqLen - 1);
    sumL2 = movsum(elementWiseMultiply(LTS_global, LTS_global), subseqLen - 1);
    sumU2_sumL2 = plusvector(sumU2, sumL2);
    sumU_sumL = plusvector(sumU, sumL);

    vector<DOUBLE> sumMASK = movsum(MASK, subseqLen - 1);
    vector<DOUBLE> norm_U_plus_norm_L(subcount,0.0);

    vector<DOUBLE> TS2 = elementWiseMultiply(a, a);

    for(int row = 1;row <= subcount;row++)
    {
        norm_U_plus_norm_L[row-1] = (sumU2_sumL2[row-1] - 2*sumU_sumL[row-1]*mu[row-1]
                                     + 2*sumMASK[row-1]*mu[row-1]*mu[row-1])*invsig[row-1]*invsig[row-1];
    }

    vector<DOUBLE> del(subcount,0);
    vector<DOUBLE> normLTS(subseqLen, 0.0), normUTS(subseqLen, 0.0);
    DOUBLE DUL2, DUL, del_ths;

    int debug_sum = 0;

    for(int row = 1;row <=subcount;row++)
    {
        vector<DOUBLE> LTS_temp(LTS_global);
        vector<DOUBLE> UTS_temp(UTS_global);

        for (int i = row; i <= row + subseqLen - 1; i++) {
            LTS_temp[i - 1] -= mu[row - 1];
            normLTS[i - row] = LTS_temp[i - 1] / sig[row - 1];
            UTS_temp[i - 1] -= mu[row - 1];
            normUTS[i - row] = UTS_temp[i - 1] / sig[row - 1];
        }

        DUL2 = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL = pow(MAX(DUL2, 0.0), 0.5);

        del_ths = ((2*bsf + DUL)*(2*bsf + DUL)+DUL*DUL)*0.5;
        del[row-1] = del_ths - norm_U_plus_norm_L[row-1];
    }
    vector<DOUBLE> dr_bwdMASK = addElementToFront(extr_vfromv(MASK, 1, subcount - 1), 0.0);
    vector<DOUBLE> dr_fwdMASK = extr_vfromv(MASK, subseqLen, MASK.size());

    vector<DOUBLE> dc_bwdTS2 = addElementToFront(extr_vfromv(TS2, 1, subcount - 1), 0.0);
    vector<DOUBLE> dc_fwdTS2 = extr_vfromv(TS2, subseqLen, TS2.size());
    DOUBLE dot_TS2_M;
    DOUBLE cov_U_plus_cov_L;
    DOUBLE dot_TS_M;
    DOUBLE M_NORM;
    DOUBLE t;
    DOUBLE dist;
    for(int diag = minlag+1;diag <= subcount;diag ++)
    {

        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(a.data() + diag - 1,
                                                           UTS_global.data(), LTS_global.data(), subseqLen));
        dot_TS_M = (elementWiseMultiply_p_sum(a.data() + diag - 1,
                                              MASK.data(), subseqLen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                               MASK.data(), subseqLen));
        int row = 1;
        int col = diag;

        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig[col - 1] * invsig[col - 1];

        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t;

        if(dist > del[row-1])
        {
            lb_profile_dummy[row-1][col-1] = true;
        }
        else if(dist < row_min_DEL[row-1])
        {
            row_min_DEL[row-1] = dist;
            MPI[row-1] = col;
        }
        
        for (row = 2; row <= subcount - diag + 1; row++) {
            col = diag + row - 1;
            cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL[row - 1] * dc_bwd[col - 1] +
                                                   (dr_fwdU_plus_dr_fwdL[row - 1]) * dc_fwd[col - 1];
            dot_TS_M = dot_TS_M - dr_bwdMASK[row - 1] * dc_bwd[col - 1] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
            dot_TS2_M = dot_TS2_M - dr_bwdMASK[row - 1] * dc_bwdTS2[col - 1] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
            M_NORM = M_NORM * invsig[col - 1]*invsig[col - 1];
            t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
                 2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
            dist = 2 * M_NORM - 2 * t ;

            if (dist > del[row - 1]) {
                lb_profile_dummy[row - 1][col - 1] = true;
            } else if (dist < row_min_DEL[row - 1]) {
                row_min_DEL[row - 1] = dist;
                MPI[row - 1] = col;
            }
        }
    }

    for(int diag = minlag+1;diag <= subcount;diag ++)
    {

        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(a.data() , UTS_global.data() + diag - 1,
                                                           LTS_global.data() + diag - 1, subseqLen));
        dot_TS_M = (elementWiseMultiply_p_sum(a.data(),
                                                    MASK.data() + diag - 1, subseqLen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                                     MASK.data() + diag - 1, subseqLen));
        int col = 1;
        int row = diag;

        M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
        M_NORM = M_NORM * invsig[col - 1] * invsig[col - 1];

        t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
             2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
        dist = 2 * M_NORM - 2 * t;

        if(dist > del[row-1])
        {
            lb_profile_dummy[row-1][col-1] = true;
        }
        else if(dist < row_min_DEL[row-1])
        {
            row_min_DEL[row-1] = dist;
            MPI[row-1] = col;
        }

        size_t temp_2 =  subcount - diag + 1;
        for (col = 2; col <= temp_2; col++) {
            row = diag + col - 1;
            cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL[row - 1] * dc_bwd[col - 1] +
                                                   (dr_fwdU_plus_dr_fwdL[row - 1]) * dc_fwd[col - 1];
            dot_TS_M = dot_TS_M - dr_bwdMASK[row - 1] * dc_bwd[col - 1] + dr_fwdMASK[row - 1] * dc_fwd[col - 1];
            dot_TS2_M = dot_TS2_M - dr_bwdMASK[row - 1] * dc_bwdTS2[col - 1] + dr_fwdMASK[row - 1] * dc_fwdTS2[col - 1];
            M_NORM = dot_TS2_M - 2 * mu[col - 1] * dot_TS_M + sumMASK[row - 1] * mu[col - 1] * mu[col - 1];
            M_NORM = M_NORM * invsig[col - 1]*invsig[col - 1];
            t = (cov_U_plus_cov_L - mu[col - 1] * sumU_sumL[row - 1] +
                 2 * mu[row - 1] * (sumMASK[row - 1] * mu[col - 1] - dot_TS_M)) * (invsig[col - 1] * invsig[row - 1]);
            dist = 2 * M_NORM - 2 * t ;

            if (dist > del[row - 1]) {
                lb_profile_dummy[col - 1][row - 1] = true;
            } else if (dist < row_min_DEL[row - 1]) {
                row_min_DEL[row - 1] = dist;
                MPI[row - 1] = col;
            }
        }
    }

    for (int row = 1; row <= subcount; row++) {
        for (int col = row; col <= MIN(size_t(row + minlag), subcount); col++) {
            lb_profile_dummy[row - 1][col - 1] = true;
        }
    }

    lb_profile_dummy = matrixOr(lb_profile_dummy, transposeMatrix_bool(lb_profile_dummy));

    for(int row = 1;row <=subcount;row++)
    {
        vector<DOUBLE> LTS_temp = LTS_global;
        vector<DOUBLE> UTS_temp = UTS_global;
        size_t temp_4 = row + subseqLen - 1;
        for ( i = row; i <= temp_4; i++) {
            LTS_temp[i - 1] -= mu[row - 1];
            normLTS[i - row] = LTS_temp[i - 1] / sig[row - 1];
            UTS_temp[i - 1] -= mu[row - 1];
            normUTS[i - row] = UTS_temp[i - 1] / sig[row - 1];
        }

        DUL2 = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL = pow(MAX(DUL2, 0.0), 0.5);

        row_min_DEL[row-1] += norm_U_plus_norm_L[row-1];
        row_min_DEL[row-1] = row_min_DEL[row-1]*2;

        MP[row-1] = (pow(row_min_DEL[row-1]- DUL2,0.5) - DUL)/2;
    }

    DOUBLE SSS2,aaa,Pratio;
    SSS2 = sum_bool_Matrix(lb_profile_dummy);
    aaa = subcount * subcount;
    Pratio = SSS2/aaa;
    auto end_time2 = std::chrono::high_resolution_clock::now();
    auto tmp_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time);

    printf("Factor:%5.2f   ths_diff: %5.2f, Selected V6 MPI purns :%5.6f , cost: %5.3fs \n",adjust_factor,ths_diff,Pratio,DOUBLE(tmp_time.count()/1000000));
    fflush(stdout);

    result.Pratio = Pratio;
    result.lb_profile_dummy = lb_profile_dummy;
    result.MPI = MPI;
    result.MP = MP;

}

