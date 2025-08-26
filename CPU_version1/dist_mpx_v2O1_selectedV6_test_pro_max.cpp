
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/typedefdouble.h"
#include "allunder/underpromax.h"
#include "alldef/allstruct.h"
using namespace std;

void dist_mpx_v2O1_selectedV6_test_pro_max(const vector<DOUBLE>& ts, int minlag,
                                                    int subseqlen, int  warpmax, DOUBLE bsf,DOUBLE adjust_factor,const vector<vector<bool>>& lb_profile_previous,PROMAX_RETURN& result)
{
    int debug_sum =  0;
    auto start_time = std::chrono::high_resolution_clock::now();

    int subcount = ts.size() - subseqlen + 1;

    vector<vector<DOUBLE>> del_test_profile(subcount,vector<DOUBLE>(subcount,0.0));
    vector<vector<bool>> lb_profile_dummy=falseMatrix(subcount,subcount);

    subcount = ts.size() - subseqlen + 1;

    vector<DOUBLE> mu = movmean(ts,0,subseqlen - 1,1);
    vector<DOUBLE> sig = movstd(ts, 0, subseqlen - 1, 1);
    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> MP(subcount, INFINITY);
    vector<DOUBLE> MPI(subcount, INFINITY);
    vector<DOUBLE> row_min_DEL(subcount, INFINITY);

    vector<DOUBLE>     UTS = movmax(ts, warpmax, warpmax);
    vector<DOUBLE>    LTS = movmin(ts, warpmax, warpmax);

    vector<vector<DOUBLE>> subs(subseqlen, vector<DOUBLE>(subcount, 0.0));
    size_t temp_1;
    for (int i = 1; i <= subcount; i++) {
        temp_1 = i + subseqlen - 1;
        for (int j = i; j <= temp_1; j++) {
            subs[j - i][i - 1] = (ts[j - 1] - mu[i - 1]) / sig[i - 1];
        }
    }

    vector<DOUBLE> ttt(vector<DOUBLE>(static_cast<int>(ts.size())-2*subseqlen,subseqlen));

    vector<DOUBLE> aa(subseqlen);
    for(int i = 0;i < subseqlen; i++)
    {
        aa[i] = i+1;
    }

    vector<DOUBLE> bb = linspace(subseqlen,1,subseqlen);
    vector<DOUBLE> mod = mergeVectors(aa,ttt,bb);

    DOUBLE SSS2 = sum_bool_Matrix(lb_profile_previous);
    int aaa = subcount * subcount;
    DOUBLE P__ratio = SSS2/aaa;

    vector<DOUBLE> sum_in_row(lb_profile_previous.size());
    int i = 0;
    for(auto value:lb_profile_previous)
    {

        sum_in_row[i] = sum_vector_bool(value);
        i++;
    }

    vector<bool> effctive_vector(sum_in_row.size());
    i =  0;
    for(auto value:sum_in_row)
    {
        if(sum_in_row[i] < subcount*P__ratio)
        {
            effctive_vector[i] = 1;
        }
        else
        {
            effctive_vector[i] = 0;
        }
        i++;
    }

    vector<DOUBLE> my_diag_sum(subcount,0);

    for(int diagID = 1;diagID <= subcount;diagID++)
    {
        my_diag_sum[diagID-1] = my_diag_sum[diagID-1] + diagID;
        for(int row = 1;row <= subcount;row ++)
        {
            int col = diagID + row - 1;
            if(col > subcount)
            {
                break;
            }
            if(lb_profile_previous[row-1][col-1])
            {
                my_diag_sum[diagID-1] = my_diag_sum[diagID-1] + 1;
            }
        }

    }

    start_time = std::chrono::high_resolution_clock::now();
    
    int len_of_table=600;
    RETURN_COM result1;

    vector<vector<SHORT>> count_table = result1.count_table;
    vector<DOUBLE> pos_UU = result1.pos_UU;
    vector<DOUBLE> pos_LL = result1.pos_LL;

    vector<bool> lb_vector(subcount);
    for(int diag_ID = subseqlen + 1;diag_ID <= subcount-1;diag_ID++)
    {

    }
    
    lb_profile_dummy = matrixOr(lb_profile_dummy, transposeMatrix_bool(lb_profile_dummy));

    debug_sum = sum_bool_Matrix(lb_profile_dummy);
    cout << debug_sum << "promax in 139" << endl;

    SSS2 = sum_bool_Matrix(lb_profile_dummy);
    aaa = subcount*subcount;
    DOUBLE Pratio = SSS2/ aaa;

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    vector<vector<bool>> lb_profile_dummy___ = lb_profile_dummy;
    vector<DOUBLE>my_diag_sum_(subcount,0.0);
    my_diag_sum = my_diag_sum_;

    for(int diagID = 1;diagID <= subcount;diagID++)
    {
        my_diag_sum[diagID-1] = my_diag_sum[diagID-1] + diagID;
        for(int row = 1;row <= subcount;row++)
        {
            int col = diagID + row - 1;
            if(col > subcount)
            {
                break;
            }
            if(lb_profile_dummy___[row-1][col-1])
            {
                my_diag_sum[diagID-1] = my_diag_sum[diagID-1]+1;
            }
        }
    }

    DOUBLE temp1 = 0;
    for(auto value:sum_in_row)
    {
        if(value < (subcount*P__ratio))
        {
            temp1++;
        }
    }
    printf("Number of exempted sequences : %f\n", 1-temp1/subcount);

    vector<DOUBLE> cnt(ts.size(),0.0);
    len_of_table = 600;
    vector<vector<SHORT>> count_table_(subseqlen,vector<SHORT>(len_of_table,0.0));
    count_table = count_table_;
    DOUBLE value,pp;

    for(i = 1;i <= subseqlen; i++)
    {
        for(int j = 1;j <=subcount;j++)
        {
            if(!effctive_vector[j-1])
            {
                continue;
            }
            value = subs[i-1][j-1];
            pp = floor(value*100 + len_of_table/2.0);
            if(pp <= 0)
            {
                pp = 1;
            }
            if(pp >= len_of_table)
            {
                pp = len_of_table;
            }

            count_table[i-1][static_cast<int>(pp)-1] = count_table[i-1][static_cast<int>(pp)-1] + 1;
        }
    }

    DOUBLE sum1 = 0;
    for(i = 1;i <= subseqlen;i++)
    {
        sum1 = 0;
        for(int j = 1;j <=len_of_table;j++)
        {
            sum1+=count_table[i-1][j-1];
            count_table[i-1][j-1] = sum1;
        }
    }

    vector<DOUBLE> normal_DIFF(ts.size(),0);
    DOUBLE LL,UU,posLL,posUU,diff2,diff_ratio,avg_diff,std_diff,ths_diff;
    vector<DOUBLE> DIFF(ts.size());
    vector<DOUBLE> MASK(DIFF.size());
    vector<bool> BACKUPMASK(UTS.size());
    vector<DOUBLE> BACKUP_DIFF(UTS.size());
    DOUBLE BACKUP_avg_diff,BACKUP_ths_diff,BACKUP_std_diff;

    for(i = 1;i <= subcount;i++) {
        if (!effctive_vector[i - 1]) {
            continue;
        }
        for (int j = 1; j <= subseqlen; j++) {
            LL = (LTS[i + j - 2] - mu[i - 1]) * invsig[i - 1];
            UU = (UTS[i + j - 2] - mu[i - 1]) * invsig[i - 1];

            posLL = floor(LL * 100 + len_of_table / 2);
            posUU = floor(UU * 100 + len_of_table / 2);
            if (posLL <= 0) {
                posLL =1 ;
            }
            if (posLL >= len_of_table) {
                posLL = len_of_table;
            }
            if (posUU <= 0) {
                posUU = 1;
            }
            if (posUU >= len_of_table) {
                posUU = len_of_table;
            }
            diff2 = count_table[j - 1][posUU - 1] - count_table[j - 1][posLL - 1];
            diff_ratio = diff2 / count_table[j - 1][len_of_table - 1];

            normal_DIFF[i + j - 2] += diff_ratio;
            cnt[i + j - 2] += 1;
        }
    }

    elementWiseDivison_vv(normal_DIFF, cnt, normal_DIFF);
    DIFF = normal_DIFF;

    avg_diff = meanWithoutNaN(DIFF);
    std_diff = stdWithOmitnan(DIFF);
    ths_diff = avg_diff + adjust_factor*std_diff;

    DOUBLE qqq,www;
    if(adjust_factor == -100)
    {
        ths_diff = 0.5;
        qqq = avg_diff;
        www = std_diff;
    }
    ths_diff = (avg_diff +0.5)/2;

    i = 0;
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
    BACKUP_DIFF = substractvector(UTS,LTS);
    BACKUP_avg_diff = Mean(BACKUP_DIFF);
    BACKUP_std_diff = stddev(BACKUP_DIFF);
    BACKUP_ths_diff = BACKUP_avg_diff + 0*BACKUP_std_diff;

    i = 0;
    for(auto value:BACKUP_DIFF)
    {
        if(value <= BACKUP_ths_diff)
        {
            BACKUPMASK[i] = 1;
        }
        else
        {
            BACKUPMASK[i] = 0;
        }
        i++;
    }

    auto normal_DIFF_size = normal_DIFF.size();
    for(i = 1;i <= normal_DIFF_size;i++)
    {
        if(isnan(normal_DIFF[i-1]))
        {
            MASK[i-1] = BACKUPMASK[i-1];
        }
    }
    elementWiseMultiply(UTS, MASK, UTS);
    elementWiseMultiply(LTS, MASK, LTS);

    DOUBLE *dr_bwdU_plus_dr_bwdL_less, *dr_fwdU_plus_dr_fwdL;
    const double* dc_bwd_less = ts.data();
    
    const double *dc_fwd = ts.data() + subseqlen - 1;

    dr_bwdU_plus_dr_bwdL_less = plusvector_p(UTS.data(), LTS.data(),subcount-1);
    
    dr_fwdU_plus_dr_fwdL = plusvector_p(UTS.data() + subseqlen - 1, LTS.data() + subseqlen - 1,UTS.size() - subseqlen + 1);

    vector<DOUBLE> sumU = movsum(UTS, subseqlen - 1);
    vector<DOUBLE> sumL = movsum(LTS, subseqlen - 1);
    vector<DOUBLE> sumU2 = movsum(elementWiseMultiply(UTS, UTS), subseqlen - 1);
    vector<DOUBLE>sumL2 = movsum(elementWiseMultiply(LTS, LTS), subseqlen - 1);
    vector<DOUBLE>sumU2_sumL2 = plusvector(sumU2, sumL2);
    vector<DOUBLE>sumU_sumL = plusvector(sumU, sumL);

    vector<DOUBLE> sumMASK = movsum(MASK, subseqlen - 1);
    vector<DOUBLE> norm_U_plus_norm_L(subcount,0.0);

    vector<DOUBLE> TS2 = elementWiseMultiply(ts,ts);

    for(int row = 1;row <= subcount;row++)
    {
        norm_U_plus_norm_L[row-1] = (sumU2_sumL2[row-1] - 2*sumU_sumL[row-1]*mu[row-1]
                                     + 2*sumMASK[row-1]*mu[row-1]*mu[row-1])*invsig[row-1]*invsig[row-1];
    }

    vector<DOUBLE> del(subcount,0);
    vector<DOUBLE> normLTS(subseqlen, 0.0), normUTS(subseqlen, 0.0);
    DOUBLE DUL2, DUL, del_ths;

    for(int row = 1;row <=subcount;row++)
    {

        for (int i = row; i <= row + subseqlen - 1; i++) {
            normLTS[i - row] = (LTS[i - 1] - mu[row - 1]) / sig[row - 1];
            normUTS[i - row] = (UTS[i - 1] - mu[row - 1]) / sig[row - 1];
        }

        DUL2 = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL = pow(MAX(DUL2, 0.0), 0.5);

        del_ths = ((2*bsf + DUL)*(2*bsf + DUL)+DUL*DUL)*0.5;
        del[row-1] = del_ths - norm_U_plus_norm_L[row-1];
    }

    DOUBLE *dr_bwdMASK_less = MASK.data();
    DOUBLE *dr_fwdMASK = MASK.data() + subseqlen - 1;
    DOUBLE *dc_bwdTS2_less = TS2.data();
    DOUBLE *dc_fwdTS2 = TS2.data() + subseqlen - 1;

    DOUBLE dot_TS2_M;
    DOUBLE cov_U_plus_cov_L;
    DOUBLE dot_TS_M;
    DOUBLE M_NORM;
    DOUBLE t;
    DOUBLE dist;
    debug_sum = sum_bool_Matrix(lb_profile_dummy);
    cout << debug_sum << "promax in 397" << endl;
    for(int diag = minlag+1;diag <= subcount;diag ++)
    {
        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() + diag - 1,
                                                                 UTS.data(),LTS.data(),subseqlen));

        dot_TS_M = (elementWiseMultiply_p_sum(ts.data()+ diag - 1,
                                                    MASK.data(),subseqlen));

        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                                     MASK.data(),subseqlen));
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
    debug_sum = sum_bool_Matrix(lb_profile_dummy);
    cout << debug_sum << "promax in 450" << endl;
    for(int diag = minlag+1;diag <= subcount;diag ++)
    {
        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data(),
                                                                 UTS.data() + diag - 1,LTS.data()+ diag - 1,subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data(),
                                                    MASK.data() + diag - 1,subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                                     MASK.data() + diag - 1,subseqlen));
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

    debug_sum = sum_bool_Matrix(lb_profile_dummy);
    
    for (int row = 1; row <= subcount; row++) {
        for (int col = row; col <= min(row + minlag, subcount); col++) {
            lb_profile_dummy[row - 1][col - 1] = true;
        }
    }

    lb_profile_dummy = matrixOr(lb_profile_dummy, transposeMatrix_bool(lb_profile_dummy));

    for(int row = 1;row <=subcount;row++)
    {
        vector<DOUBLE> LTS_temp = LTS;
        vector<DOUBLE> UTS_temp = UTS;

        for ( i = row; i <= row + subseqlen - 1; i++) {
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
    debug_sum = sum_bool_Matrix(lb_profile_dummy);
    
    SSS2 = sum_bool_Matrix(lb_profile_dummy);
    aaa = subcount * subcount;
    Pratio = SSS2/aaa;
    auto end_time2 = std::chrono::high_resolution_clock::now();
    auto tmp_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time);

    printf("Factor:%5.2f,ths_diff: %5.2f, Selected V6 plus_pro MPI purns :%5.6f , cost: %5.3fs \n",adjust_factor,ths_diff,Pratio,DOUBLE(tmp_time.count()/1000000));
    fflush(stdout);

    result.Pratio = Pratio;
    result.lb_profile_dummy = lb_profile_dummy;
    result.MPI = MPI;
    result.MP = MP;

}

