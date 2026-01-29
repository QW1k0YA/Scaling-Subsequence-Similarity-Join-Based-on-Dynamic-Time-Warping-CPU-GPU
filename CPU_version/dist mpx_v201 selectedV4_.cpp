
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
#include "alldef/elseoperation.h"
using namespace std;

void dist_mpx_v2O1_selectedV4_(const vector<DOUBLE>& ts, int minlag, int subseqlen, int warpmax, DOUBLE bsf,
                               double adjust_factor, RETURN_V4& result)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    auto subcount = ts.size() - subseqlen + 1;
    vector<vector<bool>> lb_profile_dummy;
    lb_profile_dummy = falseMatrix(subcount,subcount);
    subcount = ts.size() - subseqlen + 1;

    vector<DOUBLE> mu,sig;
    mu = movmean(ts,0,subseqlen-1,1);
    sig = movstd(ts,0,subseqlen-1,1);
    
    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);

    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);

    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> MP(subcount,INFINITY);
    vector<DOUBLE> MPI(subcount,INFINITY);
    vector<DOUBLE> row_min_DEL(subcount,INFINITY);

    vector<DOUBLE> UTS, LTS;
    if (warpmax >= 0) {
        UTS = movmax(ts, warpmax, warpmax);
        LTS = movmin(ts, warpmax, warpmax);
    }
    vector<DOUBLE> DIFF = substractvector(UTS,LTS);
    DOUBLE avg_diff = sum_vector(DIFF)/DIFF.size();
    DOUBLE std_diff = stddev(DIFF);
    DOUBLE ths_diff = avg_diff+adjust_factor*std_diff;

    vector<DOUBLE> MASK;
    for(DOUBLE value:DIFF)
    {
        if(value <= ths_diff)
        {
            MASK.push_back(1);
        }
        else
        {
            MASK.push_back(0);
        }
    }

    UTS = elementWiseMultiply(UTS,MASK);
    LTS = elementWiseMultiply(LTS,MASK);

    vector<DOUBLE> dr_bwdU, dr_bwdL, dc_bwd;
    dr_bwdU = addElementToFront(extr_vfromv(UTS, 1, subcount - 1), 0.0);
    dr_bwdL = addElementToFront(extr_vfromv(LTS, 1, subcount - 1), 0.0);
    dc_bwd = addElementToFront(extr_vfromv(ts, 1, subcount - 1), 0.0);

    vector<DOUBLE> dr_fwdU, dr_fwdL, dc_fwd;
    dr_fwdU = extr_vfromv(UTS, subseqlen, UTS.size());
    dr_fwdL = extr_vfromv(LTS, subseqlen, UTS.size());
    dc_fwd = extr_vfromv(ts, subseqlen, ts.size());

    vector<DOUBLE> dr_bwdU_plus_dr_bwdL, dr_fwdU_plus_dr_fwdL;
    dr_bwdU_plus_dr_bwdL = plusvector(dr_bwdU, dr_bwdL);
    dr_fwdU_plus_dr_fwdL = plusvector(dr_fwdU, dr_fwdL);

    vector<DOUBLE> sumU, sumL, sumU2, sumL2, sumU2_sumL2, sumU_sumL;

    sumU = movsum(UTS, subseqlen - 1);
    sumL = movsum(LTS, subseqlen - 1);
    sumU2 = movsum(elementWiseMultiply(UTS, UTS), subseqlen - 1);
    sumL2 = movsum(elementWiseMultiply(LTS, LTS), subseqlen - 1);
    sumU2_sumL2 = plusvector(sumU2, sumL2);
    sumU_sumL = plusvector(sumU, sumL);

    vector<DOUBLE> TS2 = elementWiseMultiply(ts,ts);
    vector<DOUBLE> sumMASK = movsum(MASK, subseqlen - 1);
    vector<DOUBLE> norm_U_plus_norm_L(subcount,0.0);

    for(int row = 1; row <=subcount; row++)
    {
        norm_U_plus_norm_L[row-1] = (sumU2_sumL2[row - 1] - 2*sumU_sumL[row-1]*mu[row-1]+2*sumMASK[row-1]*mu[row-1]*mu[row-1])*invsig[row-1]*invsig[row-1];
    }

    vector<DOUBLE> del(subcount,0);
    vector<DOUBLE> normLTS(subseqlen, 0.0), normUTS(subseqlen, 0.0);
    DOUBLE DUL2, DUL, del_ths;
    for (int row = 1; row <= subcount; row++) {

        vector<DOUBLE> LTS_temp(LTS);
        vector<DOUBLE> UTS_temp(UTS);

        unsigned long long temp1 = row + subseqlen - 1;
        for (int i = row; i <= temp1; i++) {
            LTS_temp[i - 1] -= mu[row - 1];
            normLTS[i - row] = LTS_temp[i - 1] / sig[row - 1];
            UTS_temp[i - 1] -= mu[row - 1];
            normUTS[i - row] = UTS_temp[i - 1] / sig[row - 1];
        }

        DUL2 = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * elementWiseMultiply_sum(normLTS,normUTS);
        DUL = sqrt(MAX(DUL2, 0.0));
        del_ths = (pow((2 * bsf + DUL), 2) + DUL * DUL) * 0.5;
        del[row - 1] = del_ths - norm_U_plus_norm_L[row-1];
    }

    vector<DOUBLE> dr_bwdMASK = addElementToFront(extr_vfromv(MASK, 1, subcount - 1), 0);
    vector<DOUBLE> dr_fwdMASK = extr_vfromv(MASK, subseqlen, MASK.size());
    vector<DOUBLE> dc_bwdTS2 = addElementToFront(extr_vfromv(TS2, 1, subcount - 1), 0);
    vector<DOUBLE> dc_fwdTS2 = extr_vfromv(TS2, subseqlen, TS2.size());

    DOUBLE cov_U_plus_cov_L,dot_TS_M,dot_TS2_M;
    for(int diag = minlag+1;diag<= subcount;diag++)
    {

        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() + diag - 1,
                                                                 UTS.data(),LTS.data(),subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data()+ diag - 1,
                                                    MASK.data(),subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data() + diag - 1,
                                                     MASK.data(),subseqlen));
        int row = 1;
        int col = diag;
        DOUBLE M_NORM = dot_TS2_M - 2*mu[col-1]*dot_TS_M + sumMASK[row-1]* mu[col-1]*mu[col-1];
        M_NORM = M_NORM*invsig[col-1]*invsig[col-1];

        DOUBLE t = (cov_U_plus_cov_L - mu[col-1]*sumU_sumL[row-1] + 2*mu[row-1]*(sumMASK[row-1]*mu[col-1] - dot_TS_M))*(invsig[col-1]*invsig[row-1]);
        DOUBLE dist = 2 *M_NORM - 2*t;

        if(dist>del[row-1])
        {
            lb_profile_dummy[row-1][col-1] = true;
        }
        else if(dist < row_min_DEL[row-1])
        {
            row_min_DEL[row-1] = dist;
            MPI[row-1] = col;
        }

        auto temp_1 = subcount - diag + 1;
        for(int row = 2;row <= temp_1;row++)
        {
            int col = diag + row -1;
            cov_U_plus_cov_L = cov_U_plus_cov_L - dr_bwdU_plus_dr_bwdL[row-1]*dc_bwd[col-1] + (dr_fwdU_plus_dr_fwdL[row-1])*dc_fwd[col-1];
            dot_TS_M = dot_TS_M - dr_bwdMASK[row-1]*dc_bwd[col-1] + dr_fwdMASK[row-1]*dc_fwd[col-1];
            dot_TS2_M = dot_TS2_M - dr_bwdMASK[row-1]*dc_bwdTS2[col-1]+dr_fwdMASK[row-1]*dc_fwdTS2[col-1];
            M_NORM = dot_TS2_M - 2*mu[col-1]*dot_TS_M + sumMASK[row-1]*mu[col-1]*mu[col-1];
            M_NORM = M_NORM*(invsig[col-1]*invsig[col-1]);
            t = (cov_U_plus_cov_L - mu[col-1]*sumU_sumL[row-1]+2*mu[row-1]*(sumMASK[row-1]*mu[col-1] - dot_TS_M))*(invsig[col-1]*invsig[row-1]);
            dist = 2*M_NORM - 2*t;
            if(dist>del[row-1])
            {
                lb_profile_dummy[row-1][col-1] = true;

            }
            else if(dist < row_min_DEL[row-1])
            {
                row_min_DEL[row-1] = dist;
                MPI[row-1] = col;
            }
        }

    }

    for(int diag = minlag + 1;diag <= subcount;diag ++)
    {

        cov_U_plus_cov_L = (elementWiseMultiply_p_plus_sum(ts.data() ,UTS.data()+ diag - 1,
                                                                 LTS.data()+ diag - 1,subseqlen));
        dot_TS_M = (elementWiseMultiply_p_sum(ts.data(),
                                                    MASK.data() + diag - 1,subseqlen));
        dot_TS2_M = (elementWiseMultiply_p_sum(TS2.data(),
                                                     MASK.data() + diag - 1,subseqlen));
        int col =1;
        int row = diag;
        DOUBLE  M_NORM=dot_TS2_M - 2 * mu[col-1] * dot_TS_M + sumMASK[row-1] * mu[col-1]*mu[col-1];
        M_NORM = M_NORM*(invsig[col-1]*invsig[col-1]);
        DOUBLE t = (cov_U_plus_cov_L - mu[col-1]*sumU_sumL[row-1] + 2*mu[row-1]*(sumMASK[row-1]*mu[col-1] - dot_TS_M))*invsig[col-1]*invsig[row-1];
        DOUBLE dist = 2*M_NORM - 2*t;

        if(dist >del[row-1])
        {
            lb_profile_dummy[col-1][row-1] = true;
        }
        else if (dist <row_min_DEL[row-1])
        {

            row_min_DEL[row-1] = dist;
            MPI[row-1] = col;
        }

        auto temp_2 =  subcount - diag + 1;
        for(col = 2;col <= temp_2;col++)
        {
            row = diag +col -1;
            cov_U_plus_cov_L = cov_U_plus_cov_L - (dr_bwdU_plus_dr_bwdL[row-1])*dc_bwd[col-1] + (dr_fwdU_plus_dr_fwdL[row-1])*dc_fwd[col-1];
            dot_TS_M = dot_TS_M - dr_bwdMASK[row-1]*dc_bwd[col-1] + dr_fwdMASK[row-1]*dc_fwd[col-1];
            dot_TS2_M = dot_TS2_M - dr_bwdMASK[row-1]*dc_bwdTS2[col-1] + dr_fwdMASK[row-1]*dc_fwdTS2[col-1];
            M_NORM = dot_TS2_M - 2*mu[col-1]*dot_TS_M +sumMASK[row-1]*mu[col-1]*mu[col-1];
            M_NORM = M_NORM*invsig[col-1]*invsig[col-1];

            t = (cov_U_plus_cov_L - mu[col-1]*sumU_sumL[row-1] + 2*mu[row-1]*(sumMASK[row-1]*mu[col-1] - dot_TS_M))*(invsig[col-1]*invsig[row-1]);
            dist = 2*M_NORM - 2*t;

            if(dist >del[row-1])
            {
                lb_profile_dummy[col-1][row-1] = true;
            }
            else if(dist< row_min_DEL[row-1])
            {
                row_min_DEL[row-1] = dist;
                MPI[row-1] = col;
            }

        }
    }

    for(size_t row = 1;row <= subcount;row++)
    {
        auto temp_3 = MIN(row + minlag, subcount);
        for(size_t col = row;col <= temp_3;col++)
        {
            lb_profile_dummy[row-1][col-1] = true;
        }
    }

    lb_profile_dummy = matrixOr(lb_profile_dummy, transposeMatrix_bool(lb_profile_dummy));

    for(size_t row = 1;row <=subcount;row++)
    {
        vector<DOUBLE> LTS_temp(LTS);
        vector<DOUBLE> UTS_temp(UTS);

        auto temp_4 = row + subseqlen - 1;
        for (size_t i = row; i <= temp_4; i++) {
            LTS_temp[i - 1] -= mu[row - 1];
            normLTS[i - row] = LTS_temp[i - 1] / sig[row - 1];
            UTS_temp[i - 1] -= mu[row - 1];
            normUTS[i - row] = UTS_temp[i - 1] / sig[row - 1];
        }

        DUL2 = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL = pow(MAX(DUL2, 0.0), 0.5);

        row_min_DEL[row-1] = row_min_DEL[row-1] + norm_U_plus_norm_L[row-1];
        row_min_DEL[row-1] = row_min_DEL[row-1]*2;

        MP[row-1] = (sqrt(row_min_DEL[row-1]-DUL2)-DUL)/2;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    std::cout << "selected cost : " << duration.count()/1000000 << " s" << std::endl;
    DOUBLE SSS2 = sum_bool_Matrix(lb_profile_dummy);
    DOUBLE aaa = subcount*subcount;
    DOUBLE Pratio = SSS2/aaa;

    printf("Factor:%5.2f Selected V4 MPI purns :%5.6f cost: %5.3fs \n", adjust_factor, Pratio,DOUBLE(duration.count())/1000000);

    result.lb_profile_dummy = lb_profile_dummy;
    result.Pratio = Pratio;
    result.MP = MP;
    result.MPI = MPI;

}

