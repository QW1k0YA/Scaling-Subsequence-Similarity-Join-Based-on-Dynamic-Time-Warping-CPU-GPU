
#include "iostream"
#include "ctime"
#include "vector"
#include "cmath"
#include "alldef/matrix.h"
#include "alldef/elseoperation.h"
#include "algorithm"
#include "chrono"
#include "allunder/underdtw.h"
using namespace std;

vector<vector<bool>> dist_mpx_v2O1_fast6(const vector<DOUBLE> &ts, int minlag, int subseqlen, int warpmax, DOUBLE bsf){

    auto start_time = std::chrono::high_resolution_clock::now();

    int subcount = ts.size() - subseqlen + 1;
    vector<vector<bool>> lb_profile_dummy = falseMatrix(subcount, subcount);

    vector<vector<DOUBLE>> ts_row = convertTo2DRowOrder(ts);

    vector<vector<DOUBLE>> ts_col = transposeMatrix_double(ts_row);
    
    subcount = ts.size() - subseqlen + 1;

    vector<DOUBLE> mu = movmean(ts, 0, subseqlen - 1, 1);
    vector<DOUBLE> sig = movstd(ts, 0, subseqlen - 1, 1);

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);

    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);

    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> UTS, LTS;
    if (warpmax >= 0) {
        UTS = movmax(ts, warpmax, warpmax);
        LTS = movmin(ts, warpmax, warpmax);
    }

    vector<DOUBLE> dr_bwdU, dr_bwdL, dc_bwd;
    dr_bwdU = addElementToFront(extr_vfromv(UTS, 1, subcount - 1), 0.0);
    dr_bwdL = addElementToFront(extr_vfromv(LTS, 1, subcount - 1), 0.0);
    dc_bwd = addElementToFront(extr_vfromv(ts, 1, subcount - 1), 0.0);

    vector<DOUBLE> dr_fwdU, dr_fwdL, dc_fwd;
    dr_fwdU = extr_vfromv(UTS, subseqlen, UTS.size());
    dr_fwdL = extr_vfromv(LTS, subseqlen, LTS.size());
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

    vector<DOUBLE> del(subcount,0);

    vector<DOUBLE> normLTS(subseqlen, 0.0), normUTS(subseqlen, 0.0);
    DOUBLE DUL2, DUL, del_ths;

    for (int row = 1; row <= subcount; row++) {

        vector<DOUBLE> LTS_temp(LTS);
        vector<DOUBLE> UTS_temp(UTS);

        for (int i = row; i <= row + subseqlen - 1; i++) {
            LTS_temp[i - 1] -= mu[row - 1];
            normLTS[i - row] = LTS_temp[i - 1] / sig[row - 1];
            UTS_temp[i - 1] -= mu[row - 1];
            normUTS[i - row] = UTS_temp[i - 1] / sig[row - 1];
        }

        DUL2 = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL = pow(MAX(DUL2, 0.0), 0.5);
        del_ths = (pow((2 * bsf + DUL), 2) + DUL * DUL) * 0.5;
        del[row - 1] = (del_ths - 2 * subseqlen) * sig[row - 1] - 2 * subseqlen * pow(mu[row - 1], 2) * invsig[row - 1];
    }

    vector<DOUBLE> norm_U_plus_norm_L_trans(subcount, 0.0);

    for (int row = 0; row < subcount; row++) {
        norm_U_plus_norm_L_trans[row] = (sumU2_sumL2[row] - 2 * sumU_sumL[row] * mu[row]) * invsig[row];
    }

    DOUBLE cov_U_plus_cov_L;

    int debug  = 0;

    for (int diag = minlag + 1; diag <= subcount; diag++) {

        vector<bool> lb_vector(subcount,0);

        debug += sum_vector_bool(lb_vector);

    }

    cout << "debug = " << debug << endl;

    size_t temp;
    for (int row = 1; row <= subcount; row++) {
        temp = MIN(row + minlag, subcount);
        for (int col = row; col <= temp; col++) {
            lb_profile_dummy[row - 1][col - 1] = true;
        }
    }

    lb_profile_dummy = matrixOr(lb_profile_dummy, transposeMatrix_bool(lb_profile_dummy));

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    std::cout << "None-selected cost : " << DOUBLE(duration.count())/1000000.0 << " s" << std::endl;

    int SSS2;

    SSS2 = sum_bool_Matrix(lb_profile_dummy);
    DOUBLE aaa = subcount *subcount;
    DOUBLE purned_pieces_dived_the_total = (SSS2+0.0) / aaa;

    printf("None-selected V6 purns :%5.6f, costs  :%5.3f s \n", purned_pieces_dived_the_total, DOUBLE(duration.count())/1000000);

    return lb_profile_dummy;
}