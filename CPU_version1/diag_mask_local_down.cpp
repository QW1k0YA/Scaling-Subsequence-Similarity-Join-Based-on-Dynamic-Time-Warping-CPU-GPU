
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define DOUBLE_BIAS 6
#define BIAS 3

void
diag_mask_local_down(const vector<DOUBLE> &ts, int minlag, int subseqlen, DOUBLE bsf, int diagID,
                     vector<bool> &lb_vector, const vector<vector<SHORT>> &count_table,
                     const vector<DOUBLE> &pos_UU, const vector<DOUBLE> &pos_LL, const vector<DOUBLE> &mu,
                     const vector<DOUBLE> &UTS_, const vector<DOUBLE> &LTS_, double **special_shared_vector,
                     double &cnt, vector<DOUBLE> &sumU2_sumL2, vector<DOUBLE> &sumU_sumL, vector<DOUBLE> &sumMASK,
                     vector<DOUBLE> &norm_U_plus_norm_L, vector<DOUBLE> &raw_DIFF_UL, vector<DOUBLE> &raw_DIFF_UL2,
                     vector<DOUBLE> &raw_DIFF_UL2_temp, vector<DOUBLE> &DUL2_raw, vector<DOUBLE> &DUL2,
                     const vector<DOUBLE> &TS2, double alpha, const vector<double> &invsig,
                     const vector<double> &invsig_2, int diag_bias)
{
    
    size_t len = ts.size();
    size_t subcount = len - subseqlen + 1;
    int diag = diagID;

    if(diagID < 0){
        cerr << "diag_ID must >= 0 in distLB";
        exit(0);
    } else if (diagID <= minlag){
        cerr << "ONLY IF diag_ID  > minlag , IT HAS MEANINGS  in distLB";
        exit(0);
    } else if(diagID > subcount){
        cerr << "diag_ID CANNOT > subcount  in distLB" ;
        exit(0);
    }

    vector<DOUBLE> normal_DIFF(len,0.0);
    size_t posLL,posUU;
    int pos;
    DOUBLE diff2,diff_ratio;
    
    for(size_t i = diag_bias;i < len;i++)
    {
        pos = i - diagID +1;
        if(pos < 0 )
        {
            continue;
        }

        posLL = pos_LL[i]; 
        posUU = pos_UU[i]; 
        diff2 = count_table[pos][posUU-1] - count_table[pos][posLL-1]; 
       
        diff_ratio = diff2/subseqlen;  
        
        normal_DIFF[i] = diff_ratio;
    }

    DOUBLE avg_diff = meanWithoutNaN_from_pos(normal_DIFF,diag_bias);
    DOUBLE ths_diff = (avg_diff + 0.5)/2 + alpha;

    substractvector_from_pos(UTS_,LTS_,raw_DIFF_UL,diag_bias);
    elementWiseMultiply_from_pos(raw_DIFF_UL,raw_DIFF_UL,raw_DIFF_UL2_temp,diag_bias);

    vector<DOUBLE> MASK(normal_DIFF.size(),0);
    int i = 0;

    for(int mpos = diag_bias;mpos < len;mpos++)
    {
        if(normal_DIFF[mpos] <= ths_diff)
        {
            MASK[mpos] = 1;
        }
        else
        {
            MASK[mpos] = 0;
        }
    }

    vector<DOUBLE> MP(subcount,INFINITY);
    vector<DOUBLE> MPI(subcount,INFINITY);
    DOUBLE *dr_bwdU_plus_dr_bwdL_less, *dr_fwdU_plus_dr_fwdL;
    const double* dc_bwd_less = ts.data() + BIAS;
    
    const double *dc_fwd = ts.data() + subseqlen - 1 - BIAS;

    vector<DOUBLE> UTS = elementWiseMultiply_from_pos(UTS_,MASK,diag_bias);
    vector<DOUBLE> LTS = elementWiseMultiply_from_pos(LTS_,MASK,diag_bias);

    dr_bwdU_plus_dr_bwdL_less = plusvector_p_from_pos(UTS.data() + BIAS, LTS.data() + BIAS,UTS.size()-1,diag_bias);
    
    dr_fwdU_plus_dr_fwdL = plusvector_p_from_pos(UTS.data() + subseqlen - BIAS - 1, LTS.data() + subseqlen - BIAS - 1,
                                                 UTS.size() - subseqlen + 1,diag_bias);

    DOUBLE sum_sumU = 0,sum_sumL = 0,sum_sumU2 = 0,sum_sumL2 = 0,sum_sumMASK = 0;
    double *U_pos = UTS.data() +BIAS;
    double *L_pos = LTS.data() +BIAS;
    double *MASK_pos = MASK.data() + BIAS;
    DOUBLE sumU_temp_fir,sumL_temp_fir,sumU_temp_sec,sumL_temp_sec;
    for(i = diag_bias;i < diag_bias + subseqlen - DOUBLE_BIAS;i++)
    { 

        sumU_temp_fir = *(U_pos + i);
        sumL_temp_fir = *(L_pos + i);

        sum_sumU += sumU_temp_fir;
        sum_sumL += sumL_temp_fir;
        sum_sumU2 += sumU_temp_fir * sumU_temp_fir;
        sum_sumL2 += sumL_temp_fir * sumL_temp_fir;
        sum_sumMASK += *(MASK_pos + i);
    }

    for(i = diag_bias;i < subcount - DOUBLE_BIAS;i++)
    {
        sumU_sumL[i] = sum_sumU + sum_sumL;
        sumU2_sumL2[i] = sum_sumU2 + sum_sumL2;
        sumMASK[i] = sum_sumMASK;
        sumU_temp_fir = *(U_pos + i);
        sumL_temp_fir = *(L_pos + i);
        sumU_temp_sec = *(U_pos + i + subseqlen -DOUBLE_BIAS);
        sumL_temp_sec = *(L_pos + i + subseqlen -DOUBLE_BIAS);

        sum_sumU = sum_sumU - sumU_temp_fir + sumU_temp_sec;
        sum_sumL = sum_sumL - sumL_temp_fir + sumL_temp_sec;
        sum_sumU2 = sum_sumU2 - sumU_temp_fir*sumU_temp_fir + sumU_temp_sec*sumU_temp_sec;
        sum_sumL2 = sum_sumL2 - sumL_temp_fir*sumL_temp_fir + sumL_temp_sec*sumL_temp_sec;
        sum_sumMASK = sum_sumMASK - *(MASK_pos + i) + *(MASK_pos + i+subseqlen - DOUBLE_BIAS);
    }

    for(size_t row = diag_bias;row < subcount;row++)
    {
        norm_U_plus_norm_L[row] = (sumU2_sumL2[row] - 2*sumU_sumL[row]*mu[row]
                                     + 2*sumMASK[row]*mu[row]*mu[row])*invsig[row]*invsig[row];
    }

    elementWiseMultiply_from_pos(raw_DIFF_UL2_temp,MASK,raw_DIFF_UL2,diag_bias);
    movsum_p_from_pos(raw_DIFF_UL2.data() + BIAS,raw_DIFF_UL2.size() - DOUBLE_BIAS, subseqlen - 1 - DOUBLE_BIAS,DUL2_raw,diag_bias);
    elementWiseMultiply_from_pos(DUL2_raw,invsig_2,DUL2,diag_bias);

    DOUBLE del_ths;
    DOUBLE *dr_bwdMASK_less = MASK.data() + BIAS;
    DOUBLE *dr_fwdMASK = MASK.data() + subseqlen - BIAS - 1;
    const double *dc_bwdTS2_less = TS2.data() + BIAS;
    const double *dc_fwdTS2 = TS2.data() + subseqlen - BIAS - 1;

    DOUBLE cov_U_plus_cov_L_sec = (elementWiseMultiply_p_plus_sum(ts.data() + 3 , UTS.data() + diag - 1 + 3,
                                                                  LTS.data()+ diag - 1 + 3,subseqlen - 6));
    DOUBLE dot_TS_M_sec = (elementWiseMultiply_p_sum(ts.data() + 3,
                                                     MASK.data() + diag - 1 + 3,subseqlen - 6));
    DOUBLE dot_TS2_M_sec = (elementWiseMultiply_p_sum(TS2.data() + 3,
                                                      MASK.data() + diag - 1 + 3,subseqlen - 6));

    DOUBLE M_NORM,t,dist,DUL_value;
    DOUBLE lb;
    if(!lb_vector[0])
    {
        int col = 0;
        int row = diag - 1;
        M_NORM = dot_TS2_M_sec - 2 * mu[col] * dot_TS_M_sec + sumMASK[row] * mu[col] * mu[col];
        M_NORM = M_NORM * invsig[col]*invsig[col];
        t = (cov_U_plus_cov_L_sec - mu[col] * sumU_sumL[row] +
             2 * mu[row] * (sumMASK[row] * mu[col] - dot_TS_M_sec)) * (invsig[col] * invsig[row]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[row];
        if(dist < DUL2[row]){
            lb = 0;
        }
        else{
            
            DUL_value = sqrt(MAX(DUL2[row],0));
            lb = 0.5*(sqrt(2*dist - DUL2[row]) - DUL_value);
        }
        if(lb > bsf || (lb*lb + special_shared_vector[col][0]+special_shared_vector[col][1] > bsf*bsf)){
            lb_vector[col] = true;
            cnt++;
        }
    }

    DOUBLE bsf_2 = bsf*bsf;
    auto temp_1 = subcount - diag + 1;
    int low_index;
    int high_index;
    int low_index_less;
    int high_index_less;

    for (low_index = 1; low_index < temp_1; low_index++) {
        high_index = diag + low_index - 1;
        low_index_less = low_index - 1;
        high_index_less = high_index - 1;

        cov_U_plus_cov_L_sec = cov_U_plus_cov_L_sec - dr_bwdU_plus_dr_bwdL_less[high_index_less] * dc_bwd_less[low_index_less] +
                               dr_fwdU_plus_dr_fwdL[high_index] * dc_fwd[low_index];
        dot_TS_M_sec = dot_TS_M_sec - dr_bwdMASK_less[high_index_less] * dc_bwd_less[low_index_less] + dr_fwdMASK[high_index] * dc_fwd[low_index];
        dot_TS2_M_sec = dot_TS2_M_sec - dr_bwdMASK_less[high_index_less] * dc_bwdTS2_less[low_index_less] + dr_fwdMASK[high_index] * dc_fwdTS2[low_index];

        if(lb_vector[low_index])
        {
            continue;
        }

        M_NORM = dot_TS2_M_sec - 2 * mu[low_index] * dot_TS_M_sec + sumMASK[high_index] * mu[low_index] * mu[low_index];
        M_NORM = M_NORM * invsig[low_index] * invsig[low_index];
        t = (cov_U_plus_cov_L_sec - mu[low_index] * sumU_sumL[high_index] +
             2 * mu[high_index] * (sumMASK[high_index] * mu[low_index] - dot_TS_M_sec)) * (invsig[low_index] * invsig[high_index]);
        dist = 2 * M_NORM - 2 * t + norm_U_plus_norm_L[high_index];

        if(dist < DUL2[high_index])
        {
            lb = 0;
        }
        else
        {
            DUL_value = sqrt(MAX(DUL2[high_index],0));
            lb = 0.5*(sqrt(2*dist - DUL2[high_index]) - DUL_value);
            if(lb > bsf || (lb*lb + special_shared_vector[low_index][0] + special_shared_vector[low_index][1] > bsf_2))
            {
                cnt++;
                lb_vector[low_index] = true;
                cnt_down++;
            }
        }
    }

    delete[] dr_bwdU_plus_dr_bwdL_less;
    delete[] dr_fwdU_plus_dr_fwdL;
}

