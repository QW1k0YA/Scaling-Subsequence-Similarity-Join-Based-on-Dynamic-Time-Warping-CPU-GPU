
#include <cstring>
#include <iomanip>
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
#include "alldef/typedefdouble.h"
#include "allunder/underdtw.h"
#include "allunder/underDIAGV3.h"
#include "allunder/underpromax.h"
#include "alldef/allclass.h"
#define DOUBLE_BIAS 6
#define BIAS 3
using namespace  std;
double cnt_up = 0;
double cnt_down = 0;
void new_dtw_motifGUI(const vector<DOUBLE> &a, int subseqLen, int maxwarp, const vector<DOUBLE> &mp_ed,
                      const vector<int> &mp_ed_index, const char *output_time_file_path) {
    int subcount = a.size() - subseqLen + 1;
    int len = a.size();
    cnt_down = 0;
    cnt_up = 0;

    vector<DOUBLE> mu(subcount), sig(subcount);
    mvmean(a, subseqLen, mu, sig);
    vector<DOUBLE> invsig(subcount);
    vector<DOUBLE> invsig_2(subcount);
    for (int i = 0; i < subcount; i++) {
        invsig[i] = 1 / sig[i];
        invsig_2[i] = invsig[i] * invsig[i];
    }

    vector<vector<DOUBLE>> my_subs(subcount, vector<DOUBLE>(subseqLen, 0.0)); 
    vector<vector<DOUBLE>> subs_U(subcount, vector<DOUBLE>(subseqLen, 0.0));
    vector<vector<DOUBLE>> subs_L(subcount, vector<DOUBLE>(subseqLen, 0.0));
    size_t temp_1;
    for (int i = 0; i < subcount; i++) {
        for (int j = 0; j < subseqLen; j++) {
            my_subs[i][j] = (a[i + j] - mu[i]) / sig[i];
        }
        lower_upper_lemire(my_subs[i], subseqLen, maxwarp, subs_L[i], subs_U[i]);
    }

    auto my_total_time = std::chrono::high_resolution_clock::now();
    printf("buffer_length: %d, maxwarp: %d\n", subseqLen, maxwarp);

    int minspacing = subseqLen;

    int k = 100;
    vector<int> min_index = min_v_k(mp_ed, k);
    
    double best_so_far = mp_ed[min_index[0]];
    std::vector<double> cb_temp(subseqLen + maxwarp + 1, 0);

    int low_index;
    int high_index;
    for (int i = 0; i < k; i++) {
        low_index = min_index[i];
        high_index = mp_ed_index[min_index[i]] - 1;

        if (LB_KK_FIRST(my_subs[low_index], my_subs[high_index], subs_U[high_index], subs_L[low_index], subseqLen,
                        best_so_far * best_so_far, cb_temp)) {
            double dtw_value = sqrt(MON_dtw(my_subs[low_index], my_subs[high_index], cb_temp,
                                            subseqLen, maxwarp, best_so_far * best_so_far));
            if (dtw_value < best_so_far) {
                best_so_far = dtw_value;
            }
        }

    }

    DOUBLE first_min = MIN(low_index, high_index);
    DOUBLE sec_min = MAX(low_index, high_index);

    cout << "the buffer_length of TS is " << a.size() << endl;
    printf("################Start calculating lower bound#################  \n");

    int debug_sum = 0;
    DOUBLE bsf = best_so_far;

    vector<bool> lb_vector(subcount, false);
    
    vector<vector<DOUBLE>> ts_row = convertTo2DRowOrder(a);
    vector<vector<DOUBLE>> ts_col = transposeMatrix_double(ts_row);

    subcount = a.size() - subseqLen + 1;

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);

    int warpmax = maxwarp;

    vector<DOUBLE> UTS(len);
    vector<DOUBLE> LTS(len);
    lower_upper_lemire(a, len, warpmax, LTS, UTS);

    vector<DOUBLE> UTS_p(subseqLen);
    vector<DOUBLE> LTS_p(subseqLen);

    vector<DOUBLE> dr_bwdU, dr_bwdL, dc_bwd;
    dr_bwdU = addElementToFront(extr_vfromv(UTS, 1, subcount - 1), 0.0);
    dr_bwdL = addElementToFront(extr_vfromv(LTS, 1, subcount - 1), 0.0);
    dc_bwd = addElementToFront(extr_vfromv(a, 1, subcount - 1), 0.0);

    vector<DOUBLE> dr_fwdU, dr_fwdL, dc_fwd;
    dr_fwdU = extr_vfromv(UTS, subseqLen, len);
    dr_fwdL = extr_vfromv(LTS, subseqLen, len);
    dc_fwd = extr_vfromv(a, subseqLen, len);

    vector<DOUBLE> dr_bwdU_plus_dr_bwdL, dr_fwdU_plus_dr_fwdL;
    dr_bwdU_plus_dr_bwdL = plusvector(dr_bwdU, dr_bwdL);
    dr_fwdU_plus_dr_fwdL = plusvector(dr_fwdU, dr_fwdL);

    const vector<DOUBLE> sumU = movsum(UTS, subseqLen - 1);
    const vector<DOUBLE> sumL = movsum(LTS, subseqLen - 1);
    const vector<DOUBLE> sumU2 = movsum(elementWiseMultiply(UTS, UTS), subseqLen - 1);
    const vector<DOUBLE> sumL2 = movsum(elementWiseMultiply(LTS, LTS), subseqLen - 1);
    const vector<DOUBLE> sumU2_sumL2 = plusvector(sumU2, sumL2);
    const vector<double> sumU_sumL = plusvector(sumU, sumL);

    vector<DOUBLE> del(subcount, 0);

    vector<DOUBLE> normLTS(subseqLen, 0.0), normUTS(subseqLen, 0.0);
    vector<DOUBLE> DUL2_fast(subcount);
    vector<DOUBLE> DUL_fast(subcount);

    DOUBLE del_ths;
    for (int row = 0; row < subcount; row++) {
        for (int i = row; i <= row + subseqLen - 1; i++) {
            normLTS[i - row] = (LTS[i] - mu[row]) * invsig[row];
            normUTS[i - row] = (UTS[i] - mu[row]) * invsig[row];
        }

        DUL2_fast[row] = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
               2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL_fast[row] = pow(MAX(DUL2_fast[row], 0.0), 0.5);

    }

        vector<DOUBLE> norm_U_plus_norm_L_trans(subcount, 0.0);

        for (int row = 0; row < subcount; row++) {
            norm_U_plus_norm_L_trans[row] = (sumU2_sumL2[row] - 2 * sumU_sumL[row] * mu[row]) * invsig[row];
        }

        int len_of_table_local = INITIAL_LEN_OF_TABLE;

        DOUBLE real_max, real_min;
        vector<DOUBLE> pos_UU(len, 1.0);
        vector<DOUBLE> pos_LL(len, 1.0);
        int len_of_cdf;
        real_max = my_subs[0][0];
        real_min = my_subs[0][0];
        for (const auto &value1: my_subs) {
            for (auto value2: value1) {
                if (value2 > real_max) real_max = value2;
                if (value2 < real_min) real_min = value2;
            }
        }
        
        double MAX_REAL_VALUE = MAX(abs(real_max), abs(real_min));
        if (MAX_REAL_VALUE > 3) {
            len_of_cdf = (MAX_REAL_VALUE - 3) * 20 + 600;
        } else {
            len_of_cdf = 600;
        }

        vector<vector<SHORT>> count_table_local(len, vector<SHORT>(static_cast<int>(len_of_cdf), 0));
        compute_shared_data_local(a, subseqLen, my_subs, subs_U, subs_L, real_max,
                                  real_min, pos_UU,
                                  pos_LL, len_of_cdf,
                                  count_table_local, MAX_REAL_VALUE); 

        DOUBLE cov_U_plus_cov_L;

        DOUBLE pruning_rate = 0;
        DOUBLE thr = 0.9;

        int minlag = subseqLen;

        vector<DOUBLE> proj(subseqLen);
        double lbk1, lbk2;

        vector<double> cb(subseqLen + warpmax + 1, 0);
        vector<double> cb1(subseqLen + warpmax);
        vector<double> cb2(subseqLen + warpmax);

        vector<DOUBLE> MP(subcount, INFINITY);
        vector<DOUBLE> MPI(subcount, INFINITY);
        vector<DOUBLE> row_min_DEL(subcount, INFINITY);

        DOUBLE max_real_value = floor(MAX(abs(real_max), abs(real_min))) * 2 * 100;
        int len_of_table_global = static_cast<int>(MAX(INITIAL_LEN_OF_TABLE,
                                                       max_real_value)); 
        vector<vector<DOUBLE>> count_table_global(subseqLen, vector<DOUBLE>(len_of_table_global, 0.0));
        
        printf("LEN_OF_TABLE_GLOBAL is %d \n", len_of_table_global);
        compute_shared_data_global(a, subseqLen, my_subs, len_of_table_global, count_table_global);

        vector<DOUBLE> MASK_global(len, 0.0);
        DOUBLE qqq, www;
        vector<DOUBLE> TS2 = elementWiseMultiply(a, a);
        compute_MASK_global(a, subseqLen, my_subs, len_of_table_global, count_table_global, subs_U, subs_L,
                            MASK_global);

        int i;
        vector<DOUBLE> UTS_global = elementWiseMultiply(UTS, MASK_global);
        vector<DOUBLE> LTS_global = elementWiseMultiply(LTS, MASK_global);

        vector<DOUBLE> dr_bwdU_global = addElementToFront(extr_vfromv(UTS_global, 1 + BIAS, UTS.size() - 1), 0.0);
        vector<DOUBLE> dr_bwdL_global = addElementToFront(extr_vfromv(LTS_global, 1 + BIAS, LTS.size() - 1), 0.0);

        vector<DOUBLE> dr_fwdU_global = extr_vfromv(UTS_global, subseqLen - BIAS, UTS_global.size());
        vector<DOUBLE> dr_fwdL_global = extr_vfromv(LTS_global, subseqLen - BIAS, LTS_global.size());

        vector<DOUBLE> dc_bwd_global = addElementToFront(extr_vfromv(a, 1 + BIAS, len), 0.0);
        vector<DOUBLE> dc_fwd_global = extr_vfromv(a, subseqLen - BIAS, a.size());

        vector<DOUBLE> dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global;
        dr_bwdU_plus_dr_bwdL_global = plusvector(dr_bwdU_global, dr_bwdL_global);
        dr_fwdU_plus_dr_fwdL_global = plusvector(dr_fwdU_global, dr_fwdL_global);

        const vector<DOUBLE> sumU_global = movsum_p(UTS_global.data() + BIAS, len - DOUBLE_BIAS,
                                                    subseqLen - 1 - DOUBLE_BIAS);
        const vector<DOUBLE> sumL_global = movsum_p(LTS_global.data() + BIAS, len - DOUBLE_BIAS,
                                                    subseqLen - 1 - DOUBLE_BIAS);
        const vector<DOUBLE> sumU2_global = movsum(
                elementWiseMultiply_p(UTS_global.data() + BIAS, UTS_global.data() + BIAS, len - DOUBLE_BIAS),
                subseqLen - 1 - DOUBLE_BIAS);
        const vector<DOUBLE> sumL2_global = movsum(
                elementWiseMultiply_p(LTS_global.data() + BIAS, LTS_global.data() + BIAS, len - DOUBLE_BIAS),
                subseqLen - 1 - DOUBLE_BIAS);
        const vector<DOUBLE> sumU2_sumL2_global = plusvector(sumU2_global, sumL2_global);
        const vector<DOUBLE> sumU_sumL_global = plusvector(sumU_global, sumL_global);

        const vector<DOUBLE> sumMASK_global = movsum_p(MASK_global.data() + BIAS, MASK_global.size() - DOUBLE_BIAS,
                                                       subseqLen - 1 - DOUBLE_BIAS);
        vector<DOUBLE> norm_U_plus_norm_L_global(subcount, 0.0);

        for (int row = 1; row <= subcount; row++) {
            norm_U_plus_norm_L_global[row - 1] =
                    (sumU2_sumL2_global[row - 1] - 2 * sumU_sumL_global[row - 1] * mu[row - 1]
                     + 2 * sumMASK_global[row - 1] * mu[row - 1] * mu[row - 1]) * invsig[row - 1] * invsig[row - 1];
        }

        vector<DOUBLE> del_global(subcount, 0);
        vector<DOUBLE> normLTS_global(subseqLen, 0.0), normUTS_global(subseqLen, 0.0);
        DOUBLE del_ths_global;
        vector<DOUBLE> DUL2_global(subcount);
        vector<DOUBLE> DUL_global(subcount);

        for (int row = 0; row < subcount; row++) {
            for (i = row; i <= row + subseqLen - 1; i++) {

                normLTS_global[i - row] = (LTS_global[i] - mu[row]) * invsig[row];
                normUTS_global[i - row] = (UTS_global[i] - mu[row]) * invsig[row];
            }

            DUL2_global[row] = pow(norm_vector(substractvector(normLTS_global, normUTS_global), 2), 2);
            
            DUL_global[row] = sqrt(MAX(DUL2_global[row], 0));

        }

        const vector<DOUBLE> dr_bwdMASK_global = addElementToFront(
                extr_vfromv(MASK_global, 1 + BIAS, MASK_global.size()), 0.0);
        const vector<DOUBLE> dr_fwdMASK_global = extr_vfromv(MASK_global, subseqLen - BIAS,
                                                             MASK_global.size());

        const vector<DOUBLE> dc_bwdTS2_global = addElementToFront(extr_vfromv(TS2, 1 + BIAS, TS2.size() - 1), 0.0);
        const vector<DOUBLE> dc_fwdTS2_global = extr_vfromv(TS2, subseqLen - BIAS, TS2.size());

        vector<DOUBLE> lb_vector_new(subcount);
        double **special_shared_vector = (double **) malloc(subcount * sizeof(double *));
        if (special_shared_vector == nullptr) {
            fprintf(stderr, "error\n");
        }

        for (i = 0; i < subcount; i++) {
            special_shared_vector[i] = (double *) malloc(2 * sizeof(double));
            if (special_shared_vector[i] == nullptr) {
                fprintf(stderr, "error\n");
                
                for (int j = 0; j < i; j++) {
                    free(special_shared_vector[j]);
                }
                free(special_shared_vector);
            }
        }
        for (i = 0; i < subcount; i++) {
            for (int j = 0; j < 2; j++) {
                special_shared_vector[i][j] = 0;
            }
        }
        double t1 = clock();

        int forsize_divide_100 = (subcount - minlag) / 100;

        int i_ = 0;
        double fast_prune_cnt = 0;
        double local_prune_cnt = 0;
        double global_prune_cnt = 0;
        double KK_prunes_cnt = 0;
        double P_prunes_cnt = 0;
        double Dtw_prunes_cnt = 0;
        long long all_cnt = 0;

        vector<DOUBLE> del_local(subcount, 0.0);
        vector<DOUBLE> raw_DIFF_UL_local(UTS.size(), 0.0);
        vector<DOUBLE> raw_DIFF_UL2_local(UTS.size(), 0.0);
        vector<DOUBLE> raw_DIFF_UL2_temp_local(UTS.size(), 0.0);
        vector<DOUBLE> DUL2_raw_local(subcount, 0.0);
        vector<DOUBLE> DUL2_local(subcount, 0.0);

        vector<DOUBLE> DUL_local(subcount);
        vector<DOUBLE> sumU_local(subcount, 0.0);
        vector<DOUBLE> sumL_local(subcount, 0.0);
        vector<DOUBLE> sumU2_local(subcount, 0.0);
        vector<DOUBLE> sumL2_local(subcount, 0.0);
        vector<DOUBLE> sumU2_sumL2_loacl(subcount, 0.0);
        vector<DOUBLE> sumU_sumL_loacl(subcount, 0.0);

        vector<DOUBLE> sumMASK_local(subcount, 0.0);
        vector<DOUBLE> norm_U_plus_norm_L_local(subcount, 0.0);

        DOUBLE fast_time = 0;
        DOUBLE global_time = 0;
        DOUBLE local_up_time = 0;
        DOUBLE local_down_time = 0;

        vector<int> diaglist(subcount - minlag);
        int index_diag = 0;

        int alpha_interval = MAX(subseqLen * 5, 1000);
        int buffer_length = MAX(subseqLen / 10, 30);
        int buffer_div_alpha = ceil(1.0 * alpha_interval / buffer_length);
        alpha_interval = buffer_div_alpha * buffer_length;

        int for_cnt = 0;
        double alpha = 0;
        alpha_generater new_alpha;
        double start_pos = -1.5;
        double end_pos = 1.5;
        vector<DOUBLE> f_x;
        vector<DOUBLE> X;
        vector<DOUBLE> bsf_X;
        vector<bool> lb_vector_for_next(subcount);
        vector<double> alpha_list = {-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};

        vector<double> alpha_table(subcount, 0.0);
        
        double best_alpha = 0;
        
        int cb_prune_cnt = 0;
        int cb_temp_cnt = 0;

        for (int i_buffer = 0; i_buffer < buffer_length; i_buffer++)

        {
            int i_table_num = 0;
            for (int diag = minlag + 1 + i_buffer; diag <= subcount; diag += buffer_length)
                
            {
                i_table_num++;
                for_cnt++;

                memset(lb_vector_new.data(), 0, sizeof(double) * subcount);
                fill(lb_vector.begin(), lb_vector.end(), false);

                double cnt_of_purn = 0;
                fast_time++;  

                diag_fast(a, subseqLen, diag, UTS, LTS, mu, sumU_sumL, invsig, norm_U_plus_norm_L_trans, del,
                          lb_vector, dr_bwdU_plus_dr_bwdL, dc_bwd, dr_fwdU_plus_dr_fwdL, dc_fwd,
                          cnt_of_purn, DUL_fast, DUL2_fast, bsf);

                double cnt1 = cnt_of_purn;
                double cnt2;
                fast_prune_cnt += cnt1;
                pruning_rate = cnt_of_purn / double(subcount - diag + 1);

                if(pruning_rate < 0.98) 

                {
                    global_time++;
                    diag_mask_global(a, subseqLen, diag, lb_vector, mu, UTS_global, LTS_global, MASK_global, TS2,
                                     sumMASK_global, invsig, sumU_sumL_global,
                                     dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global, dc_bwd_global,
                                     dc_fwd_global, dr_bwdMASK_global, dr_fwdMASK_global,
                                     dc_bwdTS2_global, dc_fwdTS2_global, DUL2_global, lb_vector_new, bsf,
                                     my_subs,
                                     special_shared_vector, cnt_of_purn, norm_U_plus_norm_L_global, DUL_global);

                }
            else
            {
                LB_KIM(bsf, subseqLen, special_shared_vector, my_subs, subcount - diag + 1, diag, lb_vector,
                       cnt_of_purn);
            }
                cnt2 = cnt_of_purn;
                global_prune_cnt += (cnt2 - cnt1);
                pruning_rate = cnt_of_purn / double(subcount - diag + 1);

                double cnt3 = cnt_of_purn;
                local_prune_cnt += (cnt3 - cnt2);
                double temp_pruning_cnt = 0;
                int best_pruning_cnt = 0;

                vector<double> temp_alpha_list = {best_alpha - 2 * ALPHA_BIAS, best_alpha - ALPHA_BIAS, best_alpha,
                                                  best_alpha + ALPHA_BIAS, best_alpha + 2 * ALPHA_BIAS};

                if (i_buffer == 0) {
                    
                    if (i_table_num % buffer_div_alpha == 0)
                        
                    {
                        for (auto temp_alpha: alpha_list) {
                            temp_pruning_cnt = 0;
                            copy(lb_vector.begin(), lb_vector.end(), lb_vector_for_next.begin());
                            diag_mask_local_up(a, minlag, subseqLen, bsf, diag, lb_vector_for_next,
                                               count_table_local, pos_UU, pos_LL, mu, UTS, LTS,
                                               special_shared_vector, temp_pruning_cnt,
                                               sumU2_sumL2_loacl, sumU_sumL_loacl, sumMASK_local,
                                               norm_U_plus_norm_L_local,
                                               raw_DIFF_UL_local, raw_DIFF_UL2_local, raw_DIFF_UL2_temp_local,
                                               DUL2_raw_local, DUL2_local, TS2, temp_alpha, invsig, invsig_2,
                                               diag - minlag - 1);

                            if (temp_pruning_cnt > best_pruning_cnt) {
                                best_pruning_cnt = temp_pruning_cnt;
                                best_alpha = temp_alpha;
                            }
                        }
                        alpha_table[i_table_num] = best_alpha;
                    } else
                        
                    {
                        if (pruning_rate < 0.85) {
                            for (auto temp_alpha: temp_alpha_list) {
                                temp_pruning_cnt = 0;
                                copy(lb_vector.begin(), lb_vector.end(), lb_vector_for_next.begin());
                                diag_mask_local_up(a, minlag, subseqLen, bsf, diag, lb_vector_for_next,
                                                   count_table_local, pos_UU, pos_LL, mu, UTS, LTS,
                                                   special_shared_vector, temp_pruning_cnt,
                                                   sumU2_sumL2_loacl, sumU_sumL_loacl, sumMASK_local,
                                                   norm_U_plus_norm_L_local,
                                                   raw_DIFF_UL_local, raw_DIFF_UL2_local, raw_DIFF_UL2_temp_local,
                                                   DUL2_raw_local, DUL2_local, TS2, temp_alpha, invsig, invsig_2,
                                                   diag - minlag - 1);

                                if (temp_pruning_cnt > best_pruning_cnt) {
                                    best_pruning_cnt = temp_pruning_cnt;
                                    best_alpha = temp_alpha;
                                }
                            }
                            alpha_table[i_table_num] = best_alpha;
                        } else
                            alpha_table[i_table_num] =
                                    2 * alpha_table[MAX(i_table_num - 1, 0)] - alpha_table[MAX(i_table_num - 2, 0)];

                    }

                }

                alpha = alpha_table[i_table_num] +
                        (alpha_table[i_table_num + 1] - alpha_table[i_table_num]) / buffer_length * (i_buffer);
                
                if (pruning_rate < 0.80) 
                {
                    local_up_time++;
                    int pos_fir = 0;
                    int pos_sec = 0;

                    diag_mask_local_up(a, minlag, subseqLen, bsf, diag, lb_vector,
                                       count_table_local, pos_UU, pos_LL, mu, UTS, LTS,
                                       special_shared_vector, cnt_of_purn,
                                       sumU2_sumL2_loacl, sumU_sumL_loacl, sumMASK_local, norm_U_plus_norm_L_local,
                                       raw_DIFF_UL_local, raw_DIFF_UL2_local, raw_DIFF_UL2_temp_local, DUL2_raw_local,
                                       DUL2_local, TS2, 0, invsig, invsig_2, diag - minlag - 1);

                }

                pruning_rate = cnt_of_purn / double(subcount - diag + 1);

                if (pruning_rate < 0.80) 
                {

                    local_down_time++;
                    diag_mask_local_down(a, minlag, subseqLen, bsf, diag, lb_vector,
                                         count_table_local, pos_UU, pos_LL, mu, UTS, LTS,
                                         special_shared_vector, cnt_of_purn,
                                         sumU2_sumL2_loacl, sumU_sumL_loacl, sumMASK_local, norm_U_plus_norm_L_local,
                                         raw_DIFF_UL_local, raw_DIFF_UL2_local, raw_DIFF_UL2_temp_local, DUL2_raw_local,
                                         DUL2_local, TS2, 0, invsig, invsig_2, diag - minlag - 1);

                }

                {
                    for (i = 1; i <= subcount - diag + 1; i++) {
                        all_cnt++;
                        if (!lb_vector[i - 1]) {

                            KK_prunes_cnt++;
                            if (LB_KK(my_subs[i - 1], subs_U[i + diag - 2].data(), subs_L[i + diag - 2].data(),
                                      subseqLen, bsf * bsf, cb1, special_shared_vector[i - 1], lbk1)) {
                                if (LB_KK(my_subs[i + diag - 2], subs_U[i - 1].data(), subs_L[i - 1].data(), subseqLen,
                                          bsf * bsf, cb2,
                                          special_shared_vector[i - 1], lbk2)) {
                                    KK_prunes_cnt--;
                                    P_prunes_cnt++;

                                    {
                                if(lbPetitjean_new(proj, my_subs[i + diag - 2], UTS_p, LTS_p, subs_U[i + diag - 2].data(), subs_L[i + diag - 2].data(),
                                                   my_subs[i - 1], subs_U[i - 1].data(), subs_L[i - 1].data(), warpmax, subseqLen,
                                                   bsf*bsf, lbk2,
                                                   cb2, mu[i - 1], sig[i - 1]))
                                        {
                                            P_prunes_cnt--;
                                            Dtw_prunes_cnt++;

                                            int m = subseqLen;
                                            if (lbk1 > lbk2) {
                                                cb[m - 1] = cb1[m - 1];
                                                for (int k = m - 2; k >= 0; k--)
                                                    cb[k] = cb[k + 1] + cb1[k];
                                            } else {
                                                cb[m - 1] = cb2[m - 1];
                                                for (int k = m - 2; k >= 0; k--)
                                                    cb[k] = cb[k + 1] + cb2[k];
                                            }

                                            cb_prune_cnt += cb_temp_cnt;
                                            double dtw_distance = sqrt(
                                                    MON_dtw(my_subs[i - 1], my_subs[i + diag - 2], cb, subseqLen,
                                                            warpmax, bsf * bsf));
    
                                            if (dtw_distance <= bsf) {
                                                first_min = i;
                                                sec_min = i + diag - 1;
                                                bsf = dtw_distance;

                                            }
                                        }

                                    }
                                }
                            }

                        }
                    }
                }

            }
        }

        cout << endl;
        cout << "cb  prunes " << cb_prune_cnt / (Dtw_prunes_cnt) << endl;

        cout << "fast times :" << fast_time << endl;
        cout << "global times :" << global_time << endl;
        cout << "local_up times :" << local_up_time << endl;
        cout << "local_down times :" << local_down_time << endl;
        cout << "DTW times " << (Dtw_prunes_cnt) << endl;

        cout << "fast prunes " << (fast_prune_cnt) / all_cnt << endl;
        cout << "global prunes " << (global_prune_cnt) / all_cnt << endl;
        cout << "local up prunes " << cnt_up / all_cnt << endl;
        cout << "local down prunes " << cnt_down / all_cnt << endl;
        cout << "KK prunes " << (KK_prunes_cnt) / all_cnt << endl;
        cout << "Pe prunes " << (P_prunes_cnt) / all_cnt << endl;
        cout << "DTW prunes " << (Dtw_prunes_cnt) / all_cnt << endl;
        cout << "first min is " << first_min << " " << "sec min is " << sec_min << endl;

        double t2 = clock();
        printf("total time is %fs   input of bsf:%5.3f  ,final distance :%f  \n", DOUBLE(t2 - t1) / CLOCKS_PER_SEC,
               best_so_far, bsf);

        FILE *op = fopen(output_time_file_path,"a");
        fprintf(op, "%lf", DOUBLE(t2 - t1)/CLOCKS_PER_SEC);
        fprintf(op, "\n");
        fclose(op);

    }

