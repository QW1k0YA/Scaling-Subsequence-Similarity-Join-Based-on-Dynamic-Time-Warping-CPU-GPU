
#include <cstring>
#include <iomanip>
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../alldef/typedefdouble.cuh"
#include "../allunder/underdtw.cuh"
#include "warp_num_per_block.h"
#include "collect_indices_kernal_after_keogh.h"
#include "process_keogh_and_dtw_kernal.cuh"
#include "device_prefix_sum.cuh"
#include <cuda_profiler_api.h>
#define DOUBLE_BIAS 6
#define BIAS 3

using namespace  std;
__device__ cudaStream_t stream_pool[STREAM_POOL_SIZE_DTW];
__global__ void init_stream_pool(cudaStream_t* streams) {
    int tid = threadIdx.x;
    if (tid < STREAM_POOL_SIZE_DTW) {
        cudaStreamCreateWithFlags(&streams[tid], cudaStreamNonBlocking);
    }
}

void new_dtw_motifGUI_malloc(const vector<FLOAT> &a, int subseqLen, int maxwarp, const vector<FLOAT> &mp_ed,
                             const vector<int> &mp_ed_index, const char *file)
{
    int subcount = a.size() - subseqLen + 1;
    int len = a.size();

    cudaSetDevice(0);
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp prop;
    cudaError_t errO = cudaGetDeviceProperties(&prop, device);
    int sm_count = prop.multiProcessorCount;

    cout << "sm = " << sm_count << endl;
    
    int THREADS_NUM = GRID_SIZE*BLOCK_SIZE;

    vector<FLOAT > mu(subcount),sig(subcount);
    mvmean(a,subseqLen,mu,sig);
    vector<FLOAT > invsig(subcount);
    vector<FLOAT > invsig_2(subcount);
    for(int i = 0;i < subcount;i ++)
    {
        invsig[i] = 1/sig[i];
        invsig_2[i] = invsig[i]*invsig[i];
    }

    FLOAT ** my_subs = (FLOAT **)malloc(subcount * sizeof(FLOAT *));
    FLOAT ** subs_U = (FLOAT **)malloc(subcount * sizeof(FLOAT *));
    FLOAT ** subs_L = (FLOAT **)malloc(subcount * sizeof(FLOAT *));

    for (int i = 0; i < subcount; ++i) {
        
        my_subs[i] = (FLOAT *)calloc(subseqLen, sizeof(FLOAT ));  
        subs_U[i] = (FLOAT *)calloc(subseqLen, sizeof(FLOAT ));
        subs_L[i] = (FLOAT *)calloc(subseqLen, sizeof(FLOAT ));
    }
    size_t temp_1;
    for(int i = 0;i < subcount;i ++)
    {
        for(int j = 0;j < subseqLen;j++)
        {
            my_subs[i][j] = (a[i + j] - mu[i]) * invsig[i];
        }
        lower_upper_lemire(my_subs[i],subseqLen,maxwarp,subs_L[i],subs_U[i]);
    }

    auto my_total_time = std::chrono::high_resolution_clock::now();
    printf("buffer_length: %d, maxwarp: %d\n",subseqLen,maxwarp);

    int minspacing = subseqLen;

    int k = 100;
    vector<int> min_index = min_v_k(mp_ed,k);

    FLOAT best_so_far = mp_ed[min_index[0]];
    std::vector<FLOAT> cb_temp(subseqLen + maxwarp + 1, 0);

    int low_index;
    int high_index;
    for(int i = 0;i < k;i++)
    {
        low_index = min_index[i];
        high_index = mp_ed_index[min_index[i]] - 1;

        if(LB_KK_FIRST(my_subs[low_index], my_subs[high_index], subs_U[high_index], subs_L[low_index], subseqLen,
                       best_so_far * best_so_far, cb_temp))
        {
            FLOAT dtw_value = sqrt(MON_dtw_host(my_subs[low_index], my_subs[high_index], cb_temp.data(),
                                                 subseqLen, maxwarp, best_so_far * best_so_far));

            if(dtw_value < best_so_far)
            {
                best_so_far = dtw_value;
            }
        }

    }

    cout << "best so far is " << best_so_far << endl;
    cout << "the buffer_length of TS is " <<  a.size() << endl;
    printf("################Start calculating lower bound#################  \n");

    int debug_sum = 0;
    FLOAT  bsf = best_so_far;

    bool * lb_vector= (bool*) malloc(subcount*sizeof(FLOAT ));
    
    subcount = a.size() - subseqLen + 1;

    int warpmax = maxwarp;

    vector<FLOAT > UTS(len);
    vector<FLOAT > LTS(len);
    lower_upper_lemire(a,len,warpmax,LTS,UTS);

    vector<FLOAT > UTS_p(subseqLen);
    vector<FLOAT > LTS_p(subseqLen);

    vector<FLOAT > dr_bwdU, dr_bwdL, dc_bwd;
    dr_bwdU = addElementToFront(extr_vfromv(UTS, 1, subcount - 1), 0.0);
    dr_bwdL = addElementToFront(extr_vfromv(LTS, 1, subcount - 1), 0.0);
    dc_bwd = addElementToFront(extr_vfromv(a, 1, subcount - 1), 0.0);

    vector<FLOAT > dr_fwdU, dr_fwdL, dc_fwd;
    dr_fwdU = extr_vfromv(UTS, subseqLen, len);
    dr_fwdL = extr_vfromv(LTS, subseqLen, len);
    dc_fwd = extr_vfromv(a, subseqLen, len);

    vector<FLOAT > dr_bwdU_plus_dr_bwdL, dr_fwdU_plus_dr_fwdL;
    dr_bwdU_plus_dr_bwdL = plusvector(dr_bwdU, dr_bwdL);
    dr_fwdU_plus_dr_fwdL = plusvector(dr_fwdU, dr_fwdL);

    const  vector<FLOAT >sumU = movsum(UTS, subseqLen - 1);
    const  vector<FLOAT >sumL = movsum(LTS, subseqLen - 1);
    const  vector<FLOAT >sumU2 = movsum(elementWiseMultiply(UTS, UTS), subseqLen - 1);
    const  vector<FLOAT >sumL2 = movsum(elementWiseMultiply(LTS, LTS), subseqLen - 1);
    const  vector<FLOAT >sumU2_sumL2 = plusvector(sumU2, sumL2);
    const  vector <FLOAT> sumU_sumL = plusvector(sumU, sumL);

    vector<FLOAT > del(subcount, 0);

    vector<FLOAT > normLTS(subseqLen, 0.0), normUTS(subseqLen, 0.0);
    vector<FLOAT > DUL2_fast(subcount);
    vector<FLOAT > DUL_fast(subcount);

    for (int row = 0; row < subcount; row++) {
        for (int i = row; i <= row + subseqLen - 1; i++) {
            normLTS[i - row] = (LTS[i] - mu[row])*invsig[row];
            normLTS[i - row] = (LTS[i] - mu[row])*invsig[row];
            normUTS[i - row] = (UTS[i] - mu[row])*invsig[row];
        }

        DUL2_fast[row] = pow(norm_vector(normLTS, 2), 2) + pow(norm_vector(normUTS, 2), 2) -
                         2 * sum_vector((elementWiseMultiply(normLTS, normUTS)));
        DUL_fast[row] = pow(MAX(DUL2_fast[row], 0.0), 0.5);

    }

    vector<FLOAT > norm_U_plus_norm_L_trans(subcount, 0.0);

    for (int row = 0; row < subcount; row++) {
        norm_U_plus_norm_L_trans[row] = (sumU2_sumL2[row] - 2 * sumU_sumL[row] * mu[row]) * invsig[row];
    }

    int len_of_table_local=INITIAL_LEN_OF_TABLE;

    FLOAT  real_max,real_min;
    vector<FLOAT > pos_UU(len, 1.0);
    vector<FLOAT > pos_LL(len, 1.0);
    int len_of_cdf;
     real_max = my_subs[0][0];
    real_min = my_subs[0][0];

    for(int row = 0;row < subcount;row++){
        for(int col = 0;col < subseqLen;col++){
            if(my_subs[row][col]> real_max && my_subs[row][col] < 100000) real_max = my_subs[row][col];
            if(my_subs[row][col] < real_min  && my_subs[row][col] > -100000) real_min = my_subs[row][col];
        }
    }

    FLOAT MAX_REAL_VALUE = MAX(abs(real_max), abs(real_min));
    if(MAX_REAL_VALUE > 3){
        len_of_cdf = (MAX_REAL_VALUE - 3)*20+600;
    }
    else{len_of_cdf = 600;
    }

    FLOAT  cov_U_plus_cov_L;

    FLOAT  pruning_rate = 0;
    FLOAT  thr = 0.9;

    int minlag = subseqLen;

    vector<FLOAT > proj(subseqLen);
    FLOAT lbk1,lbk2;

    vector<FLOAT > row_min_DEL(subcount, INFINITY);

    FLOAT  max_real_value = floor(MAX(abs(real_max), abs(real_min))) * 2 * 100;
    auto len_of_table_global = static_cast<long long>(MAX(INITIAL_LEN_OF_TABLE, max_real_value));

    vector<vector<FLOAT >> count_table_global(subseqLen, vector<FLOAT >(len_of_table_global, 0.0));

    printf("LEN_OF_TABLE_GLOBAL is %lld \n", len_of_table_global);
    compute_shared_data_global(a, subseqLen,my_subs,len_of_table_global,count_table_global);

    vector<FLOAT > MASK_global(len, 0.0);
    FLOAT  qqq,www;
    vector<FLOAT > TS2 = elementWiseMultiply(a, a);
    compute_MASK_global(a, subseqLen, len_of_table_global, count_table_global, subs_U, subs_L, MASK_global);

    int i;
    vector<FLOAT > UTS_global = elementWiseMultiply(UTS, MASK_global);
    vector<FLOAT > LTS_global = elementWiseMultiply(LTS, MASK_global);

    vector<FLOAT > dr_bwdU_global = addElementToFront(extr_vfromv(UTS_global, 1 + BIAS, UTS.size() - 1), 0.0);
    vector<FLOAT > dr_bwdL_global = addElementToFront(extr_vfromv(LTS_global, 1 + BIAS, LTS.size() - 1), 0.0);

    vector<FLOAT > dr_fwdU_global = extr_vfromv(UTS_global, subseqLen - BIAS, UTS_global.size());
    vector<FLOAT > dr_fwdL_global = extr_vfromv(LTS_global, subseqLen - BIAS, LTS_global.size());

    vector<FLOAT >dc_bwd_global = addElementToFront(extr_vfromv(a, 1 + BIAS, len), 0.0);
    vector<FLOAT >dc_fwd_global = extr_vfromv(a, subseqLen - BIAS, len);

    vector<FLOAT > dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global;
    dr_bwdU_plus_dr_bwdL_global = plusvector(dr_bwdU_global, dr_bwdL_global);
    dr_fwdU_plus_dr_fwdL_global = plusvector(dr_fwdU_global, dr_fwdL_global);

    const vector<FLOAT >sumU_global = movsum_p(UTS_global.data() + BIAS, len - DOUBLE_BIAS , subseqLen - 1 - DOUBLE_BIAS);
    const vector<FLOAT >sumL_global = movsum_p(LTS_global.data() + BIAS, len - DOUBLE_BIAS, subseqLen - 1 - DOUBLE_BIAS);
    const vector<FLOAT >sumU2_global = movsum(elementWiseMultiply_p(UTS_global.data() + BIAS, UTS_global.data() + BIAS, len - DOUBLE_BIAS), subseqLen - 1 - DOUBLE_BIAS);
    const vector<FLOAT >sumL2_global = movsum(elementWiseMultiply_p(LTS_global.data() + BIAS, LTS_global.data() + BIAS, len - DOUBLE_BIAS), subseqLen - 1 - DOUBLE_BIAS);
    const vector<FLOAT >sumU2_sumL2_global = plusvector(sumU2_global, sumL2_global);
    const vector<FLOAT >sumU_sumL_global = plusvector(sumU_global, sumL_global);

    const vector<FLOAT > sumMASK_global = movsum_p(MASK_global.data() + BIAS, MASK_global.size() - DOUBLE_BIAS, subseqLen - 1 - DOUBLE_BIAS);
    vector<FLOAT > norm_U_plus_norm_L_global(subcount, 0.0);

    const vector<FLOAT > TS2_global = elementWiseMultiply(a, a);

    for(int row = 1;row <= subcount;row++)
    {
        norm_U_plus_norm_L_global[row - 1] = (sumU2_sumL2_global[row - 1] - 2 * sumU_sumL_global[row - 1] * mu[row - 1]
                                             + 2*sumMASK_global[row-1]*mu[row-1]*mu[row-1]) * invsig[row-1] * invsig[row-1];
    }

    vector<FLOAT > del_global(subcount, 0);
    vector<FLOAT > normLTS_global(subseqLen, 0.0), normUTS_global(subseqLen, 0.0);
    FLOAT  del_ths_global;
    vector<FLOAT  >DUL2_global(subcount);
    vector<FLOAT  >DUL_global(subcount);

    for(int row = 0;row < subcount;row++)
    {
        for (i = row; i <= row + subseqLen - 1; i++) {

            normLTS_global[i - row] = (LTS_global[i] - mu[row]) *invsig[row];
            normUTS_global[i - row] = (UTS_global[i] - mu[row]) *invsig[row];
        }

        DUL2_global[row] =pow(norm_vector(substractvector(normLTS_global,normUTS_global),2),2);
        
        DUL_global[row] = sqrt(MAX(DUL2_global[row],0));

    }

    const vector<FLOAT > dr_bwdMASK_global = addElementToFront(extr_vfromv(MASK_global, 1 + BIAS, MASK_global.size()), 0.0);
    const vector<FLOAT > dr_fwdMASK_global = extr_vfromv(MASK_global, subseqLen - BIAS,
                                                        MASK_global.size());

    const vector<FLOAT > dc_bwdTS2_global = addElementToFront(extr_vfromv(TS2, 1 + BIAS , TS2.size() - 1), 0.0);
    const vector<FLOAT > dc_fwdTS2_global = extr_vfromv(TS2, subseqLen - BIAS, TS2.size());

    vector<FLOAT > lb_vector_new(subcount);
    FLOAT **special_shared_vector = (FLOAT **)malloc(subcount * sizeof(FLOAT *));
    
    for (i = 0; i < subcount; i++) {
        special_shared_vector[i] = (FLOAT *)malloc(2 * sizeof(FLOAT));
        if (special_shared_vector[i] == nullptr) {

            for (int j = 0; j < i; j++) {
                free(special_shared_vector[j]);
            }
            free(special_shared_vector);
        }
    }
    for(i = 0;i < subcount;i++)
    {
        for(int j = 0;j < 2;j++)
        {
            special_shared_vector[i][j] = 0;
        }
    }

    int forsize_divide_100 = (subcount- minlag)/100 ;

    int i_ = 0;
    FLOAT fast_prune_cnt = 0;
    FLOAT local_prune_cnt = 0;
    FLOAT global_prune_cnt = 0;
    FLOAT KK_prunes_cnt = 0;
    FLOAT P_prunes_cnt = 0;
    FLOAT Dtw_prunes_cnt = 0;
    long long all_cnt = 0;

    vector<FLOAT > del_local(subcount, 0.0);
    vector<FLOAT > raw_DIFF_UL_local(UTS.size(), 0.0);
    vector<FLOAT > raw_DIFF_UL2_local (UTS.size(), 0.0);
    vector<FLOAT > raw_DIFF_UL2_temp_local (UTS.size(), 0.0);
    vector<FLOAT > DUL2_raw_local(subcount, 0.0);
    vector<FLOAT > DUL2_local(subcount, 0.0);

    vector<FLOAT > DUL_local(subcount);
    vector<FLOAT >sumU_local (subcount, 0.0);
    vector<FLOAT >sumL_local (subcount, 0.0);
    vector<FLOAT >sumU2_local  (subcount, 0.0);
    vector<FLOAT >sumL2_local (subcount, 0.0);
    vector<FLOAT >sumU2_sumL2_loacl (subcount, 0.0);
    vector<FLOAT >sumU_sumL_loacl (subcount, 0.0);

    vector<FLOAT > sumMASK_local(subcount, 0.0);
    vector<FLOAT > norm_U_plus_norm_L_local(subcount, 0.0);

    FLOAT  fast_time = 0;
    FLOAT  global_time = 0;
    FLOAT  local_up_time = 0;
    FLOAT  local_down_time = 0;

    int index_diag = 0;
    int alpha_interval = MAX(subseqLen*5,1000);
    int buffer_length = MAX(subseqLen/10,30);
    int buffer_div_alpha = ceil(1.0*alpha_interval/buffer_length);
    alpha_interval = buffer_div_alpha*buffer_length;

    int for_cnt = 0;
    FLOAT alpha = 0;
    FLOAT start_pos = -1.5;
    FLOAT end_pos = 1.5;
    vector<FLOAT > f_x;
    vector<FLOAT > X;
    vector<FLOAT > bsf_X;
    vector<bool> lb_vector_for_next(subcount);

    vector<FLOAT> alpha_table(subcount,0.0);
    FLOAT best_alpha = 0;

    FLOAT  *d_a, *d_TS2,*d_mu, *d_sig, *d_sumU_sumL, *d_invsig, *d_norm_U_plus_norm_L_trans, *d_del;
    FLOAT  *d_dr_bwdU_plus_dr_bwdL, *d_dc_bwd, *d_dr_fwdU_plus_dr_fwdL, *d_dc_fwd;
    FLOAT  *d_UTS, *d_LTS, *d_UTS_global, *d_LTS_global, *d_MASK_global;
    FLOAT  *d_pos_UU, *d_pos_LL;
    SHORT *d_count_table_local;
    
    FLOAT *d_sumMASK_global, *d_sumU_sumL_global;
    FLOAT *d_dr_bwdU_plus_dr_bwdL_global, *d_dr_fwdU_plus_dr_fwdL_global;
    FLOAT *d_dc_bwd_global, *d_dc_fwd_global;
    FLOAT *d_dr_bwdMASK_global, *d_dr_fwdMASK_global;
    FLOAT *d_dc_bwdTS2_global, *d_dc_fwdTS2_global;
    FLOAT *d_DUL2_global, *d_norm_U_plus_norm_L_global, *d_DUL_global;
    FLOAT *d_DUL_fast,*d_DUL2_fast;
    FLOAT  *d_bsf;
    FLOAT *d_cnt_up, *d_cnt_down;

    int *d_prefix_sum;
    cudaMalloc(&d_prefix_sum, GRID_SIZE * sizeof(int)); CUERR
    int *d_num_inclusive;
    cudaMalloc(&d_num_inclusive, THREADS_NUM * sizeof(int)); CUERR

    int init_d_num_of_dtw = 0;

    size_t *h_num_of_dtw;
    h_num_of_dtw = (size_t*) malloc(sizeof(size_t));
    h_num_of_dtw[0] = 0;

    cudaMalloc(&d_a, a.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_a, a.data(), a.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_TS2, TS2.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_TS2, TS2.data(), TS2.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_mu, mu.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_mu, mu.data(), mu.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_sig, sig.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_sig, sig.data(), sig.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_sumU_sumL, sumU_sumL.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_sumU_sumL, sumU_sumL.data(), sumU_sumL.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_invsig, invsig.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_invsig, invsig.data(), invsig.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_norm_U_plus_norm_L_trans, norm_U_plus_norm_L_trans.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_norm_U_plus_norm_L_trans, norm_U_plus_norm_L_trans.data(),
               norm_U_plus_norm_L_trans.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dr_bwdU_plus_dr_bwdL, dr_bwdU_plus_dr_bwdL.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_dr_bwdU_plus_dr_bwdL, dr_bwdU_plus_dr_bwdL.data(),
               dr_bwdU_plus_dr_bwdL.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dc_bwd, dc_bwd.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_dc_bwd, dc_bwd.data(), dc_bwd.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dr_fwdU_plus_dr_fwdL, dr_fwdU_plus_dr_fwdL.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_dr_fwdU_plus_dr_fwdL, dr_fwdU_plus_dr_fwdL.data(),
               dr_fwdU_plus_dr_fwdL.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dc_fwd, dc_fwd.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_dc_fwd, dc_fwd.data(), dc_fwd.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_UTS, UTS.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_UTS, UTS.data(), UTS.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_LTS, LTS.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_LTS, LTS.data(), LTS.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_UTS_global, UTS_global.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_UTS_global, UTS_global.data(), UTS_global.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_LTS_global, LTS_global.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_LTS_global, LTS_global.data(), LTS_global.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_MASK_global, MASK_global.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_MASK_global, MASK_global.data(), MASK_global.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_pos_UU, pos_UU.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_pos_UU, pos_UU.data(), pos_UU.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_pos_LL, pos_LL.size() * sizeof(FLOAT )); CUERR
    cudaMemcpy(d_pos_LL, pos_LL.data(), pos_LL.size() * sizeof(FLOAT ), cudaMemcpyHostToDevice);

    cudaMalloc(&d_sumMASK_global, sumMASK_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_sumMASK_global, sumMASK_global.data(),
               sumMASK_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_sumU_sumL_global, sumU_sumL_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_sumU_sumL_global, sumU_sumL_global.data(),
               sumU_sumL_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dr_bwdU_plus_dr_bwdL_global, dr_bwdU_plus_dr_bwdL_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dr_bwdU_plus_dr_bwdL_global, dr_bwdU_plus_dr_bwdL_global.data(),
               dr_bwdU_plus_dr_bwdL_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dr_fwdU_plus_dr_fwdL_global, dr_fwdU_plus_dr_fwdL_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dr_fwdU_plus_dr_fwdL_global, dr_fwdU_plus_dr_fwdL_global.data(),
               dr_fwdU_plus_dr_fwdL_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dc_bwd_global, dc_bwd_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dc_bwd_global, dc_bwd_global.data(),
               dc_bwd_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dc_fwd_global, dc_fwd_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dc_fwd_global, dc_fwd_global.data(),
               dc_fwd_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dr_bwdMASK_global, dr_bwdMASK_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dr_bwdMASK_global, dr_bwdMASK_global.data(),
               dr_bwdMASK_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dr_fwdMASK_global, dr_fwdMASK_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dr_fwdMASK_global, dr_fwdMASK_global.data(),
               dr_fwdMASK_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dc_bwdTS2_global, dc_bwdTS2_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dc_bwdTS2_global, dc_bwdTS2_global.data(),
               dc_bwdTS2_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_dc_fwdTS2_global, dc_fwdTS2_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_dc_fwdTS2_global, dc_fwdTS2_global.data(),
               dc_fwdTS2_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_DUL2_global, DUL2_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_DUL2_global, DUL2_global.data(),
               DUL2_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_norm_U_plus_norm_L_global, norm_U_plus_norm_L_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_norm_U_plus_norm_L_global, norm_U_plus_norm_L_global.data(),
               norm_U_plus_norm_L_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_DUL_global, DUL_global.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_DUL_global, DUL_global.data(),
               DUL_global.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    cudaMalloc(&d_DUL_fast, DUL_fast.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_DUL_fast, DUL_fast.data(),
               DUL_fast.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);
    cudaMalloc(&d_DUL2_fast, DUL2_fast.size() * sizeof(FLOAT)); CUERR
    cudaMemcpy(d_DUL2_fast, DUL2_fast.data(),
               DUL2_fast.size() * sizeof(FLOAT), cudaMemcpyHostToDevice);

    int num_threads = THREADS_NUM;

    bool* d_lb_vector;
    cudaMalloc(&d_lb_vector, num_threads * STEP_LENGTH * sizeof(bool));
    bool* h_lb_vector = new bool[num_threads * STEP_LENGTH];
    memset(h_lb_vector,0,sizeof(bool)*num_threads * STEP_LENGTH);
    cudaMemcpy(d_lb_vector, h_lb_vector, num_threads * STEP_LENGTH * sizeof(bool), cudaMemcpyHostToDevice);

    int* d_indices;
    cudaMalloc(&d_indices, num_threads * STEP_LENGTH * sizeof(int)); CUERR
    int* d_diag;
    cudaMalloc(&d_diag, num_threads * STEP_LENGTH * sizeof(int)); CUERR

    int* h_indices = new int[num_threads * STEP_LENGTH];
    memset( h_indices,0,sizeof(int)*num_threads * STEP_LENGTH);
    cudaMemcpy(d_diag,  h_indices, num_threads * STEP_LENGTH * sizeof(int), cudaMemcpyHostToDevice);
    int* h_indices2 = new int[num_threads * STEP_LENGTH];
    memset( h_indices2,0,sizeof(int)*num_threads * STEP_LENGTH);
    cudaMemcpy(d_indices,  h_indices, num_threads * STEP_LENGTH * sizeof(int), cudaMemcpyHostToDevice);

    FLOAT*  d_lb_vector_new;
    cudaMalloc(&d_lb_vector_new, num_threads * STEP_LENGTH * sizeof(FLOAT));

    cout << "best so far is " << bsf << endl;
    FLOAT* h_bsf = (FLOAT *)malloc(BSF_POOL * sizeof(FLOAT));
    for(int i = 0; i < BSF_POOL; i++)
    {
        h_bsf[i] = bsf;
    }
    FLOAT*  d_bsf_global = nullptr;

    cudaMalloc(&d_bsf_global, BSF_POOL * sizeof(FLOAT));
    cudaMemcpy(d_bsf_global, h_bsf, BSF_POOL * sizeof(FLOAT), cudaMemcpyHostToDevice);

    FLOAT ** d_my_subs = nullptr;
    FLOAT ** d_subs_U = nullptr;
    FLOAT ** d_subs_L = nullptr;
    FLOAT ** d_shared_special_vector = nullptr;

    cudaMalloc(&d_my_subs, subcount * sizeof(FLOAT *)); CUERR
    cudaMalloc(&d_subs_U, subcount * sizeof(FLOAT *)); CUERR
    cudaMalloc(&d_subs_L, subcount * sizeof(FLOAT *)); CUERR
    cudaMalloc(&d_shared_special_vector, subcount *sizeof(FLOAT *)); CUERR

    FLOAT ** tmp_my_subs = (FLOAT **)malloc(subcount * sizeof(FLOAT *));
    FLOAT ** tmp_subs_U = (FLOAT **)malloc(subcount * sizeof(FLOAT *));
    FLOAT ** tmp_subs_L = (FLOAT **)malloc(subcount * sizeof(FLOAT *));
    FLOAT ** tmp_d_shared_special_vector = (FLOAT **)malloc(subcount * sizeof(FLOAT *));

    for (int i = 0; i < subcount; ++i) {
        
        cudaMalloc(&tmp_my_subs[i], subseqLen * sizeof(FLOAT )); CUERR
        cudaMalloc(&tmp_subs_U[i], subseqLen * sizeof(FLOAT )); CUERR
        cudaMalloc(&tmp_subs_L[i], subseqLen * sizeof(FLOAT )); CUERR
        cudaMalloc(&tmp_d_shared_special_vector[i],  2 * sizeof(FLOAT )); CUERR

        cudaMemcpy(tmp_my_subs[i], my_subs[i],
                   subseqLen * sizeof(FLOAT ), cudaMemcpyHostToDevice);
        cudaMemcpy(tmp_subs_U[i], subs_U[i],
                   subseqLen * sizeof(FLOAT ), cudaMemcpyHostToDevice);
        cudaMemcpy(tmp_subs_L[i], subs_L[i],
                   subseqLen * sizeof(FLOAT ), cudaMemcpyHostToDevice);
        cudaMemcpy(tmp_d_shared_special_vector[i], special_shared_vector[i],
                   2 * sizeof(FLOAT ), cudaMemcpyHostToDevice);
    }

    cudaMemcpy(d_my_subs, tmp_my_subs,
               subcount * sizeof(FLOAT *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_subs_U, tmp_subs_U,
               subcount * sizeof(FLOAT *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_subs_L, tmp_subs_L,
               subcount * sizeof(FLOAT *), cudaMemcpyHostToDevice);
    cudaMemcpy(d_shared_special_vector, tmp_d_shared_special_vector,
               2 * sizeof(FLOAT *), cudaMemcpyHostToDevice);

    for (int i = 0; i < subcount; ++i) {
        free(my_subs[i]);
        free(subs_U[i]);
        free(subs_L[i]);
    }

    const int block_size = BLOCK_SIZE;       
    const int num_blocks = GRID_SIZE; 

    cudaError_t error = cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024 * 1024);
    if (error != cudaSuccess) {

    }

    int step = STEP_LENGTH;

    if (errO != cudaSuccess) {
        std::cerr << "Failed to get device properties: " << cudaGetErrorString(errO) << std::endl;
    }
    int numBlocks;
    int blockSize = 32*WARP_NUMS;

    cudaMemcpy(&bsf, d_bsf_global, sizeof(FLOAT), cudaMemcpyDeviceToHost);

    cout << "initial bsf = " << bsf << endl;

    FLOAT t_first = clock();
    double t_total_dtw = 0;
    double t_total_diag = 0;
    int d_start_pos,d_end_pos;
    int t_for_cnt = ceil(1.0*subcount/step);

    size_t dtw_cnt = 0;
    int DTW_SHARED_MEM_SIZE = 2*subseqLen*sizeof(FLOAT)*WARP_NUMS;

    bool profileflag = true;

    int bl_size = ceil(maxwarp/31.0);
    if(bl_size == 2)
    {
        bl_size = 3;
    }
    else if(bl_size == 4)
    {
        bl_size = 5;
    }
    for(int t_index = 0;t_index < t_for_cnt;t_index ++)
    {
        d_start_pos = t_index * step;
        d_end_pos = d_start_pos + step;

        for (int diag = minlag + 1; diag <= subcount; diag += THREADS_NUM)
        {
            if (d_start_pos > subcount - diag + 1) {
                continue;
            }

            FLOAT t1 = clock();

            GLOBAL_DIAG<<<num_blocks, block_size>>>(
                    minlag,
                    subcount,
                    subseqLen,
                    len,
                    warpmax,
                    d_a,          
                    d_mu,
                    d_sumU_sumL,
                    d_invsig,
                    d_norm_U_plus_norm_L_trans,
                    d_dr_bwdU_plus_dr_bwdL,
                    d_dc_bwd,
                    d_dr_fwdU_plus_dr_fwdL,
                    d_dc_fwd,
                    d_UTS,
                    d_LTS,
                    d_UTS_global,
                    d_LTS_global,
                    d_MASK_global,
                    d_TS2,
                    d_lb_vector,
                    d_lb_vector_new,
                    d_sumMASK_global,
                    d_sumU_sumL_global,
                    d_dr_bwdU_plus_dr_bwdL_global,
                    d_dr_fwdU_plus_dr_fwdL_global,
                    d_dc_bwd_global,
                    d_dc_fwd_global,
                    d_dr_bwdMASK_global,
                    d_dr_fwdMASK_global,
                    d_dc_bwdTS2_global,
                    d_dc_fwdTS2_global,
                    d_DUL2_global,
                    d_DUL_global,
                    d_norm_U_plus_norm_L_global,
                    d_my_subs,
                    d_shared_special_vector,
                    d_bsf_global, d_DUL_fast,
                    d_DUL2_fast, diag, d_start_pos, d_end_pos);

            cudaError_t syncErr = cudaDeviceSynchronize();
            if (syncErr != cudaSuccess) {
                printf("fail in diag %s\n", cudaGetErrorString(syncErr));
                printf("diag = %d st = %d end = %d\n",diag,d_start_pos,d_end_pos);
            }

            cudaMemcpy(d_prefix_sum, &init_d_num_of_dtw, sizeof(int), cudaMemcpyHostToDevice);

            FLOAT t2 = clock();
            
            calculate_num_for_each_block<<<GRID_SIZE, BLOCK_SIZE>>>(d_lb_vector, subcount,
                                                                    d_start_pos, d_end_pos,
                                                                    diag, d_prefix_sum, d_num_inclusive);
            cudaDeviceSynchronize();

            device_prefix_sum_block<<<1, GRID_SIZE>>>(d_prefix_sum);

            cudaDeviceSynchronize();

            collect_indices<<<GRID_SIZE, BLOCK_SIZE>>>(d_diag, d_indices, d_lb_vector, subcount,
                                                       d_start_pos, d_end_pos,
                                                       diag, d_prefix_sum, d_num_inclusive);

            cudaMemcpy(h_num_of_dtw, &d_prefix_sum[GRID_SIZE - 1], sizeof(int), cudaMemcpyDeviceToHost);
            
            dtw_cnt+= h_num_of_dtw[0];
            syncErr = cudaDeviceSynchronize();
            if (syncErr != cudaSuccess) {
                printf("fail in collect %s\n", cudaGetErrorString(syncErr));
            }

            FLOAT t3 = clock();

            if(subseqLen < 1024)
            {

                process_keogh_and_dtw_kernel_for_a_Parallelogram<<<h_num_of_dtw[0], 32, DTW_SHARED_MEM_SIZE>>>
                        (d_my_subs, d_subs_L, d_subs_U, subseqLen, d_bsf_global, subcount,
                         warpmax, d_indices, d_diag);

            }
            else
            {

                process_keogh_and_dtw_kernel_for_a_Parallelogram_without_shared_memory_and_nomalized<<<h_num_of_dtw[0], 32>>>
                        (d_my_subs, d_UTS, d_LTS, d_mu, d_invsig, subseqLen, d_bsf_global, subcount,
                         warpmax, d_indices, d_diag, bl_size);

            }
            syncErr = cudaDeviceSynchronize();
            if (syncErr != cudaSuccess) {
                printf("fail in process %s\n", cudaGetErrorString(syncErr));
                printf("diag = %d st = %d end = %d\n",diag,d_start_pos,d_end_pos);
            }

            syn_bsfs<<<1,BSF_POOL>>>(d_bsf_global);

            cudaMemcpy(&bsf, d_bsf_global, sizeof(FLOAT), cudaMemcpyDeviceToHost);

            FLOAT t4 = clock();
            t_total_dtw += (t4 - t3) ;
            t_total_diag += (t2 - t1) ;

        }
    }

    cudaMemcpy(&bsf, d_bsf_global, sizeof(FLOAT), cudaMemcpyDeviceToHost);

    printf("dtw cnt = %llu\n",dtw_cnt);

    FLOAT t_last = clock();
    printf("dtw time : %fs,diag time : %fs \n",t_total_dtw/CLOCKS_PER_SEC,t_total_diag/CLOCKS_PER_SEC);
    printf("total time is %fs   input of bsf:%5.3f  ,final distance :%f  \n", FLOAT(t_last - t_first)/CLOCKS_PER_SEC, best_so_far, bsf);

    FILE *op = fopen(file,"a");
    fprintf(op, "%lf",  FLOAT(t_last - t_first)/CLOCKS_PER_SEC);
    fprintf(op, "\n");
    fclose(op);

    for (int i = 0; i < subcount; ++i) {
        cudaFree(tmp_my_subs[i]);
        cudaFree(tmp_subs_U[i]);
        cudaFree(tmp_subs_L[i]);
    }
    cudaFree(d_my_subs);
    cudaFree(d_subs_U);
    cudaFree(d_subs_L);
    cudaFree(d_bsf_global);

    free(tmp_my_subs);
    free(tmp_subs_U);
    free(tmp_subs_L);

    cudaFree(d_a);
    cudaFree(d_TS2);
    cudaFree(d_mu);
    cudaFree(d_sig);
    cudaFree(d_sumU_sumL);
    cudaFree(d_invsig);
    cudaFree(d_norm_U_plus_norm_L_trans);

    cudaFree(d_dr_bwdU_plus_dr_bwdL);
    cudaFree(d_dc_bwd);
    cudaFree(d_dr_fwdU_plus_dr_fwdL);
    cudaFree(d_dc_fwd);
    cudaFree(d_UTS);
    cudaFree(d_LTS);
    cudaFree(d_UTS_global);
    cudaFree(d_LTS_global);
    cudaFree(d_MASK_global);
    cudaFree(d_pos_UU);
    cudaFree(d_pos_LL);
    cudaFree(d_DUL2_global);
    cudaFree(d_DUL_global);
    cudaFree(d_DUL2_fast);
    cudaFree(d_DUL_fast);

    cudaFree(d_lb_vector);
    cudaFree(d_lb_vector_new);
    cudaFree(d_indices);
    cudaFree(d_diag);
    cudaFree(d_prefix_sum);
    cudaFree(d_num_inclusive);

    cout << endl;

    free(h_num_of_dtw);
    free(lb_vector);

    for (int j = 0; j < subcount; j++) {
        free(special_shared_vector[j]);
    }
    free(special_shared_vector);

    free(my_subs);
    free(subs_U);
    free(subs_L);

    delete[] h_lb_vector;
    delete[] h_indices;
}

