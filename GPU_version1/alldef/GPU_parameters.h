
#ifndef GPU_DTW_LB_GPU_PARAMETERS_H
#define GPU_DTW_LB_GPU_PARAMETERS_H
#include "typedefdouble.cuh"
#include <driver_types.h>
#include "vector"

#define BLOCK_SIZE 64
#define GRID_SIZE (6*30)

#define STEP_LENGTH 1025
#define STREAM_POOL_SIZE_DTW 8

#define STREAM_POOL_SIZE_KEOGH 16
#define BSF_POOL 32

#define CB_NUM 4
#define CB_LEN 256

#ifdef __CUDACC__
#define CUERR {                                                            \
        cudaError_t err;                                                       \
        if ((err = cudaGetLastError()) != cudaSuccess) {                       \
            std::cout << "CUDA error: " << cudaGetErrorString(err) << " : "    \
                      << __FILE__ << ", line " << __LINE__ << std::endl;       \
            exit(1);                                                           \
        }                                                                      \
    }
#endif

using namespace std;

__device__ FLOAT  elementWiseMultiply_p_plus_sum(const FLOAT  *vector1, const FLOAT  *vector2, const FLOAT  *vector3, long long len);
__device__ void diag_fast(const FLOAT  *ts, int subseqlen, int diag_ID, const FLOAT *UTS,const FLOAT *LTS,
                          const FLOAT *mu,  const FLOAT *sumU_sumL, const FLOAT *invsig,
                          const FLOAT *norm_U_plus_norm_L_trans,  const FLOAT *del,bool *lb_vector,
                          const FLOAT *dr_bwdU_plus_dr_bwdL, const FLOAT *dc_bwd,
                          const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_fwd, FLOAT &cnt,int EPOCH);
__device__ void diag_fast_segment(const FLOAT *ts, int subseqlen, int diag_ID, const FLOAT *UTS, const FLOAT *LTS,
                                  const FLOAT *mu, const FLOAT *sumU_sumL, const FLOAT *invsig,
                                  const FLOAT *norm_U_plus_norm_L_trans, bool *lb_vector,
                                  const FLOAT *dr_bwdU_plus_dr_bwdL, const FLOAT *dc_bwd,
                                  const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_fwd, FLOAT &cnt,
                                  int start_pos, int end_pos, FLOAT bsf,
                                  const FLOAT *DUL, const FLOAT *DUL2);
__device__ bool LB_KK(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen, FLOAT threshold_2,
                      FLOAT *cb, FLOAT special_shared_vector[], FLOAT &lbk);

__device__ void
diag_mask_global(const FLOAT *ts, int subseqlen, int diagID, bool *lb_vector, const FLOAT *mu, const FLOAT *UTS,
                 const FLOAT *LTS, const FLOAT *MASK, const FLOAT *TS2, const FLOAT *sumMASK,
                 const FLOAT *invsig, const FLOAT *sumU_sumL, const FLOAT *dr_bwdU_plus_dr_bwdL,
                 const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_bwd, const FLOAT *dc_fwd,
                 const FLOAT *dr_bwdMASK, const FLOAT *dr_fwdMASK, const FLOAT *dc_bwdTS2,
                 const FLOAT *dc_fwdTS2, const FLOAT *DUL2, FLOAT *lb_vector_new, FLOAT  bsf, FLOAT **subs,
                 FLOAT **special_shared_vector, FLOAT &cnt, const FLOAT *DUL,
                 const FLOAT *norm_U_plus_norm_L_global, int start_pos, int end_pos);
__device__ void atomicMinFloat(FLOAT*  address, FLOAT val);
__global__ void traversal_dtw_array(const FLOAT *dtw_array, FLOAT *bsf, int size);
__device__ FLOAT MON_dtw(const FLOAT *lines, const FLOAT *cols, int l, int w, FLOAT bsf, FLOAT *buffers, int tid);

__global__ void
GLOBAL_DIAG(int minlag, int subcount, int subseqLen, int len, int warpmax, const FLOAT *a, const FLOAT *mu,
            const FLOAT *sumU_sumL, const FLOAT *invsig, const FLOAT *norm_U_plus_norm_L_trans,
            const FLOAT *dr_bwdU_plus_dr_bwdL, const FLOAT *dc_bwd, const FLOAT *dr_fwdU_plus_dr_fwdL,
            const FLOAT *dc_fwd, const FLOAT *UTS, const FLOAT *LTS, const FLOAT *UTS_global,
            const FLOAT *LTS_global, const FLOAT *MASK_global, const FLOAT *TS2, bool *d_lb_vector,
            FLOAT *d_lb_vector_new, const FLOAT *sumMASK_global, const FLOAT *sumU_sumL_global,
            const FLOAT *dr_bwdU_plus_dr_bwdL_global, const FLOAT *dr_fwdU_plus_dr_fwdL_global,
            const FLOAT *dc_bwd_global, const FLOAT *dc_fwd_global, const FLOAT *dr_bwdMASK_global,
            const FLOAT *dr_fwdMASK_global, const FLOAT *dc_bwdTS2_global, const FLOAT *dc_fwdTS2_global,
            const FLOAT *DUL2_global, const FLOAT *DUL_global, const FLOAT *norm_U_plus_norm_L_global,
            FLOAT **my_subs, FLOAT **special_shared_vector, const FLOAT *d_bsf_global, const FLOAT *DUL_fast,
            const FLOAT *DUL2_fast, int diag, int start_pos, int end_pos);
__device__ FLOAT dtw(const FLOAT *A, const FLOAT *B, int m, int r, FLOAT threshold_2,   FLOAT *cost, FLOAT *cost_prev);
__global__ void test(FLOAT *dtw_array, int subcount, FLOAT *global_bsf);
__global__ void test_array(FLOAT *dtw_array,FLOAT *dtw_array_debug,int subcount);
__device__ void early_exit_DTW511_blockfunc( 
        FLOAT  * Subject,
        FLOAT  * cQuery,
        FLOAT  &Dist,
        int num_features,
        FLOAT  threshold);
__device__ void early_exit_DTW63_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold, int w);
__device__ void early_exit_DTW127_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold, int w);
__device__ void early_exit_DTW255_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold, int w);
__device__ void early_exit_DTW511_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold, int w);
__device__ void
early_exit_DTW1023_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold, int w);
__device__ void
early_exit_DTW2047_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold, int w);
__device__ void
early_exit_DTW_256_to_511_with_band(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int num_features, FLOAT threshold,
                                    int w, FLOAT *special_shared_memory);
__device__ void early_exit_DTW511_without_band( 
        FLOAT  * Subject,
        FLOAT  * cQuery,
        FLOAT  &Dist,
        int num_features,
        FLOAT  threshold);
__global__ void process_CDTW
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int w, const int *indices, const int *diag_of_indices);
__global__ void process_dtw_kernel
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount, int diag,
         int *d_counter, int ran_idx, int w, FLOAT *cb1, FLOAT *cb2);
__global__ void process_lbkeogh_kernel
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount,
         int diag, int *d_counter, int ran_idx, FLOAT **my_L, FLOAT **my_U, bool *lb_vector,
         FLOAT *cb1, FLOAT *cb2);
__global__ void process_dtw_kernel_debug
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount, int diag,
         int *d_counter, int ran_idx, int w);
__global__ void process_dtw_kernel_256_to_511
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount,
         int diag, int *d_counter, int ran_idx, int w, FLOAT *special_shared_vector);
__global__ void process_MON_dtw_kernel
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount,
         int diag, int *d_counter, int w, int ran_idx);
__global__ void
collect_indices_kernel(int *d_indices, bool *lb_vector, int subcount, int start_pos, int end_pos, int diag,
                       int *d_counter);
__global__ void syn_bsfs(FLOAT *bsf);

__device__ void
DTW_stairs(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold_2, int w, FLOAT *q, FLOAT *t, FLOAT cb1,
           FLOAT cb2);
__device__ void
DTW_stairs_without_shared(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold_2, int w, const FLOAT cb[]);
__device__ void
DTW_stairs_for_block_without_shared(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold, int w,
                                    int bl_size, FLOAT cb[], FLOAT threshold_2);
__device__ void CDTw_for_block(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int M, FLOAT threshold, int w,
                               FLOAT *sh_Query, FLOAT *sh_Subject,FLOAT *sh_Penalty);
__device__ void
DTW_stairs_test(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold_2, int w, FLOAT *q, FLOAT *t,
                FLOAT cb1, FLOAT cb2, int *cb_prune_time, int low_index, int high_index);
__device__ void
DTW_stairs_over_warp(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold, int w, FLOAT *q, FLOAT *t,FLOAT* shared_up,FLOAT* shared_down);
__device__ void
DTW_stairs_for_block(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold, int w, FLOAT *q, FLOAT *t,
                     int bl_size, FLOAT threshold_2);
__device__ FLOAT lb_keogh_warp(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen);
__device__ FLOAT lb_keogh_warp_with_nomalise(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen,const FLOAT* mu,const FLOAT* invsig);
__device__ FLOAT lb_keogh_stride(const FLOAT *q, const FLOAT *U, const FLOAT *L, long long seqlen);
__device__ FLOAT lb_keogh_stride_with_nomalise(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen,
                                               const FLOAT mu, const FLOAT invsig);
#endif 
