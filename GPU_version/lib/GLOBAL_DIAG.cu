

#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "GPU_parameters.h"
#include "underdtw.cuh"
#include "dtw_511.cuh"
#include <curand_kernel.h>
#include "ECG_511_EARLY_EXIT_DTW.cuh"
using namespace std;
#define DOUBLE_BIAS 6
#define BIAS 3
#define MIN(x,y) ((x)<(y)?(x):(y))

/**
 * Atomic operation to update the minimum float value at a memory address.
 * This function atomically updates the value at the given address to be the
 * minimum of the current value and the provided value, ensuring thread safety
 * in concurrent GPU operations.
 *
 * @param address: Pointer to the memory location to update
 * @param val: Value to compare against the current value
 */
__device__ void atomicMinFloat(FLOAT*  address, FLOAT val) {
        unsigned int* addr_as_ui = (unsigned int*)address;
        unsigned int old = *addr_as_ui, assumed;
        do {
            assumed = old;
            FLOAT current_val = __int_as_float(assumed);
            if (val >= current_val) break;
            old = atomicCAS(addr_as_ui, assumed, __float_as_int(val));
        } while (assumed != old);
    }

/**
 * CUDA kernel for global diagonal computation with multiple pruning layers.
 * This kernel performs the first stage of DTW distance computation by applying
 * multiple levels of lower bound pruning (fast, global) to eliminate unnecessary
 * full DTW calculations on the GPU.
 *
 * @param minlag: Minimum lag between subsequences
 * @param subcount: Total number of subsequences
 * @param subseqLen: Length of each subsequence
 * @param len: Length of the original time series
 * @param warpmax: Maximum warping window size
 * @param a: Original time series data
 * @param mu: Mean values for normalization
 * @param sumU_sumL: Sum of upper and lower envelopes
 * @param invsig: Inverse of standard deviation for normalization
 * @param norm_U_plus_norm_L_trans: Normalized sum of upper and lower envelopes
 * @param dr_bwdU_plus_dr_bwdL: Backward differences of upper and lower envelopes
 * @param dc_bwd: Backward differences of time series
 * @param dr_fwdU_plus_dr_fwdL: Forward differences of upper and lower envelopes
 * @param dc_fwd: Forward differences of time series
 * @param UTS: Upper envelope of the time series
 * @param LTS: Lower envelope of the time series
 * @param UTS_global: Global upper envelope
 * @param LTS_global: Global lower envelope
 * @param MASK_global: Global mask for pruning
 * @param TS2: Squared time series values
 * @param d_lb_vector: Boolean vector for lower bound pruning flags
 * @param d_lb_vector_new: Updated lower bound pruning flags
 * @param sumMASK_global: Sum of global mask values
 * @param sumU_sumL_global: Global sum of upper and lower envelopes
 * @param dr_bwdU_plus_dr_bwdL_global: Global backward differences of envelopes
 * @param dr_fwdU_plus_dr_fwdL_global: Global forward differences of envelopes
 * @param dc_bwd_global: Global backward differences of time series
 * @param dc_fwd_global: Global forward differences of time series
 * @param dr_bwdMASK_global: Global backward differences of mask
 * @param dr_fwdMASK_global: Global forward differences of mask
 * @param dc_bwdTS2_global: Global backward differences of squared time series
 * @param dc_fwdTS2_global: Global forward differences of squared time series
 * @param DUL2_global: Global squared distance between upper and lower envelopes
 * @param DUL_global: Global distance between upper and lower envelopes
 * @param norm_U_plus_norm_L_global: Global normalized sum of envelopes
 * @param my_subs: Array of normalized subsequences
 * @param special_shared_vector: Special shared vectors for optimization
 * @param d_bsf_global: Best-so-far distance value
 * @param DUL_fast: Fast distance between upper and lower envelopes
 * @param DUL2_fast: Squared fast distance between upper and lower envelopes
 * @param diag: Current diagonal being processed
 * @param start_pos: Starting position for processing
 * @param end_pos: Ending position for processing
 */
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
            const FLOAT *DUL2_fast, int diag, int start_pos, int end_pos)
{

    int tid = threadIdx.x + blockDim.x*blockIdx.x;

    diag = diag + tid;
    end_pos = MIN(end_pos, subcount - diag + 1);
    
    if (start_pos > subcount - diag + 1) {
        return;
    }

    bool *lb_vector = &d_lb_vector[tid*STEP_LENGTH];
    FLOAT *lb_vector_new = &d_lb_vector_new[tid*STEP_LENGTH];

    FLOAT bsf = *d_bsf_global;

    __syncthreads();

    FLOAT cnt_of_purn = 0;
    
    diag_fast_segment(a, subseqLen, diag, UTS, LTS, mu, sumU_sumL, invsig, norm_U_plus_norm_L_trans,
              lb_vector, dr_bwdU_plus_dr_bwdL, dc_bwd, dr_fwdU_plus_dr_fwdL, dc_fwd,
              cnt_of_purn,start_pos,end_pos,bsf,DUL_fast,DUL2_fast);


    diag_mask_global(a, subseqLen, diag, lb_vector, mu, UTS_global, LTS_global,
                     MASK_global, TS2,sumMASK_global, invsig, sumU_sumL_global,
                     dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global, dc_bwd_global,
                     dc_fwd_global, dr_bwdMASK_global, dr_fwdMASK_global,
                     dc_bwdTS2_global, dc_fwdTS2_global, DUL2_global, lb_vector_new, bsf,
                     my_subs,special_shared_vector, cnt_of_purn, DUL_global, norm_U_plus_norm_L_global,
                     start_pos,end_pos);
}

__global__ void GLOBAL_DIAG_3TASKS(
    int minlag, int subcount, int subseqLen, int len, int warpmax,
    const FLOAT* a, const FLOAT* mu, const FLOAT* sumU_sumL, const FLOAT* invsig,
    const FLOAT* norm_U_plus_norm_L_trans, const FLOAT* dr_bwdU_plus_dr_bwdL,
    const FLOAT* dc_bwd, const FLOAT* dr_fwdU_plus_dr_fwdL, const FLOAT* dc_fwd,
    const FLOAT* UTS, const FLOAT* LTS, const FLOAT* UTS_global, const FLOAT* LTS_global,
    const FLOAT* MASK_global, const FLOAT* TS2, const FLOAT* sumMASK_global,
    const FLOAT* sumU_sumL_global, const FLOAT* dr_bwdU_plus_dr_bwdL_global,
    const FLOAT* dr_fwdU_plus_dr_fwdL_global, const FLOAT* dc_bwd_global,
    const FLOAT* dc_fwd_global, const FLOAT* dr_bwdMASK_global, const FLOAT* dr_fwdMASK_global,
    const FLOAT* dc_bwdTS2_global, const FLOAT* dc_fwdTS2_global,
    const FLOAT* DUL2_global, const FLOAT* DUL_global, const FLOAT* norm_U_plus_norm_L_global,
    FLOAT** my_subs, FLOAT** special_shared_vector, const FLOAT* d_bsf_global,
    const FLOAT* DUL_fast, const FLOAT* DUL2_fast,
    
    Task t0, Task t1, Task t2,
    
    bool* lb0, FLOAT* lb_new0,
    bool* lb1, FLOAT* lb_new1,
    bool* lb2, FLOAT* lb_new2
) {
    if (threadIdx.x != 0) return; 

    FLOAT bsf = *d_bsf_global;
    FLOAT cnt_of_purn = 0;
    
    {
        int diag = t0.diag;
        int start_pos = t0.start_pos;
        int end_pos = MIN(t0.end_pos, subcount - diag + 1);
        if (start_pos <= subcount - diag + 1) {
            diag_fast_segment(a, subseqLen, diag, UTS, LTS, mu, sumU_sumL, invsig,
                norm_U_plus_norm_L_trans, lb0, dr_bwdU_plus_dr_bwdL, dc_bwd,
                dr_fwdU_plus_dr_fwdL, dc_fwd, cnt_of_purn, start_pos, end_pos, bsf, DUL_fast, DUL2_fast);

            diag_mask_global(a, subseqLen, diag, lb0, mu, UTS_global, LTS_global,
                MASK_global, TS2, sumMASK_global, invsig, sumU_sumL_global,
                dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global,
                dc_bwd_global, dc_fwd_global, dr_bwdMASK_global, dr_fwdMASK_global,
                dc_bwdTS2_global, dc_fwdTS2_global, DUL2_global, lb_new0, bsf,
                my_subs, special_shared_vector, cnt_of_purn, DUL_global,
                norm_U_plus_norm_L_global, start_pos, end_pos);
        }
    }

    
    {
        int diag = t1.diag;
        int start_pos = t1.start_pos;
        int end_pos = MIN(t1.end_pos, subcount - diag + 1);
        if (start_pos <= subcount - diag + 1) {
            diag_fast_segment(a, subseqLen, diag, UTS, LTS, mu, sumU_sumL, invsig,
                norm_U_plus_norm_L_trans, lb1, dr_bwdU_plus_dr_bwdL, dc_bwd,
                dr_fwdU_plus_dr_fwdL, dc_fwd, cnt_of_purn, start_pos, end_pos, bsf, DUL_fast, DUL2_fast);

            diag_mask_global(a, subseqLen, diag, lb1, mu, UTS_global, LTS_global,
                MASK_global, TS2, sumMASK_global, invsig, sumU_sumL_global,
                dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global,
                dc_bwd_global, dc_fwd_global, dr_bwdMASK_global, dr_fwdMASK_global,
                dc_bwdTS2_global, dc_fwdTS2_global, DUL2_global, lb_new1, bsf,
                my_subs, special_shared_vector, cnt_of_purn, DUL_global,
                norm_U_plus_norm_L_global, start_pos, end_pos);
        }
    }

    
    {
        int diag = t2.diag;
        int start_pos = t2.start_pos;
        int end_pos = MIN(t2.end_pos, subcount - diag + 1);
        if (start_pos <= subcount - diag + 1) {
            diag_fast_segment(a, subseqLen, diag, UTS, LTS, mu, sumU_sumL, invsig,
                norm_U_plus_norm_L_trans, lb2, dr_bwdU_plus_dr_bwdL, dc_bwd,
                dr_fwdU_plus_dr_fwdL, dc_fwd, cnt_of_purn, start_pos, end_pos, bsf, DUL_fast, DUL2_fast);

            diag_mask_global(a, subseqLen, diag, lb2, mu, UTS_global, LTS_global,
                MASK_global, TS2, sumMASK_global, invsig, sumU_sumL_global,
                dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global,
                dc_bwd_global, dc_fwd_global, dr_bwdMASK_global, dr_fwdMASK_global,
                dc_bwdTS2_global, dc_fwdTS2_global, DUL2_global, lb_new2, bsf,
                my_subs, special_shared_vector, cnt_of_purn, DUL_global,
                norm_U_plus_norm_L_global, start_pos, end_pos);
        }
    }
}

__global__ void
GLOBAL_DIAG_merge3(int minlag, int subcount, int subseqLen, int len, int warpmax, const FLOAT *a, const FLOAT *mu,
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
            const FLOAT *DUL2_fast, int diag, int start_pos, int end_pos)
{

    int tid = threadIdx.x + blockDim.x*blockIdx.x;

    diag = diag + tid;
    end_pos = MIN(end_pos, subcount - diag + 1);
    
    if (start_pos > subcount - diag + 1) {
        return;
    }

    bool *lb_vector = &d_lb_vector[tid*STEP_LENGTH];
    FLOAT *lb_vector_new = &d_lb_vector_new[tid*STEP_LENGTH];

    FLOAT bsf = *d_bsf_global;

    __syncthreads();

    FLOAT cnt_of_purn = 0;

    diag_fast_segment(a, subseqLen, diag, UTS, LTS, mu, sumU_sumL, invsig, norm_U_plus_norm_L_trans,
              lb_vector, dr_bwdU_plus_dr_bwdL, dc_bwd, dr_fwdU_plus_dr_fwdL, dc_fwd,
              cnt_of_purn,start_pos,end_pos,bsf,DUL_fast,DUL2_fast);

    diag_mask_global(a, subseqLen, diag, lb_vector, mu, UTS_global, LTS_global,
                     MASK_global, TS2,sumMASK_global, invsig, sumU_sumL_global,
                     dr_bwdU_plus_dr_bwdL_global, dr_fwdU_plus_dr_fwdL_global, dc_bwd_global,
                     dc_fwd_global, dr_bwdMASK_global, dr_fwdMASK_global,
                     dc_bwdTS2_global, dc_fwdTS2_global, DUL2_global, lb_vector_new, bsf,
                     my_subs,special_shared_vector, cnt_of_purn, DUL_global, norm_U_plus_norm_L_global,
                     start_pos,end_pos);

}

