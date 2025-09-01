

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