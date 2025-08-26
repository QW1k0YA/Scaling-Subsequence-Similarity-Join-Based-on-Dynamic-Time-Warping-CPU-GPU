
#include <cuda_runtime.h>
#include "GPU_parameters.h"
#include "matrix.cuh"
__global__ void process_lbkeogh_kernel
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount,
         int diag, int *d_counter, int ran_idx, FLOAT **my_L, FLOAT **my_U, bool *lb_vector,
         FLOAT *cb1, FLOAT *cb2)
{
    
    int *indices_local = &d_indices[subcount*(diag_offset - diag)];
    bool *lb_local = &lb_vector[(diag_offset - diag)*subcount];
    
    FLOAT threshold = bsf[ran_idx];
    FLOAT threshold_2 = threshold * threshold;
    int bid = blockIdx.x;
    int tid = threadIdx.x%32;
    int step = (subseqLen/32);
    
    int d_counter_temp = d_counter[diag_offset - diag];
    for (int i = bid; i < d_counter_temp; i += gridDim.x ) {
        int original_index = indices_local[i];
        FLOAT*  sub1 = my_subs[original_index];
        FLOAT*  sub2 = my_subs[original_index + diag_offset - 1];
        FLOAT* U1 = my_U[original_index];
        FLOAT* L1 = my_L[original_index];
        FLOAT* U2 = my_U[original_index + diag_offset - 1];
        FLOAT* L2 = my_L[original_index + diag_offset - 1];

        FLOAT partial_dist1;
        partial_dist1 = lb_keogh_warp(sub1+ tid *step, U2 + tid*step, L2+ tid*step, step);
        FLOAT partial_dist2;
        partial_dist2 = lb_keogh_warp(sub2+ tid*step, U1+ tid*step, L1+ tid*step, step);

        #pragma unroll
        for (int mask = 1; mask <= 16; mask <<= 1) {
            
            FLOAT tmp1 = __shfl_down_sync(0xFFFFFFFF, partial_dist1, mask);
            FLOAT tmp2 = __shfl_down_sync(0xFFFFFFFF, partial_dist2, mask);
            if (tid + mask < 32) {
                partial_dist1 += tmp1;
                partial_dist2 += tmp2;
            }
        }
        cb1[tid] = partial_dist1;
        cb2[tid] = partial_dist2;

        if (tid == 0 && (partial_dist1 > threshold_2 || partial_dist2 > threshold_2)) {
              
            lb_local[i] = true;

        }
    }
}
