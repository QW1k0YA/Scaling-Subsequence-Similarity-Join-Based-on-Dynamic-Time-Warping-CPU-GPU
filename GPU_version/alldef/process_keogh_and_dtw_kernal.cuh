
#include "GPU_parameters.h"
#ifndef GPU_DTW_PROCESS_KEOGH_AND_DTW_KERNAL_CUH
#define GPU_DTW_PROCESS_KEOGH_AND_DTW_KERNAL_CUH
__global__ void process_keogh_and_dtw_kernel
        (FLOAT **my_subs,  FLOAT **my_L, FLOAT **my_U,int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount, int diag,
         int *d_counter, int ran_idx, int w);
__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int subcount, int w, const int *indices, const int *diag_of_indices);
__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram_keogh_prune
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int subcount, int w, const int *indices, const int *diag_of_indices,
         int *num_prune);
__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram_without_shared_memory
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int subcount, int w, const int *indices, const int *diag_of_indices,
         int num_of_dtw);
__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram_without_shared_memory_and_nomalized
        (FLOAT **my_subs, const FLOAT *UTS,
         const FLOAT *LTS, const FLOAT *mu,
         const FLOAT *invsig, int subseqLen,
         FLOAT *bsf, int subcount, int w,
         const int *indices,
         const int *diag_of_indices,
         int bl_size);
#endif 
