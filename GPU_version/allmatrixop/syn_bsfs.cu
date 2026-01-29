
#include <cuda_runtime.h>
#include "GPU_parameters.h"
#include "matrix.cuh"

__global__ void syn_bsfs(FLOAT *bsf)
{
    __shared__ FLOAT bsf_array[BSF_POOL];
    int tid = threadIdx.x;
    bsf_array[tid] = bsf[tid];

    for (int stride = blockDim.x/2; stride > 0; stride >>= 1) {
        if (tid < stride && tid + stride < BSF_POOL) {
            bsf_array[tid] = fminf(bsf_array[tid],bsf_array[tid + stride]);
        }
        __syncthreads();
    }

    bsf[tid] = bsf_array[0];

}

/**
 * CUDA kernel for synchronizing best-so-far (bsf) values across three streams.
 * This function combines the best-so-far distance values from three different
 * computation streams by finding the minimum value across all streams and
 * updating all three arrays with this minimum value to maintain consistency.
 *
 * @param bsf: First array of best-so-far distance values
 * @param bsf1: Second array of best-so-far distance values
 * @param bsf2: Third array of best-so-far distance values
 */
__global__ void syn_bsfs_3(FLOAT *bsf,FLOAT *bsf1,FLOAT *bsf2)
{
    __shared__ FLOAT bsf_array[BSF_POOL];
    int tid = threadIdx.x;
    if (blockIdx.x == 0) bsf_array[tid] = bsf[tid];
    if (blockIdx.x == 1) bsf_array[tid] = bsf1[tid];
    if (blockIdx.x == 2) bsf_array[tid] = bsf2[tid];

    for (int stride = blockDim.x/2; stride > 0; stride >>= 1) {
        if (tid < stride && tid + stride < BSF_POOL) {
            bsf_array[tid] = fminf(bsf_array[tid],bsf_array[tid + stride]);
        }
        __syncthreads();
    }

    if (blockIdx.x == 0) bsf[tid] = bsf_array[0];
    if (blockIdx.x == 1) bsf1[tid] = bsf_array[0];
    if (blockIdx.x == 2) bsf2[tid] = bsf_array[0];

    FLOAT min_bsf = MIN(MIN(bsf[0],bsf1[0]),bsf2[0]);
    bsf[tid] = min_bsf;
    bsf1[tid] = min_bsf;
    bsf2[tid] = min_bsf;

}