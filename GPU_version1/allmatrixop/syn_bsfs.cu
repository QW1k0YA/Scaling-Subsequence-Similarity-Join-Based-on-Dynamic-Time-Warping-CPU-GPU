
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