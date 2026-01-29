
#include "device_prefix_sum.cuh"
#include <cub/cub.cuh>
#include "GPU_parameters.h"
void device_prefix_sum(int *d_array, int n) {
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, d_array, d_array, n);
    cudaMalloc(&d_temp_storage, temp_storage_bytes);
    cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, d_array, d_array, n);
    cudaFree(d_temp_storage);
}

/**
 * CUDA kernel for performing block-level inclusive prefix sum operation.
 * This function computes the inclusive prefix sum of integer values within a block
 * using the CUB library's block scan primitive, which is used for indexing and
 * counting operations in the DTW computation pipeline.
 *
 * @param d_out: Input/output array for prefix sum computation
 */
__global__ void device_prefix_sum_block(int *d_out)
{

    typedef cub::BlockScan<int, GRID_SIZE> BlockScan;
    __shared__ typename BlockScan::TempStorage temp_storage;
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int thread_data = d_out[tid];
    BlockScan(temp_storage).InclusiveSum(thread_data, thread_data);
    __syncthreads();
    d_out[tid] = thread_data;

}