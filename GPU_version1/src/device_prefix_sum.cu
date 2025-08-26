
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