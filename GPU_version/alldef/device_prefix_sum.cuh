
#ifndef GPU_DTW_DEVICE_PREFIX_SUM_CUH
#define GPU_DTW_DEVICE_PREFIX_SUM_CUH
void device_prefix_sum(int *d_array, int n);
__global__ void device_prefix_sum_block(int *d_out);
#endif 
