
#ifndef GPU_DTW_COLLECT_INDICES_KERNAL_AFTER_KEOGH_H
#define GPU_DTW_COLLECT_INDICES_KERNAL_AFTER_KEOGH_H
__global__ void
collect_indices_kernel_after_keogh(int *d_indices, bool *lb_vector, int subcount,
                                   int start_pos, int end_pos, int diag,
                                   int *d_counter);
__global__ void
calculate_num_for_each_block(bool *lb_vector, int subcount, int start_pos, int end_pos, int diag, int *global_offset,
                             int *num_inclusive);
__global__ void
collect_indices(int *d_diag, int *d_indices, bool *lb_vector, int subcount, int start_pos, int end_pos, int diag,
                const int *global_offset, const int *num_inclusive);
#endif 
