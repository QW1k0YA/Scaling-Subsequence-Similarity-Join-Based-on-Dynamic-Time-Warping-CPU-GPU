
#include "matrix.cuh"
#include "../alldef/collect_indices_kernal_after_keogh.h"
#include <cub/cub.cuh>
__global__ void
collect_indices_kernel_after_keogh(int *d_indices, bool *lb_vector, int subcount,
                                   int start_pos, int end_pos, int diag,
                                   int *d_counter) {
    int tid =blockIdx.x * blockDim.x + threadIdx.x;
    bool *lb_local = &lb_vector[tid*subcount];
    int* indices_local = &d_indices[tid*subcount];
    diag = diag + tid;

    if(start_pos > subcount - diag + 1) return;
    end_pos = MIN(end_pos, subcount - diag + 1);

    d_counter[tid] = 0;
    for(int i = start_pos;i < end_pos;i++)
    {
        if ( !lb_local[i]) {
            indices_local[d_counter[tid]++] = i;
        }
        else
        {
            lb_local[i] = false;
        }
    }
}

__global__ void
calculate_num_for_each_block(bool *lb_vector, int subcount, int start_pos, int end_pos, int diag, int *global_offset,
                             int *num_inclusive) {
    typedef cub::BlockScan<int, BLOCK_SIZE> BlockScan;
    __shared__ typename BlockScan::TempStorage temp_storage;

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    bool *lb_local = &lb_vector[tid * STEP_LENGTH];
    diag = diag + tid;

    end_pos = MIN(end_pos, subcount - diag + 1);

    int num_in_the_thread = 0;
    int lb_pos = 0;
    for (int i = start_pos; i < end_pos; i++) {

        if (!lb_local[lb_pos]) {
            num_in_the_thread++;
        }

        lb_pos++;
    }

    int init_value = 0;
    __syncthreads();
    BlockScan(temp_storage).InclusiveSum(num_in_the_thread, num_inclusive[tid],init_value);

    if (threadIdx.x == blockDim.x - 1) {
        global_offset[blockIdx.x] = num_inclusive[tid];
    }

}

__global__ void
collect_indices(int *d_diag, int *d_indices, bool *lb_vector, int subcount, int start_pos, int end_pos, int diag,
                const int *global_offset, const int *num_inclusive) {
    
    int tid =blockIdx.x * blockDim.x + threadIdx.x;
    bool *lb_local = &lb_vector[tid*STEP_LENGTH];
    int bid = blockIdx.x;

    int global_off = (bid > 0) ? global_offset[bid - 1] : 0;
    int block_off = (threadIdx.x >0) ? num_inclusive[tid - 1] : 0;
    diag = diag + tid;
    if(start_pos > subcount - diag + 1) return;
    end_pos = MIN(end_pos, subcount - diag + 1);

    int* indices_local = &d_indices[block_off + global_off];
    int* diag_local = &d_diag[block_off + global_off];
    unsigned int pos_for_indices = 0;
    unsigned int pos_for_diag_local = 0;

    int lb_pos = 0;
    for(int i = start_pos;i < end_pos;i++)
    {

        if ( !lb_local[lb_pos]) {
            indices_local[pos_for_indices++] = i;
            diag_local[pos_for_diag_local++] = diag;
        }
        else
        {
            lb_local[lb_pos] = false;
        }
        lb_pos++;
    }
}

