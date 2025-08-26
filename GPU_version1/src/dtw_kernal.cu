
#include <cuda_runtime.h>
#include "GPU_parameters.h"
#include "matrix.cuh"
#include "warp_num_per_block.h"

#define THREAD_NUM_PER_WARP 32

__global__ void process_dtw_kernel
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount, int diag,
         int *d_counter, int ran_idx, int w, FLOAT *cb1, FLOAT *cb2)
{
    
    int *indices_local = &d_indices[subcount*(diag_offset - diag)];
    extern __shared__ FLOAT shared_mem[];
    
    FLOAT* q = &shared_mem[0];                  
    FLOAT* t = &shared_mem[subseqLen*WARP_NUMS];

    FLOAT threshold = bsf[ran_idx];
    FLOAT threshold_2 = threshold * threshold;
    int bid = blockIdx.x;

    for (int i = bid; i < d_counter[diag_offset - diag]; i += gridDim.x*WARP_NUMS ) {
        int i_block = threadIdx.x/32;
        int original_index = indices_local[i*WARP_NUMS+i_block];

        {
            FLOAT*  sub1 = my_subs[original_index];
            FLOAT*  sub2 = my_subs[original_index + diag_offset - 1];

            if(original_index + diag_offset - 1 >= subcount){
                printf("dtw_kernal 67 %d %d\n",original_index ,diag_offset - 1);
            }
            FLOAT result = 0.0f;

            int bl_size;
            bl_size= ceil(w/31.0);

            if(w < 32)
            {

            }
            else{
                DTW_stairs_for_block(sub1, sub2, result, subseqLen, threshold_2, w, q, t,bl_size);
            }

            __syncthreads();

            if(threadIdx.x%32 == 16)
            {

                if(result < threshold)
                {
                    atomicMinFloat(&bsf[ran_idx], result);
                    __threadfence();
                }
            }

        }

    }
}

__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int subcount, int w, const int *indices, const int *diag_of_indices,
         int num_of_dtw)
{
    extern __shared__ FLOAT shared_mem[];
    
    FLOAT* q = &shared_mem[0];                  
    FLOAT* t = &shared_mem[subseqLen*WARP_NUMS];

    FLOAT cb1;
    FLOAT cb2;
    
    int bid = blockIdx.x;
    int tid = threadIdx.x%32;
    int step = (subseqLen/32);
    FLOAT threshold = bsf[bid%BSF_POOL]; 
    FLOAT threshold_2 = threshold * threshold;

    int i = bid;

    {
        int original_index = indices[i];
        int diag_offset = diag_of_indices[i];

        if(diag_offset != 0) 
        {

            FLOAT*  sub1 = my_subs[original_index];
            FLOAT*  sub2 = my_subs[original_index + diag_offset - 1];

            FLOAT* U1 = my_U[original_index];
            FLOAT* L1 = my_L[original_index];
            FLOAT* U2 = my_U[original_index + diag_offset - 1];
            FLOAT* L2 = my_L[original_index + diag_offset - 1];

            FLOAT partial_dist1 = 0;
        partial_dist1 = lb_keogh_warp(sub1+ tid *step, U2 + tid*step, L2+ tid*step, step);
            FLOAT partial_dist2 = 0;
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
            cb1 = partial_dist1;
            cb2 = partial_dist2;

            bool flag = 1;
            if (tid == 0 && (partial_dist1 > threshold_2 || partial_dist2 > threshold_2)) {
                
                flag = false; 

            }
            flag = __shfl_sync(0x1F, flag, 0);
            __syncwarp();

            if(flag)
            {

                FLOAT result = 0.0f;

                int bl_size;
                bl_size= ceil(w/31.0);

                if(w < 32)
                {
                    DTW_stairs(sub1, sub2, result, subseqLen, threshold_2, w, q, t, cb1, cb2);
                }
                else{
                    DTW_stairs_for_block(sub1, sub2, result, subseqLen, threshold_2, w, q, t,bl_size);
                }

                __syncwarp();
                if(threadIdx.x%32 == 16)
                {
                    if(result < threshold)
                    {

                        atomicMinFloat(&bsf[bid%BSF_POOL], result);

                        __threadfence();
                    }
                }
            }
        }

    }

}

__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram_without_shared_memory
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int subcount, int w, const int *indices, const int *diag_of_indices,
         int num_of_dtw)
{
    extern __shared__ FLOAT shared_mem[];
    
    FLOAT* q = &shared_mem[0];                  
    FLOAT* t = &shared_mem[subseqLen*WARP_NUMS];

    FLOAT cb1;
    FLOAT cb2;
    
    int bid = blockIdx.x;
    int tid = threadIdx.x%32;
    int step = (subseqLen/32);
    FLOAT threshold = bsf[bid%BSF_POOL]; 
    FLOAT threshold_2 = threshold * threshold;

    int i = bid;

    {
        int original_index = indices[i];
        int diag_offset = diag_of_indices[i];

        if(diag_offset != 0) 
        {

            FLOAT*  sub1 = my_subs[original_index];
            FLOAT*  sub2 = my_subs[original_index + diag_offset - 1];

            FLOAT* U1 = my_U[original_index];
            FLOAT* L1 = my_L[original_index];
            FLOAT* U2 = my_U[original_index + diag_offset - 1];
            FLOAT* L2 = my_L[original_index + diag_offset - 1];

            FLOAT partial_dist1 = 0;
            partial_dist1 = lb_keogh_warp(sub1+ tid *step, U2 + tid*step, L2+ tid*step, step);
            FLOAT partial_dist2 = 0;
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
            cb1 = partial_dist1;
            cb2 = partial_dist2;

            bool flag = 1;
            if (tid == 0 && (partial_dist1 > threshold_2 || partial_dist2 > threshold_2)) {
                
                flag = false; 

            }
            flag = __shfl_sync(0x1F, flag, 0);
            __syncwarp();

            if(flag)
            {

                FLOAT result = 0.0f;

                int bl_size;
                bl_size= ceil(w/31.0);

                if(w < 32)
                {
                    DTW_stairs_without_shared(sub1, sub2, result, subseqLen, threshold_2, w, q, t);
                }
                else{
                    DTW_stairs_for_block_without_shared(sub1, sub2, result, subseqLen, threshold_2, w, q, t,bl_size);
                }

                __syncwarp();
                if(threadIdx.x%32 == 16)
                {

                    if(result < threshold)
                    {

                        atomicMinFloat(&bsf[bid%BSF_POOL], result);

                        __threadfence();
                    }
                }
            }
        }

    }

}

__global__ void process_keogh_and_dtw_kernel_for_a_Parallelogram_keogh_prune
        (FLOAT **my_subs, FLOAT **my_L, FLOAT **my_U, int subseqLen, FLOAT *bsf,
         int subcount, int w, const int *indices, const int *diag_of_indices,
         int *num_prune)
{
    extern __shared__ FLOAT shared_mem[];
    
    FLOAT* q = &shared_mem[0];                  
    FLOAT* t = &shared_mem[subseqLen*WARP_NUMS];

    FLOAT cb1;
    FLOAT cb2;
    
    int bid = blockIdx.x;
    int tid = threadIdx.x%32;
    int step = (subseqLen/32);
    FLOAT threshold = bsf[bid%BSF_POOL]; 
    FLOAT threshold_2 = threshold * threshold;

    int i = bid;

    {
        int original_index = indices[i];
        int diag_offset = diag_of_indices[i];

        if(diag_offset != 0) 
        {

            FLOAT*  sub1 = my_subs[original_index];
            FLOAT*  sub2 = my_subs[original_index + diag_offset - 1];

            FLOAT* U1 = my_U[original_index];
            FLOAT* L1 = my_L[original_index];
            FLOAT* U2 = my_U[original_index + diag_offset - 1];
            FLOAT* L2 = my_L[original_index + diag_offset - 1];

            FLOAT partial_dist1 = 0;
            partial_dist1 = lb_keogh_warp(sub1+ tid *step, U2 + tid*step, L2+ tid*step, step);
            FLOAT partial_dist2 = 0;
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

            cb1 = partial_dist1;
            cb2 = partial_dist2;

            bool flag = 1;
            if (tid == 0 && (partial_dist1 > threshold_2 || partial_dist2 > threshold_2)) {

                flag = false; 

            }

            flag = __shfl_sync(0x1F, flag, 0);

            __syncwarp();
            if(flag)
            {

                FLOAT result = 0.0f;

                int bl_size;
                bl_size= ceil(w/31.0);

                if(w < 32)
                {

                    DTW_stairs(sub1, sub2, result, subseqLen, threshold_2, w, q, t, cb1, cb2);
                }
                else{
                    DTW_stairs_for_block(sub1, sub2, result, subseqLen, threshold_2, w, q, t,bl_size);
                }

                __syncwarp();
                if(threadIdx.x%32 == 16)
                {
                    if(result < threshold)
                    {
                        atomicMinFloat(&bsf[bid%BSF_POOL], result);

                        __threadfence();
                    }
                }
            }
        }

    }

}

__global__ void process_MON_dtw_kernel
        (FLOAT **my_subs, int *d_indices, int subseqLen, int diag_offset, FLOAT *bsf, int subcount,
         int diag, int *d_counter, int w, int ran_idx)
{
    
    extern __shared__ FLOAT mon_buffers[];
    
    int *indices_local = &d_indices[subcount*(diag_offset - diag)];

    FLOAT threshold = bsf[ran_idx];
    FLOAT threshold_2 = threshold * threshold;
    int tid = blockIdx.x*blockDim.x + threadIdx.x;
    int stride = blockDim.x*gridDim.x;

    FLOAT *local_buffer = &mon_buffers[threadIdx.x*(2*subseqLen + 2)];

    for (int i = tid; i < d_counter[diag_offset - diag]; i += stride ) {

        int original_index = indices_local[i];

        indices_local[i] = 0;
        
        FLOAT*  sub1 = my_subs[original_index];
        FLOAT*  sub2 = my_subs[original_index + diag_offset - 1];

        if(original_index + diag_offset - 1 >= subcount){
            printf("dtw_kernal 67 %d %d\n",original_index ,diag_offset - 1);
        }
        
        FLOAT result = 0.0f;
        result = sqrt(MON_dtw(sub1, sub2, subseqLen, w, threshold_2, local_buffer,tid));
        if(result < 5)
        {
            printf("oi !!!%d %d %f\n",original_index,original_index + diag_offset - 1,result);
        }

        if(threadIdx.x == blockDim.x - 1)
        {
            if(diag_offset == 6544 && original_index == 2396)
            {
               printf("dtw_value = %f \n" ,result);
            }

            if(result < threshold)
            {
                atomicMinFloat(&bsf[ran_idx], result);
                __threadfence();
            }
        }

    }
}