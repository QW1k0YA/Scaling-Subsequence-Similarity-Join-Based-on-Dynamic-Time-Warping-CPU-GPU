
#include <cuda_runtime.h>
#include "GPU_parameters.h"
#include "matrix.cuh"
#define REGISTER_NUM 2
#define THREAD_NUM_PER_WARP 32

__device__ void
DTW_stairs_for_block_without_shared(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold, int w,
                                    int bl_size, FLOAT cb[], FLOAT threshold_2) {
    
    bl_size = REGISTER_NUM;

    int num_per_bl = bl_size*bl_size;
    int tid = threadIdx.x;
    
    int num_tid = 32;
    size_t vote;
    
    FLOAT *q =cQuery;
    FLOAT *t = Subject;
    
    __syncthreads();
    int row_bias;
    int col_bias;
    int mid_tid;
    
    {

        mid_tid = 16;
        row_bias = (mid_tid-tid)*bl_size;
        col_bias = (tid - mid_tid)*bl_size;
    }

    int i_temp = row_bias;
    int j_temp = col_bias;

    FLOAT DTW_FIR[REGISTER_NUM*REGISTER_NUM];
    FLOAT DTW_SEC[REGISTER_NUM*REGISTER_NUM];
    for(int i = 0;i < num_per_bl;i++)
    {
        DTW_FIR[i] = INFINITY;
        DTW_SEC[i] = INFINITY;
        
    }
    if (tid == mid_tid) DTW_FIR[num_per_bl-1] = 0;

    FLOAT DTW_UP[REGISTER_NUM];
    FLOAT DTW_DOWN[REGISTER_NUM];
    
    for(int i = 0;i < bl_size;i++)
    {
        DTW_UP[i] = INFINITY;
        DTW_DOWN[i] = INFINITY;
    }

    FLOAT q0,t0,t1;
    int cb_index ;
    FLOAT cb_temp;
    int switch_for_stair = 0;
    bool mask;
    FLOAT d;
    bool flag_pruning = false;

    int target_j_temp = (m/2) - (m/2)%bl_size - 1;

    int w_bias = ceil((w + 1.0)/bl_size);
    for(int step = 0; step < w_bias; step++){
        mask = switch_for_stair%2;
        switch_for_stair++;
        
        if(!mask)
        {
            for(int up_i = 0; up_i < bl_size; up_i++)
            {
                DTW_UP[up_i] = __shfl_up_sync(0xFFFFFFFF, DTW_SEC[num_per_bl - bl_size + up_i], 1, 32);
                
            }

            d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
            DTW_FIR[0] = d + MIN(DTW_UP[0],MIN(DTW_SEC[bl_size - 1],DTW_FIR[num_per_bl-1]));

            i_temp++;

            for(int i = 1;i < bl_size;i++)
            {

                d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                DTW_FIR[i] = d + MIN(DTW_UP[i],MIN(DTW_UP[i-1],DTW_FIR[i-1]));

                i_temp++;
            }
            
            j_temp++;
            i_temp-=(bl_size);
            
            for(int i = 1;i < bl_size;i++)
            
            {

                d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                DTW_FIR[i*bl_size] =  d + MIN(DTW_SEC[(i+1)*bl_size - 1],MIN(DTW_SEC[(i)*bl_size - 1],DTW_FIR[(i-1)*bl_size]));

                i_temp++;
                
                for(int j = 1;j < bl_size;j++)
                
                {

                    d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                    DTW_FIR[i*bl_size + j] =  d + MIN(DTW_FIR[i*bl_size + j - 1],
                                                      MIN(DTW_FIR[(i-1)*bl_size + j - 1],DTW_FIR[(i-1)*bl_size + j]));

                    i_temp++;
                }
                j_temp++;
                i_temp-=(bl_size);
                
            }

        }
        else
        {
            for(int down_i = 0; down_i < bl_size; down_i++)
            {
                DTW_DOWN[down_i] = __shfl_down_sync(0xFFFFFFFF, DTW_FIR[(down_i+1)*bl_size - 1], 1, 32);
                
            }

            d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
            DTW_SEC[0] = d + MIN(DTW_DOWN[0],MIN(DTW_FIR[num_per_bl - bl_size],
                                                                         DTW_SEC[num_per_bl - 1]));

            j_temp++;
            
            for(int i = 1;i < bl_size;i++)
            {

                d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                DTW_SEC[i*bl_size] = d + MIN(DTW_SEC[(i-1)*bl_size],MIN(DTW_DOWN[i-1],DTW_DOWN[i]));

                j_temp++;
            }
            
            i_temp++;
            j_temp-=(bl_size);
            
            for(int i = 1;i < bl_size;i++)
            {

                d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                DTW_SEC[i] =  d +
                        MIN(DTW_SEC[i - 1], MIN(DTW_FIR[num_per_bl - bl_size + i],DTW_FIR[num_per_bl - bl_size + i-1]));

                i_temp++;
            }
            i_temp-=(bl_size - 1);
            j_temp++;
            
            for(int i = 1;i < bl_size;i++)
                
            {
                for(int j = 1;j < bl_size;j++)
                    
                {

                    d = (i_temp >= 0 && j_temp >= 0) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                    DTW_SEC[i*bl_size + j] =  d + MIN(DTW_SEC[i*bl_size + j - 1],
                                                                              MIN(DTW_SEC[(i-1)*bl_size + j - 1],DTW_SEC[(i-1)*bl_size + j]));

                    i_temp++;
                }
                j_temp++;
                i_temp-=(bl_size - 1);
                
            }

            j_temp-=(bl_size);
            i_temp +=(bl_size - 1);

        }

    }

    __syncthreads();

    FLOAT x,y,z;
    
    for(int step = w_bias; step < ceil(2.0*m/bl_size) - w_bias - 1; step++){
        mask = switch_for_stair%2;
        switch_for_stair++;
        
        if(!mask)
        {
            for(int up_i = 0; up_i < bl_size; up_i++)
            {
                DTW_UP[up_i] = __shfl_up_sync(0xFFFFFFFF, DTW_SEC[num_per_bl - bl_size + up_i], 1, 32);
                
            }

            d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);

            DTW_FIR[0] = d + MIN(DTW_UP[0],MIN(DTW_SEC[bl_size - 1],DTW_FIR[num_per_bl-1]));

            i_temp++;
            
            for(int i = 1;i < bl_size;i++)
            {
                d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);
                DTW_FIR[i] = d + MIN(DTW_UP[i],MIN(DTW_UP[i-1],DTW_FIR[i-1]));

                i_temp++;
            }

            j_temp++;
            i_temp-=(bl_size);
            
            for(int i = 1;i < bl_size;i++)
                
            {

               d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);

                DTW_FIR[i*bl_size] =  d +
                        MIN(DTW_SEC[(i+1)*bl_size - 1],MIN(DTW_SEC[(i)*bl_size - 1],DTW_FIR[(i-1)*bl_size]));

                i_temp++;
                
                for(int j = 1;j < bl_size;j++)
                    
                {

                    d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);

                    DTW_FIR[i*bl_size + j] =  d +
                            MIN(DTW_FIR[i*bl_size + j - 1],MIN(DTW_FIR[(i-1)*bl_size + j - 1],DTW_FIR[(i-1)*bl_size + j]));

                    bool ifcb = (tid == 16) && (j_temp == target_j_temp) && (j == bl_size - 1) && (i == bl_size - 1) ;
                    ifcb = __shfl_sync(0x1F, ifcb, 16);
                    if(ifcb)
                    {
                        cb_index =(w + j_temp)/CB_LEN+1;
                        cb_temp =cb[cb_index];
                        bool all_greater = true;
                        for (int index_of_fir = 0; index_of_fir < bl_size; ++index_of_fir) {
                            if (DTW_FIR[index_of_fir] <= threshold_2 - cb_temp) {
                                all_greater = false;
                                break;
                            }
                        }

                        for (int index_of_fir = 1; index_of_fir < bl_size; ++index_of_fir) {
                            int left_index = index_of_fir*bl_size;
                            if (DTW_FIR[left_index] <= threshold_2 - cb_temp) {
                                all_greater = false;
                                break;
                            }
                        }

                        vote = __ballot_sync(0xFFFFFFFF, all_greater);
                        if (vote == 0xFFFFFFFF) {
                            if(tid == 16)
                            {
                                Dist = INFINITY;

                            }
                            flag_pruning = true;
                            return;
                        }
                    }
                    i_temp++;
                }
                j_temp++;
                i_temp-=(bl_size);
                
            }
        }
            
        else
        {
            for(int down_i = 0; down_i < bl_size; down_i++)
            {
                DTW_DOWN[down_i] = __shfl_down_sync(0xFFFFFFFF, DTW_FIR[(down_i+1)*bl_size - 1], 1, 32);
                
            }

            d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY  :DIST(q[i_temp],t[j_temp]);

            DTW_SEC[0] = d + MIN(DTW_DOWN[0],MIN(DTW_FIR[num_per_bl - bl_size],
                                                                         DTW_SEC[num_per_bl - 1]));

            j_temp++;

            for(int i = 1;i < bl_size;i++)
            {
                d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);

                DTW_SEC[i*bl_size] = d + MIN(DTW_SEC[(i-1)*bl_size],MIN(DTW_DOWN[i-1],DTW_DOWN[i]));

                j_temp++;
            }
            
            i_temp++;
            j_temp-=(bl_size);
            
            for(int i = 1;i < bl_size;i++)
            {
                d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);

                DTW_SEC[i] =  d +
                              MIN(DTW_SEC[i - 1], MIN(DTW_FIR[num_per_bl - bl_size + i],DTW_FIR[num_per_bl - bl_size + i-1]));

                i_temp++;
            }
            i_temp-=(bl_size - 1);
            j_temp++;
            
            for(int i = 1;i < bl_size;i++)
                
            {
                for(int j = 1;j < bl_size;j++)
                    
                {
                    d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY : DIST(q[i_temp],t[j_temp]);

                    DTW_SEC[i*bl_size + j] = d + MIN(DTW_SEC[i*bl_size + j - 1],
                                                                              MIN(DTW_SEC[(i-1)*bl_size + j - 1],DTW_SEC[(i-1)*bl_size + j]));

                    i_temp++;
                }
                j_temp++;
                i_temp-=(bl_size - 1);
                
            }
            j_temp-=bl_size;
            i_temp +=(bl_size - 1);
        }

    }

    __syncthreads();

    if(!flag_pruning)
    {
    for(int step = ceil(2.0*m/bl_size) - w_bias - 1; step <  ceil(2.0*m/bl_size); step++){
        mask = switch_for_stair%2;
        switch_for_stair++;
        
        if(!mask)
        {
            for(int up_i = 0; up_i < bl_size; up_i++)
            {
                DTW_UP[up_i] = __shfl_up_sync(0xFFFFFFFF, DTW_SEC[num_per_bl - bl_size + up_i], 1, 32);
                
            }

            d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

            DTW_FIR[0] = d + MIN(DTW_UP[0],MIN(DTW_SEC[bl_size - 1],DTW_FIR[num_per_bl-1]));

            if(tid == mid_tid && i_temp == m - 1 && j_temp == m - 1)
            {
                Dist = sqrt(DTW_FIR[0]);

            }
            i_temp++;
            
            for(int i = 1;i < bl_size;i++)
            {
                d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

                DTW_FIR[i] = d + MIN(DTW_UP[i],MIN(DTW_UP[i-1],DTW_FIR[i-1]));

                if(tid == mid_tid && i_temp == m - 1 && j_temp == m - 1)
                {
                    Dist = sqrt(DTW_FIR[i]);

                }
                i_temp++;
            }
            
            j_temp++;
            i_temp-=(bl_size);
            
            for(int i = 1;i < bl_size;i++)
                
            {
                d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

                DTW_FIR[i*bl_size] =  d + MIN(DTW_SEC[(i+1)*bl_size - 1],MIN(DTW_SEC[(i)*bl_size - 1],DTW_FIR[(i-1)*bl_size]));

                if(tid == mid_tid && i_temp == m - 1  && j_temp == m - 1 )
                {
                    Dist = sqrt(DTW_FIR[i*bl_size]);

                }
                i_temp++;
                
                for(int j = 1;j < bl_size;j++)
                    
                {
                    d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

                    DTW_FIR[i*bl_size + j] =  d + MIN(DTW_FIR[i*bl_size + j - 1],
                                                                              MIN(DTW_FIR[(i-1)*bl_size + j - 1],DTW_FIR[(i-1)*bl_size + j]));

                    if(tid == mid_tid && i_temp == m - 1 && j_temp == m - 1)
                    {
                        Dist = sqrt(DTW_FIR[i*bl_size + j]);

                    }
                    i_temp++;
                }
                j_temp++;
                i_temp-=(bl_size);
                
            }

        }
            
        else
        {
            for(int down_i = 0; down_i < bl_size; down_i++)
            {
                DTW_DOWN[down_i] = __shfl_down_sync(0xFFFFFFFF, DTW_FIR[(down_i+1)*bl_size - 1], 1, 32);
                
            }
            d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

            DTW_SEC[0] = d + MIN(DTW_DOWN[0],MIN(DTW_FIR[num_per_bl - bl_size],
                                                                         DTW_SEC[num_per_bl - 1]));

            j_temp++;
            
            for(int i = 1;i < bl_size;i++)
            {
                d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

                DTW_SEC[i*bl_size] = d + MIN(DTW_SEC[(i-1)*bl_size],
                                                                     MIN(DTW_DOWN[i-1],DTW_DOWN[i]));

                j_temp++;
            }
            
            i_temp++;
            j_temp-=(bl_size);
            
            for(int i = 1;i < bl_size;i++)
            {
                d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

                DTW_SEC[i] =  d +
                              MIN(DTW_SEC[i - 1], MIN(DTW_FIR[num_per_bl - bl_size + i],DTW_FIR[num_per_bl - bl_size + i-1]));

                i_temp++;
            }
            j_temp++;
            i_temp-=(bl_size - 1);
            
            for(int i = 1;i < bl_size;i++)
                
            {
                for(int j = 1;j < bl_size;j++)
                    
                {
                    d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;

                    DTW_SEC[i*bl_size + j] =  d + MIN(DTW_SEC[i*bl_size + j - 1],
                                                                              MIN(DTW_SEC[(i-1)*bl_size + j - 1],DTW_SEC[(i-1)*bl_size + j]));

                    i_temp++;
                }
                j_temp++;
                i_temp-=(bl_size - 1);
                
            }
            j_temp-=bl_size;
            i_temp +=(bl_size - 1);

        }
    }
    }
}

__device__ void
DTW_stairs_without_shared(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold_2, int w, const FLOAT cb[]) {

    int tid = threadIdx.x%THREAD_NUM_PER_WARP;
    
    int num_tid = THREAD_NUM_PER_WARP;
    
    FLOAT *q =cQuery;
    FLOAT *t = Subject;

    FLOAT DTW_FIR;
    
    if (tid == 16) {
        DTW_FIR = 0;  
    } else {
        DTW_FIR = INFINITY;
    }
    FLOAT DTW_SEC = INFINITY;

    FLOAT DTW_UP = 0;
    FLOAT DTW_DOWN = 0;
    
    int row_bias = 16 - tid;
    int col_bias = tid - 16;

    int i_temp = row_bias;
    int j_temp = col_bias;
    size_t vote;
    
    bool flag_pruning = 0;

    FLOAT d;
    int switch_for_stair = 0;
    bool mask;
    for(int i = 0;i < w;i++){
        mask = switch_for_stair%2;
        switch_for_stair++;
        
        if(!mask)
        {
            DTW_UP = __shfl_up_sync(0xFFFFFFFF, DTW_SEC, 1, 32);
            d = ((i_temp >= 0) && (j_temp >= 0)) ? DIST(q[i_temp],t[j_temp]) : INFINITY;
            DTW_FIR = d + MIN(DTW_FIR,MIN(DTW_SEC,DTW_UP));

            j_temp++;
        }
        
        else
        {
            DTW_DOWN = __shfl_down_sync(0xFFFFFFFF, DTW_FIR, 1, 32);
            
            d = ((i_temp >= 0) && (j_temp >= 0)) ? DIST(q[i_temp],t[j_temp]) : INFINITY;
            DTW_SEC = d + MIN(DTW_FIR,MIN(DTW_SEC,DTW_DOWN));

            i_temp++;
        }
    }

    int cb_index ;
    FLOAT cb_temp;
    for(int i = w;i < 2*m - 1 - w;i++){

        mask = switch_for_stair%2;
        switch_for_stair++;
        
        if(!mask)
        {
            DTW_UP = __shfl_up_sync(0xFFFFFFFF, DTW_SEC, 1, 32);
            
            d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY  :DIST(q[i_temp],t[j_temp]);

            if(i_temp - j_temp ==  w)
            {
                DTW_FIR = d + MIN(DTW_FIR,DTW_SEC);
            }
            else
            {
                DTW_FIR = d + MIN(DTW_FIR,MIN(DTW_SEC,DTW_UP));
            }

            bool ifcb = (tid == 16) && (j_temp == m/2) ;
            ifcb = __shfl_sync(0x1F, ifcb, 16);
            if(ifcb)
            {
                cb_index =(w + j_temp)/CB_LEN+1;
                cb_temp =cb[cb_index];
                vote = __ballot_sync(0xFFFFFFFF, DTW_FIR > threshold_2 - cb_temp);
                if (vote == 0xFFFFFFFF) {
                    if(tid == 16)
                    {
                        Dist = INFINITY;

                    }
                    return;
                    flag_pruning = true;
                }
            }

            j_temp++;

        }
            
        else
        {
            DTW_DOWN = __shfl_down_sync(0xFFFFFFFF, DTW_FIR, 1, 32);
            
            d = (i_temp - j_temp > w ||  j_temp - i_temp> w) ? INFINITY  :DIST(q[i_temp],t[j_temp]);

            if(j_temp - i_temp ==  w)
            {
                DTW_SEC = d + MIN(DTW_FIR,DTW_SEC);
            }
            else
            {
                DTW_SEC = d + MIN(DTW_FIR,MIN(DTW_SEC,DTW_DOWN));
            }

            bool ifcb = (tid == 16) && (j_temp == m/2);
            ifcb = __shfl_sync(0x1F, ifcb, 16);
            if(ifcb)
            {
                cb_index =(w + j_temp)/CB_LEN+1;
                cb_temp =cb[cb_index];
                vote = __ballot_sync(0xFFFFFFFF, DTW_SEC > threshold_2 - cb_temp);
                if (vote == 0xFFFFFFFF) {
                    if(tid == 16)
                    {
                        Dist = INFINITY;

                    }

                    return;
                    flag_pruning = true;
                }
            }

            i_temp++;

        }

    }

    if(!flag_pruning)
    {
        for(int i = 2*m  - 1- w;i < 2*m - 1;i++){
            mask = switch_for_stair%2;
            switch_for_stair++;
            
            if(!mask)
            {
                DTW_UP = __shfl_up_sync(0xFFFFFFFF, DTW_SEC, 1, 32);
                
                d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                DTW_FIR = d + MIN(DTW_FIR,MIN(DTW_SEC,DTW_UP));

                j_temp++;
            }
                
            else
            {
                DTW_DOWN = __shfl_down_sync(0xFFFFFFFF, DTW_FIR, 1, 32);
                
                d = (i_temp < m && j_temp < m) ?  DIST(q[i_temp],t[j_temp]) : INFINITY;
                DTW_SEC = d + MIN(DTW_FIR,MIN(DTW_SEC,DTW_DOWN));

                i_temp++;
            }

        }

        if(tid == 16)
        {
            Dist = sqrt(DTW_FIR);

        }
    }

}

