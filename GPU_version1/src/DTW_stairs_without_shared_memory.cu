
#include <cuda_runtime.h>
#include "GPU_parameters.h"
#include "matrix.cuh"
#define REGISTER_NUM 4
#define THREAD_NUM_PER_WARP 32

__device__ void
DTW_stairs_for_block_without_shared(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold, int w, FLOAT *q, FLOAT *t,int bl_size) {
    
    int num_per_bl = bl_size*bl_size;
    int tid = threadIdx.x;
    
    int num_tid = 32;
    
    q = cQuery;
    t = Subject;
    
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

    int switch_for_stair = 0;
    bool mask;
    FLOAT d;

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
                if(i_temp <0 || j_temp < 0)
                {
                    printf("fuck in 482 i_temp = %d,j_temp = %d,tid = %d,step = %d\n",i_temp,j_temp,tid,step);
                }

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
                    if( i_temp <0 || j_temp < 0)
                    {
                        printf("fuck in 499 i_temp = %d,j_temp = %d,tid = %d,step = %d\n",i_temp,j_temp,tid,step);
                    }

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

__device__ void
DTW_stairs_without_shared(FLOAT *Subject, FLOAT *cQuery, FLOAT &Dist, int m, FLOAT threshold_2, int w, FLOAT *q, FLOAT *t) {

    int tid = threadIdx.x%THREAD_NUM_PER_WARP;
    
    int num_tid = THREAD_NUM_PER_WARP;
    
    q =cQuery;
    t = Subject;

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
    FLOAT cb_temp1;
    FLOAT cb_temp2;
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

