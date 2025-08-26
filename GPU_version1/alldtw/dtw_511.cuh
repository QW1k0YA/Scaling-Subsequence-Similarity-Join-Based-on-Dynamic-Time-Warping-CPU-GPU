
#ifndef SUB_WARP_FULLDTW_ECG
#define SUB_WARP_FULLDTW_ECG

__global__ void dtw_511(FLOAT *Subject, FLOAT *Dist, int num_features, FLOAT *cQuery) {

    const int blid = blockIdx.x;
    const int thid = threadIdx.x;
    const int lane = num_features+1;
    const int base = blid*num_features;
    const int WARP_SIZE = 32;
    const int l = thid;

    FLOAT  penalty_left = INFINITY;
    FLOAT  penalty_diag = 0;  
    FLOAT  penalty_here0 = INFINITY; 
    FLOAT  penalty_here1 = INFINITY; 
    FLOAT  penalty_here2 = INFINITY; 
    FLOAT  penalty_here3 = INFINITY; 
    FLOAT  penalty_here4 = INFINITY; 
    FLOAT  penalty_here5 = INFINITY; 
    FLOAT  penalty_here6 = INFINITY; 
    FLOAT  penalty_here7 = INFINITY; 
    FLOAT  penalty_here8 = INFINITY; 
    FLOAT  penalty_here9 = INFINITY; 
    FLOAT  penalty_here10 = INFINITY; 
    FLOAT  penalty_here11 = INFINITY; 
    FLOAT  penalty_here12 = INFINITY; 
    FLOAT  penalty_here13 = INFINITY; 
    FLOAT  penalty_here14 = INFINITY; 
    FLOAT  penalty_here15 = INFINITY; 
    FLOAT  penalty_temp0;
    FLOAT  penalty_temp1;

    if (thid == 0) {
        penalty_left = INFINITY;
        penalty_diag = INFINITY;
        penalty_here0 = INFINITY; 
        penalty_here1 = INFINITY; 
        penalty_here2 = INFINITY; 
        penalty_here3 = INFINITY; 
        penalty_here4 = INFINITY; 
        penalty_here5 = INFINITY; 
        penalty_here6 = INFINITY; 
        penalty_here7 = INFINITY; 
        penalty_here8 = INFINITY; 
        penalty_here9 = INFINITY; 
        penalty_here10 = INFINITY; 
        penalty_here11 = INFINITY; 
        penalty_here12 = INFINITY; 
        penalty_here12 = INFINITY; 
        penalty_here13 = INFINITY; 
        penalty_here14 = INFINITY; 
        penalty_here15 = INFINITY; 
    }

    const FLOAT  subject_value0 = l == 0 ? 0 : Subject[base + 16 * l - 1];

    const FLOAT  subject_value1 = Subject[base + 16 * l];
    const FLOAT  subject_value2 = Subject[base + 16 * l + 1];
    const FLOAT  subject_value3 = Subject[base + 16 * l + 2];
    const FLOAT  subject_value4 = Subject[base + 16 * l + 3];
    const FLOAT  subject_value5 = Subject[base + 16 * l + 4];
    const FLOAT  subject_value6 = Subject[base + 16 * l + 5];
    const FLOAT  subject_value7 = Subject[base + 16 * l + 6];
    const FLOAT  subject_value8 = Subject[base + 16 * l + 7];
    const FLOAT  subject_value9 = Subject[base + 16 * l + 8];
    const FLOAT  subject_value10 = Subject[base + 16 * l + 9];
    const FLOAT  subject_value11 = Subject[base + 16 * l + 10];
    const FLOAT  subject_value12 = Subject[base + 16 * l + 11];
    const FLOAT  subject_value13 = Subject[base + 16 * l + 12];
    const FLOAT  subject_value14 = Subject[base + 16 * l + 13];
    const FLOAT  subject_value15 = Subject[base + 16 * l + 14];

    int counter = 1;
    FLOAT  query_value = INFINITY;
    FLOAT  new_query_value = cQuery[thid];
    if (thid == 0) query_value = new_query_value;
    if (thid == 0) penalty_here1 = 0; 
    
    new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
    
    penalty_temp0 = penalty_here0;
    penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
    
    penalty_temp1 = INFINITY;
    penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
    penalty_temp0 = penalty_here2;
    penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
    penalty_temp1 = penalty_here3;
    penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));
    penalty_temp0 = penalty_here4;
    penalty_here4 = (query_value-subject_value4) * (query_value-subject_value4) + min(penalty_here3, min(penalty_here4, penalty_temp1));
    penalty_temp1 = penalty_here5;
    penalty_here5 = (query_value-subject_value5) * (query_value-subject_value5) + min(penalty_here4, min(penalty_here5, penalty_temp0));
    penalty_temp0 = penalty_here6;
    penalty_here6 = (query_value-subject_value6) * (query_value-subject_value6) + min(penalty_here5, min(penalty_here6, penalty_temp1));
    penalty_temp1 = penalty_here7;
    penalty_here7 = (query_value-subject_value7) * (query_value-subject_value7) + min(penalty_here6, min(penalty_here7, penalty_temp0));

    penalty_temp0 = penalty_here8;
    penalty_here8 = (query_value-subject_value8) * (query_value-subject_value8) + min(penalty_here7, min(penalty_here8, penalty_temp1));
    penalty_temp1 = penalty_here9;
    penalty_here9 = (query_value-subject_value9) * (query_value-subject_value9) + min(penalty_here8, min(penalty_here9, penalty_temp0));
    penalty_temp0 = penalty_here10;
    penalty_here10 = (query_value-subject_value10) * (query_value-subject_value10) + min(penalty_here9, min(penalty_here10, penalty_temp1));
    penalty_temp1 = penalty_here11;
    penalty_here11 = (query_value-subject_value11) * (query_value-subject_value11) + min(penalty_here10, min(penalty_here11, penalty_temp0));
    penalty_temp0 = penalty_here12;
    penalty_here12 = (query_value-subject_value12) * (query_value-subject_value12) + min(penalty_here11, min(penalty_here12, penalty_temp1));
    penalty_temp1 = penalty_here13;
    penalty_here13 = (query_value-subject_value13) * (query_value-subject_value13) + min(penalty_here12, min(penalty_here13, penalty_temp0));

    penalty_temp0 = penalty_here14;
    penalty_here14 = (query_value-subject_value14) * (query_value-subject_value14) + min(penalty_here13, min(penalty_here14, penalty_temp1));
    penalty_here15 = (query_value-subject_value15) * (query_value-subject_value15) + min(penalty_here14, min(penalty_here15, penalty_temp0));

    query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
    if (thid == 0) query_value = new_query_value;
    new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
    counter++;

    penalty_diag = penalty_left;
    penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here15, 1, 32);

    if (thid == 0) penalty_left = INFINITY;

    for (int k = 3; k < lane - WARP_SIZE-1; k++) {

        const int i = k-l;
        
        penalty_temp0 = penalty_here0;
        penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
        
        penalty_temp1 = penalty_here1;
        penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
        penalty_temp0 = penalty_here2;
        penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
        penalty_temp1 = penalty_here3;
        penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));
        penalty_temp0 = penalty_here4;
        penalty_here4 = (query_value-subject_value4) * (query_value-subject_value4) + min(penalty_here3, min(penalty_here4, penalty_temp1));
        penalty_temp1 = penalty_here5;
        penalty_here5 = (query_value-subject_value5) * (query_value-subject_value5) + min(penalty_here4, min(penalty_here5, penalty_temp0));
        penalty_temp0 = penalty_here6;
        penalty_here6 = (query_value-subject_value6) * (query_value-subject_value6) + min(penalty_here5, min(penalty_here6, penalty_temp1));
        penalty_temp1 = penalty_here7;
        penalty_here7 = (query_value-subject_value7) * (query_value-subject_value7) + min(penalty_here6, min(penalty_here7, penalty_temp0));

        penalty_temp0 = penalty_here8;
        penalty_here8 = (query_value-subject_value8) * (query_value-subject_value8) + min(penalty_here7, min(penalty_here8, penalty_temp1));
        penalty_temp1 = penalty_here9;
        penalty_here9 = (query_value-subject_value9) * (query_value-subject_value9) + min(penalty_here8, min(penalty_here9, penalty_temp0));
        penalty_temp0 = penalty_here10;
        penalty_here10 = (query_value-subject_value10) * (query_value-subject_value10) + min(penalty_here9, min(penalty_here10, penalty_temp1));
        penalty_temp1 = penalty_here11;
        penalty_here11 = (query_value-subject_value11) * (query_value-subject_value11) + min(penalty_here10, min(penalty_here11, penalty_temp0));
        penalty_temp0 = penalty_here12;
        penalty_here12 = (query_value-subject_value12) * (query_value-subject_value12) + min(penalty_here11, min(penalty_here12, penalty_temp1));
        penalty_temp1 = penalty_here13;
        penalty_here13 = (query_value-subject_value13) * (query_value-subject_value13) + min(penalty_here12, min(penalty_here13, penalty_temp0));

        penalty_temp0 = penalty_here14;
        penalty_here14 = (query_value-subject_value14) * (query_value-subject_value14) + min(penalty_here13, min(penalty_here14, penalty_temp1));
        penalty_here15 = (query_value-subject_value15) * (query_value-subject_value15) + min(penalty_here14, min(penalty_here15, penalty_temp0));

        if (counter%32 == 0) new_query_value = cQuery[i+2*thid-1];
        query_value = __shfl_up_sync(0xFFFFFFFF, query_value, 1, 32);
        if (thid == 0) query_value = new_query_value;
        new_query_value = __shfl_down_sync(0xFFFFFFFF, new_query_value, 1, 32);
        
        counter++;

        penalty_diag = penalty_left;
        penalty_left = __shfl_up_sync(0xFFFFFFFF, penalty_here15, 1, 32);
        
        if (thid == 0) penalty_left = INFINITY;

    }
    penalty_temp0 = penalty_here0;
    penalty_here0 = (query_value-subject_value0) * (query_value-subject_value0) + min(penalty_left, min(penalty_here0, penalty_diag));
    penalty_temp1 = penalty_here1;
    penalty_here1 = (query_value-subject_value1) * (query_value-subject_value1) + min(penalty_here0, min(penalty_here1, penalty_temp0));
    penalty_temp0 = penalty_here2;
    penalty_here2 = (query_value-subject_value2) * (query_value-subject_value2) + min(penalty_here1, min(penalty_here2, penalty_temp1));
    penalty_temp1 = penalty_here3;
    penalty_here3 = (query_value-subject_value3) * (query_value-subject_value3) + min(penalty_here2, min(penalty_here3, penalty_temp0));
    penalty_temp0 = penalty_here4;
    penalty_here4 = (query_value-subject_value4) * (query_value-subject_value4) + min(penalty_here3, min(penalty_here4, penalty_temp1));
    penalty_temp1 = penalty_here5;
    penalty_here5 = (query_value-subject_value5) * (query_value-subject_value5) + min(penalty_here4, min(penalty_here5, penalty_temp0));
    penalty_temp0 = penalty_here6;
    penalty_here6 = (query_value-subject_value6) * (query_value-subject_value6) + min(penalty_here5, min(penalty_here6, penalty_temp1));
    penalty_temp1 = penalty_here7;
    penalty_here7 = (query_value-subject_value7) * (query_value-subject_value7) + min(penalty_here6, min(penalty_here7, penalty_temp0));

    penalty_temp0 = penalty_here8;
    penalty_here8 = (query_value-subject_value8) * (query_value-subject_value8) + min(penalty_here7, min(penalty_here8, penalty_temp1));
    penalty_temp1 = penalty_here9;
    penalty_here9 = (query_value-subject_value9) * (query_value-subject_value9) + min(penalty_here8, min(penalty_here9, penalty_temp0));
    penalty_temp0 = penalty_here10;
    penalty_here10 = (query_value-subject_value10) * (query_value-subject_value10) + min(penalty_here9, min(penalty_here10, penalty_temp1));
    penalty_temp1 = penalty_here11;
    penalty_here11 = (query_value-subject_value11) * (query_value-subject_value11) + min(penalty_here10, min(penalty_here11, penalty_temp0));
    penalty_temp0 = penalty_here12;
    penalty_here12 = (query_value-subject_value12) * (query_value-subject_value12) + min(penalty_here11, min(penalty_here12, penalty_temp1));
    penalty_temp1 = penalty_here13;
    penalty_here13 = (query_value-subject_value13) * (query_value-subject_value13) + min(penalty_here12, min(penalty_here13, penalty_temp0));

    penalty_temp0 = penalty_here14;
    penalty_here14 = (query_value-subject_value14) * (query_value-subject_value14) + min(penalty_here13, min(penalty_here14, penalty_temp1));
    penalty_here15 = (query_value-subject_value15) * (query_value-subject_value15) + min(penalty_here14, min(penalty_here15, penalty_temp0));

    if(thid == blockDim.x-1)
    {
        Dist[blid] = sqrt(penalty_here15);
        if(Dist[blid] < 6)
        printf("%f \n",Dist[blid]);
    }

}
#endif

