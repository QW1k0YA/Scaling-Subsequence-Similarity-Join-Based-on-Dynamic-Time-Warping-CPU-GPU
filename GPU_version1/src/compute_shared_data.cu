
#include "../alldef/matrix.cuh"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;

int FIND_POS_LOCAL_f(FLOAT input,FLOAT max_abs_value,int length)
{
    int pos;
    if(input >  3)
    {
        pos = 300 + (input - 3)*10;
    }
    else if(input < -3)
    {
        pos =  - 300  + (input + 3)*10;
    }
    else
    {
        pos = input * 100;
    }

    pos+= length/2;
    return  pos;
}

void
compute_shared_data_local(const vector<FLOAT > &ts, int subseqlen, const vector<vector<FLOAT >> &subs,
                          const vector<vector<FLOAT >> &UU, const vector<vector<FLOAT >> &LL, FLOAT  &real_min,
                          FLOAT  &real_max, vector<FLOAT> &pos_UU, vector<FLOAT> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, FLOAT MAX_REAL_VALUE)
{

    int subcount = ts.size() - subseqlen + 1;
    auto start_time = std::chrono::high_resolution_clock::now();
    FLOAT  max_possible_value = floor(sqrt(subseqlen));

    vector<FLOAT > cnt(ts.size(), 0.0);
    vector<FLOAT > strange_count(ts.size(), 0.0);
    int pp;
    FLOAT  t;
    size_t ts_size = ts.size();
    for(int j = 0;j <subcount;j++)
    {

        for(int i = 0;i < subseqlen;i++)
        {
            
            FLOAT  value = subs[j][i];
            
            pp=FIND_POS_LOCAL_f(value,MAX_REAL_VALUE,len_of_cdf);
            
            if(pp < 0 )
            {
                pp = 0;
            }
            if(pp >= len_of_cdf)
            {
                pp = len_of_cdf - 1;
            }

            count_table_cdf[i+j][pp] = count_table_cdf[i+j][pp] + 1;
            
        }
    }

    for(int i = 0;i < ts_size;i++)
    {
        int sum1 = 0;
        for(int j = 0;j < len_of_cdf;j++)
        {
            sum1 = sum1 + count_table_cdf[i][j];
            count_table_cdf[i][j] = sum1;
        }
    }

    vector<FLOAT > sum_UU(ts.size(), 0.0);
    vector<FLOAT > sum_LL(ts.size(), 0.0);

    for(int i = 0;i < subcount;i++)
    {
        for(int j = 0;j < subseqlen;j++)
        {

            sum_UU[i+j] = sum_UU[i+j]+UU[i][j];
            sum_LL[i+j] = sum_LL[i+j]+LL[i][j];

            cnt[i+j] = cnt[i+j] + 1;
        }
    }
    sum_UU = elementWiseDivison_vv(sum_UU,cnt);
    sum_LL = elementWiseDivison_vv(sum_LL,cnt);

    FLOAT  posUU,posLL;

    for(int i = 0;i < ts_size;i++)
    {
        FLOAT  LL_temp = sum_LL[i];
        FLOAT  UU_temp = sum_UU[i];

        posLL = FIND_POS_LOCAL_f(LL_temp,MAX_REAL_VALUE,len_of_cdf);
        posUU = FIND_POS_LOCAL_f(UU_temp,MAX_REAL_VALUE,len_of_cdf);
        if(posLL < 0)
        {
            posLL = 1;
        }
        if(posLL > len_of_cdf)
        {
            posLL =  len_of_cdf;
        }
        if(posUU < 0)
        {
            posLL = 1;
        }
        if(posLL >  len_of_cdf )
        {
            posLL =  len_of_cdf;
        }
        if(posUU < 0)
        {
            posUU = 1;
        }
        if(posUU >  len_of_cdf )
        {
            posUU = len_of_cdf;
        }
        pos_UU[i] = posUU;
        pos_LL[i] = posLL;
    }
    printf("computing the Shared DATA: V1, LEN_OF_TABLE_LOCAL is %d \n", len_of_cdf);

}

void
compute_shared_data_local(const vector<FLOAT > &ts, int subseqlen, FLOAT **subs,
                          FLOAT **UU, FLOAT **LL, FLOAT  &real_min,
                          FLOAT  &real_max, vector<FLOAT> &pos_UU, vector<FLOAT> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, FLOAT MAX_REAL_VALUE)
{

    int subcount = ts.size() - subseqlen + 1;
    auto start_time = std::chrono::high_resolution_clock::now();
    FLOAT  max_possible_value = floor(sqrt(subseqlen));

    vector<FLOAT > cnt(ts.size(), 0.0);
    vector<FLOAT > strange_count(ts.size(), 0.0);
    int pp;
    FLOAT  t;
    size_t ts_size = ts.size();
    for(int j = 0;j <subcount;j++)
    {

        for(int i = 0;i < subseqlen;i++)
        {
            
            FLOAT  value = subs[j][i];
            
            pp=FIND_POS_LOCAL_f(value,MAX_REAL_VALUE,len_of_cdf);
            
            if(pp < 0 )
            {
                pp = 0;
            }
            if(pp >= len_of_cdf)
            {
                pp = len_of_cdf - 1;
            }

            count_table_cdf[i+j][pp] = count_table_cdf[i+j][pp] + 1;
            
        }
    }

    for(int i = 0;i < ts_size;i++)
    {
        int sum1 = 0;
        for(int j = 0;j < len_of_cdf;j++)
        {
            sum1 = sum1 + count_table_cdf[i][j];
            count_table_cdf[i][j] = sum1;
        }
    }

    vector<FLOAT > sum_UU(ts.size(), 0.0);
    vector<FLOAT > sum_LL(ts.size(), 0.0);

    for(int i = 0;i < subcount;i++)
    {
        for(int j = 0;j < subseqlen;j++)
        {

            sum_UU[i+j] = sum_UU[i+j]+UU[i][j];
            sum_LL[i+j] = sum_LL[i+j]+LL[i][j];

            cnt[i+j] = cnt[i+j] + 1;
        }
    }
    sum_UU = elementWiseDivison_vv(sum_UU,cnt);
    sum_LL = elementWiseDivison_vv(sum_LL,cnt);

    FLOAT  posUU,posLL;

    for(int i = 0;i < ts_size;i++)
    {
        FLOAT  LL_temp = sum_LL[i];
        FLOAT  UU_temp = sum_UU[i];

        posLL = FIND_POS_LOCAL_f(LL_temp,MAX_REAL_VALUE,len_of_cdf);
        posUU = FIND_POS_LOCAL_f(UU_temp,MAX_REAL_VALUE,len_of_cdf);
        if(posLL < 0)
        {
            posLL = 1;
        }
        if(posLL > len_of_cdf)
        {
            posLL =  len_of_cdf;
        }
        if(posUU < 0)
        {
            posLL = 1;
        }
        if(posLL >  len_of_cdf )
        {
            posLL =  len_of_cdf;
        }
        if(posUU < 0)
        {
            posUU = 1;
        }
        if(posUU >  len_of_cdf )
        {
            posUU = len_of_cdf;
        }
        pos_UU[i] = posUU;
        pos_LL[i] = posLL;
    }
    printf("computing the Shared DATA: V1, LEN_OF_TABLE_LOCAL is %d \n", len_of_cdf);

}

