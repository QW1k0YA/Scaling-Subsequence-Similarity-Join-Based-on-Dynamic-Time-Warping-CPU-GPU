
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;

int FIND_POS_LOCAL_f(double input,double max_abs_value,int length)
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
compute_shared_data_local(const vector<DOUBLE> &ts, int subseqlen, const vector<vector<DOUBLE>> &subs,
                          const vector<vector<DOUBLE>> &UU, const vector<vector<DOUBLE>> &LL, DOUBLE &real_min,
                          DOUBLE &real_max, vector<double> &pos_UU, vector<double> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, double MAX_REAL_VALUE)
{

    int subcount = ts.size() - subseqlen + 1;
    auto start_time = std::chrono::high_resolution_clock::now();
    DOUBLE max_possible_value = floor(sqrt(subseqlen));

    vector<DOUBLE> cnt(ts.size(),0.0);
    vector<DOUBLE> strange_count(ts.size(),0.0);
    int pp;
    DOUBLE t;
    size_t ts_size = ts.size();
    for(int j = 0;j <subcount;j++)
    {

        for(int i = 0;i < subseqlen;i++)
        {
            
            DOUBLE value = subs[j][i];
            
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

    vector<DOUBLE> sum_UU(ts.size(),0.0);
    vector<DOUBLE> sum_LL(ts.size(),0.0);

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

    DOUBLE posUU,posLL;

    for(int i = 0;i < ts_size;i++)
    {
        DOUBLE LL_temp = sum_LL[i];
        DOUBLE UU_temp = sum_UU[i];

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

void compute_shared_data_origin(const vector<DOUBLE>& ts,int subseqlen,int warpmax,const vector<DOUBLE> &mu,const vector<DOUBLE>& sig,const vector<vector<DOUBLE>>&subs,RETURN_COM& result)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    DOUBLE max_possible_value = floor(sqrt(subseqlen));

    DOUBLE real_max = maxinmatrix(subs);
    DOUBLE real_min = mininmatrix(subs);
    DOUBLE max_real_value = floor(MAX(abs(real_max), abs(real_min))) * 2 * 100;

    int len_of_table = static_cast<int>(MAX(600.0, max_real_value));
    auto len = ts.size();
    auto subcount = len - subseqlen + 1;

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> UTS = movmax(ts,warpmax,warpmax);
    vector<DOUBLE> LTS = movmin(ts,warpmax,warpmax);

    vector<DOUBLE> cnt(ts.size(),0.0);

    vector<vector<SHORT>> count_table(ts.size(),vector<SHORT>(static_cast<int>(len_of_table),0));
    vector<DOUBLE> strange_count(ts.size(),0.0);
    int pp;
    DOUBLE t;
    size_t ts_size = ts.size();
    for(int i = 1;i <= subseqlen;i++)
    {
        for(int j = 1;j <=subcount;j++)
        {
            DOUBLE value = subs[i-1][j-1];
            pp = static_cast<int>(floor(value*100 + len_of_table/2));
            if(pp <= 0 )
            {
                pp = 1;
            }
            if(pp > len_of_table)
            {
                pp = len_of_table;
            }

            t = MAX(abs(value * 100) - 200, 0.0);
            strange_count[i+j-2] = strange_count[i+j-2] + t/100;
            count_table[i+j-2][pp-1] = count_table[i+j-2][pp-1] + 1.0;

        }
    }

    for(int i = 1;i <= ts_size;i++)
    {
        int sum1 = 0;
        for(int j = 1;j <=len_of_table;j++)
        {
            sum1 = sum1 + count_table[i-1][j-1];
            count_table[i-1][j-1] = sum1;
        }
    }

    for(int i = 1;i <= ts_size;i++)
    {
        strange_count[i-1] = strange_count[i-1]/count_table[i-1][len_of_table-1];
    }

    vector<bool> attached_mask(strange_count.size());

    int i = 0;
    size_t strange_count_size = strange_count.size();
    for(i = 0;i < strange_count_size;i++)
    {
        if(strange_count[i] >= 1)
        {
            attached_mask[i] = true;
        }
        else
        {
            attached_mask[i] = false;
        }
        i++;
    }

    vector<DOUBLE> sum_UU(ts.size(),0.0);
    vector<DOUBLE> sum_LL(ts.size(),0.0);
    DOUBLE LL,UU;
    for(int i = 1;i <= subcount;i++)
    {
        for(int j = 1;j <= subseqlen;j++)
        {
            LL = (LTS[i+j-2] - mu[i-1])*invsig[i-1];
            UU = (UTS[i+j-2] - mu[i-1])*invsig[i-1];

            sum_UU[i+j-2] = sum_UU[i+j-2]+UU;
            sum_LL[i+j-2] = sum_LL[i+j-2]+LL;
            cnt[i+j-2] = cnt[i+j-2] + 1;
        }
    }
    sum_UU = elementWiseDivison_vv(sum_UU,cnt);
    sum_LL = elementWiseDivison_vv(sum_LL,cnt);

    vector<DOUBLE> pos_UU(ts.size(),1.0);
    vector<DOUBLE> pos_LL(ts.size(),1.0);

    DOUBLE posUU;
    DOUBLE posLL;

    for(i = 1;i <= ts_size;i++)
    {
        LL = sum_LL[i-1];
        UU = sum_UU[i-1];
        posLL = floor(LL*100 + len_of_table/2);
        posUU = floor(UU*100 + len_of_table/2);
        if(posLL <= 0)
        {
            posLL = 1;
        }
        if(posLL >= len_of_table)
        {
            posLL = len_of_table;
        }
        if(posUU <= 0)
        {
            posLL = 1;
        }
        if(posLL >= len_of_table)
        {
            posLL = len_of_table;
        }
        if(posUU <= 0)
        {
            posUU = 1;
        }
        if(posUU >= len_of_table)
        {
            posUU = len_of_table;
        }
        pos_UU[i-1] = posUU;
        pos_LL[i-1] = posLL;
    }
    printf("computing the Shared DATA: V1, LEN_OF_TABLE_LOCAL is %d \n", len_of_table);

    result.len_of_table = len_of_table;
    result.pos_UU = pos_UU;
    result.pos_LL = pos_LL;
    result.count_table = count_table;
    result.UTS = UTS;
    result.LTS = LTS;
    result.attached_mask = attached_mask;

}

