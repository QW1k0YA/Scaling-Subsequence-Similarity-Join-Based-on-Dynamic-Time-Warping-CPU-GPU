
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allunder/underfast.h"
#include "alldef/allstruct.h"
using namespace std;

void compute_shared_dataV3(const vector<DOUBLE>& ts,int subseqlen,int warpmax,const vector<DOUBLE> &mu,const vector<DOUBLE>& sig,const vector<vector<DOUBLE>>&subs,RETURN_COM& result)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    DOUBLE max_possible_value = ceil(sqrt(subseqlen));
    DOUBLE real_max = maxinmatrix(subs);
    DOUBLE real_min = mininmatrix(subs);
    DOUBLE max_abs = ceil(MAX(abs(real_max), abs(real_min) + 0.1));

    DOUBLE max_real_value,len_of_table,base_len,extend_len;
    if(max_abs <= 5)
    {
        max_real_value = max_abs*2*100;
        len_of_table = MAX(600.0, max_real_value);
        base_len = len_of_table;
        extend_len = 0;
    }
    else
    {
        extend_len = (max_abs-5)*2*10;
        base_len= 1000 ;
        len_of_table=base_len+extend_len;
    }

    vector<vector<SHORT>> count_table(ts.size(),vector<SHORT>(len_of_table,0.0));
    auto len = ts.size();
    auto subcount = len - subseqlen +1;

    vector<vector<DOUBLE>> sig_row = convertTo2DRowOrder(sig);
    vector<vector<DOUBLE>> invsig_m = elementWiseDivision(1.0, sig_row);
    vector<DOUBLE> invsig = extractVector(invsig_m, 0, 1);

    vector<DOUBLE> UTS = movmax(ts,warpmax,warpmax);
    vector<DOUBLE> LTS = movmin(ts,warpmax,warpmax);

    vector<DOUBLE> cnt(ts.size(),0.0);

    DOUBLE value;
    int pp;
    for(int i = 1;i<= subseqlen;i++)
    {
        for(int j = 1;j <=subcount;j++)
        {
            value = subs[i-1][j-1];
            if(abs(value)<=5)
            {
                pp=floor(value*100+ len_of_table/2);
            }
            else
            {
                if(value < 0)
                {
                    pp=floor(value*10)-450+ len_of_table/2;
                }
                else
                {
                    pp=floor(value*10)+450+ len_of_table/2;
                }
            }
            count_table[i+j-2][pp-1] = count_table[i+j-2][pp-1]+1;
            if((pp<0)||(pp>len_of_table))
            {
                cerr<<"wrong index IN compute_shared_dataV3";
                if(pp<0)
                {
                    pp =1;
                }
                if(pp>=len_of_table)
                {
                    pp = static_cast<int>(len_of_table);
                }
            }
        }
    }

    for(int i = 1;i<=len;i++)
    {
        for(int j = 2;j <=len_of_table;j++)
        {
            count_table[i-1][j-1] = count_table[i-1][j-1] + count_table[i-1][j-2];
        }
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

    vector<DOUBLE> pos_UU(ts.size(),0.0);
    vector<DOUBLE> pos_LL(ts.size(),0.0);

    DOUBLE posLL;
    DOUBLE posUU;

    auto ts_size = ts.size();
    for(int i = 1;i <= ts_size;i++)
    {
        LL = sum_LL[i-1];
        UU = sum_UU[i-1];

        if(abs(LL)<=5)
        {
            posLL=floor(LL*100+ len_of_table/2);;
        }
        else
        {
            if(LL<0)
            {
                posLL=floor(LL*10)-450+ len_of_table/2;
            }
            else
            {
                posLL=floor(LL*10)+450+ len_of_table/2;
            }
        }

        if(abs(UU)<=5)
        {
            posUU=floor(UU*100+ len_of_table/2);
        }
        else
        {
            if(UU<0)
            {
                posUU=floor(UU*10)-450+ len_of_table/2;
            }
            else
            {
                posUU=floor(UU*10)+450+ len_of_table/2;
            }
        }

        if(posLL <= 0 )
        {
            posLL = 1;
        }
        if(posLL > len_of_table)
        {
            posLL = len_of_table;
        }
        if(posUU<=0)
        {
            posUU = 1;
        }
        if(posUU>len_of_table)
        {
            posUU = len_of_table;
        }
        pos_UU[i-1] = posUU;
        pos_LL[i-1] = posLL;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    printf("computing the Shared DATA: V3, LEN_OF_TABLE is %f, cost: %5.3f \n", len_of_table,DOUBLE(duration.count()/1000000));

    result.len_of_table = len_of_table;
    result.pos_UU = pos_UU;
    result.pos_LL = pos_LL;
    result.count_table = count_table;
    result.UTS = UTS;
    result.LTS = LTS;

}

