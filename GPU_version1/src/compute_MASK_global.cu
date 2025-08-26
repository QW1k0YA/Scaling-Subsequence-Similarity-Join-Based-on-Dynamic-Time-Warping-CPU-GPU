
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../allunder/underdtw.cuh"

using namespace std;
void compute_MASK_global(const vector<FLOAT > &a, int subseqLen, int len_of_table_global,
                         vector<vector<FLOAT >> &count_table_global, const vector<vector<FLOAT >> &subs_U,
                         const vector<vector<FLOAT >> &subs_L, vector<FLOAT > &MASK_global)
{
    int subcount = a.size() - subseqLen + 1;
    int adjust_factor = 0;
    vector<FLOAT > ttt(vector<FLOAT >(static_cast<int>(a.size()) - 2 * subseqLen, subseqLen));
    vector<FLOAT > aa_global(subseqLen);
    for(int i = 0; i < subseqLen; i++)
    {
        aa_global[i] = i + 1;
    }

    vector<FLOAT > bb_global = linspace(subseqLen, 1, subseqLen);
    vector<FLOAT > mod = mergeVectors(aa_global, ttt, bb_global);

    vector<FLOAT > normal_DIFF(a.size(), 0.0);
    vector<FLOAT > cnt(a.size(), 0.0);
    FLOAT  UU,LL,posLL,posUU,diff2;
    FLOAT pad=0.0; 
    int i;
    FLOAT  value;
    int pp;
    for(i = 0;i < subcount;i++)
    {
        for(int j =0; j < subseqLen; j++)
        {
            LL = subs_L[i][j]-pad;
            UU = subs_U[i][j]+pad;

            posLL = FIND_POS_GLOBAL(LL,len_of_table_global);
            posUU = FIND_POS_GLOBAL(UU,len_of_table_global);

            if(posLL <= 0 )
            {
                posLL = 1;
            }
            if(posLL > len_of_table_global)
            {
                posLL = len_of_table_global;
            }
            if(posUU<=0)
            {
                posUU = 1;
            }
            if(posUU > len_of_table_global)
            {
                posUU = len_of_table_global;
            }

            diff2 = count_table_global[j][posUU - 1] - count_table_global[j][posLL - 1];

            normal_DIFF[i+j] += diff2 / subcount;
            cnt[i+j] += 1;
        }
    }

    normal_DIFF = elementWiseDivison_vv(normal_DIFF,mod);
    vector<FLOAT > DIFF = normal_DIFF;

    FLOAT  avg_diff,std_diff,ths_diff;
    avg_diff = Mean(DIFF);
    std_diff = stddev(DIFF);
    ths_diff = avg_diff + adjust_factor*std_diff;

    i = 0;

    for(auto value:DIFF)
    {
        if(value <= ths_diff)
        {
            MASK_global[i++] = 1;
        }
        else
        {
            MASK_global[i++] = 0;
        }
    }

}

void compute_MASK_global(const vector<FLOAT > &a, int subseqLen, int len_of_table_global,
                         vector<vector<FLOAT >> &count_table_global, FLOAT **subs_U,
                         FLOAT **subs_L, vector<FLOAT > &MASK_global)
{
    int subcount = a.size() - subseqLen + 1;
    int adjust_factor = 0;
    vector<FLOAT > ttt(vector<FLOAT >(static_cast<int>(a.size()) - 2 * subseqLen, subseqLen));
    vector<FLOAT > aa_global(subseqLen);
    for(int i = 0; i < subseqLen; i++)
    {
        aa_global[i] = i + 1;
    }

    vector<FLOAT > bb_global = linspace(subseqLen, 1, subseqLen);
    vector<FLOAT > mod = mergeVectors(aa_global, ttt, bb_global);

    vector<FLOAT > normal_DIFF(a.size(), 0.0);
    vector<FLOAT > cnt(a.size(), 0.0);
    FLOAT  UU,LL,posLL,posUU,diff2;
    FLOAT pad=0.0; 
    int i;
    FLOAT  value;
    int pp;
    for(i = 0;i < subcount;i++)
    {
        for(int j =0; j < subseqLen; j++)
        {
            LL = subs_L[i][j]-pad;
            UU = subs_U[i][j]+pad;

            posLL = FIND_POS_GLOBAL(LL,len_of_table_global);
            posUU = FIND_POS_GLOBAL(UU,len_of_table_global);

            if(posLL <= 0 )
            {
                posLL = 1;
            }
            if(posLL > len_of_table_global)
            {
                posLL = len_of_table_global;
            }
            if(posUU<=0)
            {
                posUU = 1;
            }
            if(posUU > len_of_table_global)
            {
                posUU = len_of_table_global;
            }

            diff2 = count_table_global[j][posUU - 1] - count_table_global[j][posLL - 1];

            normal_DIFF[i+j] += diff2 / subcount;
            cnt[i+j] += 1;
        }
    }

    normal_DIFF = elementWiseDivison_vv(normal_DIFF,mod);
    vector<FLOAT > DIFF = normal_DIFF;

    FLOAT  avg_diff,std_diff,ths_diff;
    avg_diff = Mean(DIFF);
    std_diff = stddev(DIFF);
    ths_diff = avg_diff + adjust_factor*std_diff;

    i = 0;

    for(auto value:DIFF)
    {
        if(value <= ths_diff)
        {
            MASK_global[i++] = 1;
        }
        else
        {
            MASK_global[i++] = 0;
        }
    }

}
