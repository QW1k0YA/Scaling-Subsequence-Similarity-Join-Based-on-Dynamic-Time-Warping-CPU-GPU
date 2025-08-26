
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allunder/underdtw.h"

using namespace std;
void compute_MASK_global(const vector<DOUBLE> &a, int subseqLen,const vector<vector<DOUBLE>> &my_subs,
                         int len_of_table_global, vector<vector<DOUBLE>> &count_table_global,
                         const vector<vector<DOUBLE>> &subs_U,const vector<vector<DOUBLE>> &subs_L,
                         vector<DOUBLE>& MASK_global)
{
    int subcount = a.size() - subseqLen + 1;
    int adjust_factor = 0;
    vector<DOUBLE> ttt(vector<DOUBLE>(static_cast<int>(a.size()) - 2 * subseqLen, subseqLen));
    vector<DOUBLE> aa_global(subseqLen);
    for(int i = 0; i < subseqLen; i++)
    {
        aa_global[i] = i + 1;
    }

    vector<DOUBLE> bb_global = linspace(subseqLen, 1, subseqLen);
    vector<DOUBLE> mod = mergeVectors(aa_global, ttt, bb_global);
    
    vector<DOUBLE> normal_DIFF(a.size(), 0.0);
    vector<DOUBLE> cnt(a.size(), 0.0);
    DOUBLE UU,LL,posLL,posUU,diff2;
    double pad=0.0; 
    int i;
    DOUBLE value;
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
    vector<DOUBLE> DIFF = normal_DIFF;

    DOUBLE avg_diff,std_diff,ths_diff;
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