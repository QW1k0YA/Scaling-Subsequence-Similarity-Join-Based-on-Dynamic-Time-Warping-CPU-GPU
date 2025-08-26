
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../allunder/underdtw.cuh"

using namespace std;

void compute_shared_data_global(const vector<FLOAT > &ts, int subseqLen, const vector<vector<FLOAT >> &my_subs,
                                int len_of_table_global, vector<vector<FLOAT >> &count_table_global)
{
    int subcount = ts.size() - subseqLen + 1;

    int i;
    FLOAT  value;
    int pp;
    for(int j = 0; j < subcount; j++)
    {
        for(i = 0; i < subseqLen; i++) {
            value = my_subs[j][i];

            pp=FIND_POS_GLOBAL(value,len_of_table_global);
            if (pp <= 0) {
                pp = 1;
            }
            if (pp > len_of_table_global) {
                pp = len_of_table_global;
            }
            count_table_global[i][pp - 1] = count_table_global[i][pp - 1] + 1;
        }
    }
    FLOAT  sum1;
    for(i = 0; i < subseqLen; i++)
    {
        sum1 = 0;
        for(int j = 0; j < len_of_table_global; j++)
        {
            sum1 += count_table_global[i][j];
            count_table_global[i][j] = sum1;
        }
    }

}

void compute_shared_data_global(const vector<FLOAT > &ts, int subseqLen, FLOAT* * &my_subs,
                                int len_of_table_global, vector<vector<FLOAT >> &count_table_global)
{
    int subcount = ts.size() - subseqLen + 1;

    int i;
    FLOAT  value;
    int pp;
    for(int j = 0; j < subcount; j++)
    {
        for(i = 0; i < subseqLen; i++) {
            value = my_subs[j][i];

            pp=FIND_POS_GLOBAL(value,len_of_table_global);
            if (pp <= 0) {
                pp = 1;
            }
            if (pp > len_of_table_global) {
                pp = len_of_table_global;
            }
            count_table_global[i][pp - 1] = count_table_global[i][pp - 1] + 1;
        }
    }
    FLOAT  sum1;
    for(i = 0; i < subseqLen; i++)
    {
        sum1 = 0;
        for(int j = 0; j < len_of_table_global; j++)
        {
            sum1 += count_table_global[i][j];
            count_table_global[i][j] = sum1;
        }
    }

}

