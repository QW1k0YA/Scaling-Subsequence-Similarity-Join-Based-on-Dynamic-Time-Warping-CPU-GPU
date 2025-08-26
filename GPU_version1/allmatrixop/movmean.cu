
#include <iostream>
#include <vector>
#include "../alldef/matrix.cuh"

#include "cmath"
using namespace std;

void mvmean(const std::vector<FLOAT > &a, int l, vector<FLOAT > &miu, vector<FLOAT > &si)
{
    size_t len_a = a.size();
    miu.reserve(len_a);
    si.reserve(len_a);
    FLOAT sum1 = 0;
    FLOAT sum2 = 0;
    for(int i = 0;i < l;i++)
    {
        sum1+= a[i];
        sum2+= a[i]*a[i];
    }

    long long ll = len_a - l +1;
    for(int i = 0;i < ll - 1;i++)
    {
        miu[i] = sum1/l;
        si[i] = sqrt(sum2/l-miu[i]*miu[i]);
        sum1 = sum1 - a[i] + a[i+l];
        sum2 = sum2 - a[i]*a[i] + a[i+l]*a[i+l];
    }
    miu[ll-1] = sum1/l;
    si[ll-1] = sqrt(sum2/l-miu[ll-1]*miu[ll-1]);

}
void mvmean_miu(const std::vector<FLOAT >& a, int len_a, int l, vector<FLOAT >& miu)
{

    FLOAT sum1 = 0;
    FLOAT sum2 = 0;
    for(int i = 0;i < l;i++)
    {
        sum1+= a[i];
    }

    long long ll = len_a - l +1;
    for(int i = 0;i < ll;i++)
    {
        miu[i] = sum1/l;
        sum1 = sum1 - a[i] + a[i+l];
    }

}

std::vector<FLOAT > movmean(const std::vector<FLOAT >& ts, int a, int b, bool c) {

    size_t ts_size = ts.size();
    if(c)
    {
        std::vector<FLOAT > result;

        if (ts.empty() || a < 0 || b < 0 || a + b == 0) {
            std::cerr << "Invalid input parameters. in movmean" << std::endl;
            return result;
        }

        for (size_t i = a; i < ts_size - b; ++i) {
            FLOAT  sum = 0.0;
            for (int j = -a; j <= b; ++j) {
                sum += ts[i + j];
            }
            result.push_back(sum / (a + b + 1));
        }

        return result;

    }
    else
    {
        std::vector<FLOAT > result;

        if (ts.empty() || a < 0 || b < 0 || a + b == 0) {
            std::cerr << "Invalid input parameters.in movmean" << std::endl;
            return result;
        }

        for (size_t i = 0; i < ts_size; ++i) {
            FLOAT  sum = 0.0;
            int count = 0;

            for (int j = -a; j <= b; ++j) {
                if (i + j >= 0 && i + j < ts_size) {
                    sum += ts[i + j];
                    count++;
                }
            }

            if (count > 0) {
                result.push_back(sum / count);
            }
            else {
                
                result.push_back(ts[i]);
            }
        }

        return result;
    }
}
