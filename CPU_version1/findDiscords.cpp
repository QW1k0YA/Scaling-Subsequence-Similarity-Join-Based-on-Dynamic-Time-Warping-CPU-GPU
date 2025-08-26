
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
#include "alldef/typedefdouble.h"
#include "allunder/undermpx_v2.h"
using namespace std;

vector<DOUBLE> findDiscords(const vector<DOUBLE>& matrixProfile,int exclusionLen)
{
    int discordCount = 3;
    vector<DOUBLE> discordIdx(discordCount,0);
    vector<DOUBLE> idx = sortInd(matrixProfile);

    vector<DOUBLE> matrixProfile_ = matrixProfile;

    DOUBLE f;
    DOUBLE exclRangeBegin, exclRangeEnd;

    for(int i = 1; i <= discordCount; i++)
    {
        
        vector<DOUBLE> temp1;
        vector<bool> temp2,temp3,temp4;
        for(auto value:idx)
        {
            temp1.push_back(matrixProfile_[value - 1]);
        }

        temp2 = isinfinite(temp1);

        for(auto value:temp1)
        {
            if(value > -1)
            {
                temp3.push_back(1);
            }
            else
            {
                temp3.push_back(0);
            }
        }

        for(int k = 1;k <= temp2.size();k++)
        {
            if(!temp2[k-1] && temp3[k-1])
            {
                temp4.push_back(1);
            }
            else
            {
                temp4.push_back(0);
            }
        }

        if(!vectorisempty(findNonZero(temp4)))
        {
            f = findNonZero(temp4)[0];
        }

        if(vectorisempty(findNonZero(temp4)) || matrixProfile_[f-1] == -1)
        {
            for(int m = i;m <= discordIdx.size();m++)
            {
                discordIdx[i-1] = NAN;
            }
            break;
        }

        discordIdx[i-1] = idx[f-1];
        exclRangeBegin = max(1.0,discordIdx[i-1] - exclusionLen + 1);
        exclRangeEnd = MIN(DOUBLE (matrixProfile_.size()), discordIdx[i - 1] + exclusionLen - 1);
        for(int n = exclRangeBegin;n <= exclRangeEnd;n++)
        {
            matrixProfile_[n-1] = NAN;
        }
    }

    return discordIdx;

}

