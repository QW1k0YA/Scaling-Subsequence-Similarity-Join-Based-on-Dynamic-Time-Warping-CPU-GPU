
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../alldef/typedefdouble.cuh"
#include "../allunder/undermpx_v2.cuh"
#include "numeric"
#include <thread>
#include "../allunder/undermpx_v2.cuh"
#include <functional>

using namespace std;

void findMotifs(const vector<FLOAT >& timeSeries, const vector<FLOAT >& mu, const vector<FLOAT >& invnorm,
                const vector<FLOAT >& matrixProfile, const vector<int> &profileIndex, int subseqLen, int exclusionLen, vector <vector<int>> &result)
{
    FLOAT  motifCount = 3;
    FLOAT  radius = 2;
    FLOAT  neighborCount = 10;
    vector<FLOAT > matrixProfile_ = matrixProfile;

    int padLen;

    std::function< vector<FLOAT >(int)> crosscov_;
    if(subseqLen < thread::hardware_concurrency()*128 )
    {
        auto crosscov = [=](int idx) -> vector<FLOAT >{
            vector<FLOAT > temp1,temp2,temp3,temp4;
            for(int i = idx + subseqLen -1;i >= idx ;i--)
            {
                temp1.push_back(timeSeries[i-1] - mu[idx-1]);
            }

            temp2 = elementwiseMultiply_nv(invnorm[idx-1],temp1);
            return convolve_valid(timeSeries,temp2);
        };
        crosscov_ = crosscov;
    }
    else
    {
        padLen = pow(2, nextpow2(timeSeries.size()));
        auto crosscov = [=](int idx) -> vector<FLOAT >{
            vector<FLOAT > temp1,temp2;
            vector<complex<FLOAT >> temp3,temp4,temp5,temp6;
            for(int i = idx;i <= idx + subseqLen -1 ;i++)
            {
                temp1.push_back(timeSeries[i-1] - mu[idx-1]);
            }
            temp2 = elementwiseMultiply_nv(invnorm[idx-1],temp1);

            temp3 = fft(temp2,padLen);
            temp4 = fft(timeSeries,padLen);

            temp5 = conjugate(temp4);

            temp6 = elementWiseMultiply_complex(temp3,temp4);
            return ifft(temp6, true);

        };
        crosscov_ = crosscov;
    }

    FLOAT  corr = 0;
    int motIdx = 1;
    for(int i = 1;i <= motifCount;i++)
    {

        int j = 1;
        corr = matrixProfile[j-1];

        for(auto value:matrixProfile)
        {
            if(value > corr || (isnan(corr)&&(!isnan(value))))
            {
                corr = value;
                motIdx = j;
            }
            j++;
        }

        vector<FLOAT > corr_(1, corr);
        FLOAT  exclRangeBegin,exclRangeEnd;

        if((isinfinite(corr_)[0]) || (abs(corr + 1) < EPS))
        {

            break;
        }

        result[0][i - 1] = MIN(motIdx, profileIndex[motIdx - 1]);
        result[1][i - 1] = max(motIdx, profileIndex[motIdx - 1]);

        auto corrProfile = crosscov_(motIdx);
        corrProfile = min_nv_Includenan(1, elementWiseMultiply(
                extr_vfromv(corrProfile, 1, timeSeries.size() - subseqLen + 1), invnorm ));

        size_t matrixProfile_size = matrixProfile.size();
        for(int k = 0;k < matrixProfile_size;k++)
        {
            if(isnan(matrixProfile[k]))
            {
                corrProfile[k] = NAN;
            }
        }

        if(exclusionLen > 0)
        {
            for(int j = 1;j <= 2;j++)
            {
                exclRangeBegin = max(1, result[j - 1][i - 1] - exclusionLen + 1);

                exclRangeEnd = MIN(FLOAT (matrixProfile.size()), result[j - 1][i - 1] + exclusionLen - 1);

                for(int l = exclRangeBegin;l <= exclRangeEnd;l++)
                {
                    corrProfile[l-1] = NAN;
                }

            }
        }

        FLOAT  neighborCorr = corrProfile[0];
        int neighbor = 0;

        size_t corrProfile_size = corrProfile.size();

        for(int j = 3;j <= neighborCount + 2;j++)
        {

            neighborCorr = corrProfile[0];
            neighbor = 1;

            for(int o = 1;o <= corrProfile_size;o++)
            {

                if(corrProfile[o-1] > neighborCorr || (isnan(neighborCorr) && !isnan(corrProfile[o-1])) )
                {
                    neighborCorr = corrProfile[o-1];
                    neighbor = o;
                }

            }
            
            vector<FLOAT > neighborCorr_v(1, neighborCorr);

            if(isinfinite(neighborCorr_v)[0] || ((1 - neighborCorr) >= (radius * (1 - corr_[0]))))
            {

                break;
            }

            result[j - 1][i - 1] = neighbor;

            if(exclusionLen > 0)
            {
                exclRangeBegin = max(1.0, FLOAT (neighbor - exclusionLen + 1));
                exclRangeEnd = MIN(FLOAT (matrixProfile.size()), FLOAT (neighbor + exclusionLen - 1));

                for(int l = exclRangeBegin;l <= exclRangeEnd;l++)
                {
                    corrProfile[l-1] = NAN;
                }
            }
        }

        for(int k = 0;k < corrProfile_size;k++)
        {
            if(isnan(corrProfile[k]))
            {

                matrixProfile_[k] = NAN;
            }
        }

    }

    result = result;
}

