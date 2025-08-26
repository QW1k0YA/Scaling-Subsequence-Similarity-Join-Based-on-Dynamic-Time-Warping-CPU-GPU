
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"

#include "../allunder/undermpx_v2.cuh"

using namespace  std;

RETURN_MPX mpx_v2(const vector<FLOAT >& timeSeries, int minlag, int subseqlen){

    int subcount =timeSeries.size() - subseqlen + 1;

    vector<FLOAT > timeSeries_ = timeSeries;

    vector<FLOAT > nanmap = findNonZero(isinfinite(movsum(timeSeries_, subseqlen - 1)));
    vector<FLOAT > nanIDX = findNonZero(isNaN(timeSeries));

    for(auto value : nanIDX)
    {
        timeSeries_[value - 1] = 0;
    }

    vector<FLOAT > mu = moving_mean(timeSeries_, subseqlen);
    vector<FLOAT > mus = moving_mean(timeSeries_, subseqlen - 1);

    vector<FLOAT > invnorm(subcount, 0.0);

    vector<FLOAT > timeSeries_1 = timeSeries_;

    for(int i = 1;i <= subcount;i++)
    {
        vector<FLOAT > temp_ts2;
        for(int j = i;j <= i + subseqlen -1;j++)
        {
            temp_ts2.push_back(timeSeries_[j-1] - mu[i-1]);
        }

        invnorm[i-1] = 1/ norm_vector(temp_ts2,2);
    }

    for(auto value:nanmap)
    {
        invnorm[value-1] = NAN;
    }

    int j = 0;
    for(auto value: isinfinite(invnorm))
    {
        if(value == 1)
        {
            invnorm[j] = NAN;
        }
        j++;
    }

    vector<FLOAT > dr_bwd = addElementToFront(substractvector(extr_vfromv(timeSeries_, 1, subcount - 1),
                                                             extr_vfromv(mu, 1, subcount - 1)), 0);
    vector<FLOAT > dc_bwd = addElementToFront(substractvector(extr_vfromv(timeSeries_, 1, subcount - 1),
                                                             extr_vfromv(mus, 2, subcount)), 0);
    vector<FLOAT > dr_fwd = substractvector(extr_vfromv(timeSeries_, subseqlen, timeSeries_.size()),
                                           extr_vfromv(mu, 1, subcount));
    vector<FLOAT > dc_fwd = substractvector(extr_vfromv(timeSeries_, subseqlen, timeSeries_.size()),
                                           extr_vfromv(mus, 1, subcount));

    vector<FLOAT > matrixProfile(subcount, -1);
    for(auto value:nanmap)
    {
        matrixProfile[value-1] = NAN;
    }

    vector<int> matrixProfileIdx(subcount,NAN);

    FLOAT  cov_,corr_;

    for(int diag = minlag + 1;diag <= subcount;diag ++)
    {
        vector<FLOAT > temp_ts,temp_ts1;
        for(int i = diag;i <= diag + subseqlen -1;i++)
        {
            temp_ts.push_back(timeSeries_[i-1]-mu[diag-1]);
        }
        for(int i = 1;i <= subseqlen;i++)
        {
            temp_ts1.push_back(timeSeries_[i-1]-mu[0]);
        }

        cov_ = elementWiseMultiply_sum(temp_ts,temp_ts1);

        for(int row = 1; row <= subcount - diag + 1; row++)
        {
            int col = diag +row -1;
            if(row > 1)
            {
                cov_ = cov_ - dr_bwd[row-1]*dc_bwd[col-1] + dr_fwd[row-1]*dc_fwd[col-1];
            }

            corr_ = cov_ * invnorm[row-1]*invnorm[col-1];

            if(corr_ > matrixProfile[row-1])
            {
                matrixProfile[row-1] = corr_;
                matrixProfileIdx[row-1] = col;
            }
            if(corr_ > matrixProfile[col-1])
            {
                matrixProfile[col-1] = corr_;
                matrixProfileIdx[col-1] = row;
            }
        }
    }

    vector<FLOAT > discordsIdx = findDiscords(matrixProfile, minlag);

    vector<vector<int>> motifIdx(12, vector<int>(3,NAN));;
    findMotifs(timeSeries_,mu,invnorm,matrixProfile,matrixProfileIdx,subseqlen,minlag,motifIdx);

    vector<FLOAT > ttt;
    vector<FLOAT > matrixProfmaile;

    for(auto v:matrixProfile)
    {
        ttt.push_back(2*subseqlen*(1-v));
    }
    for(auto value: max_nv_Includenan(0,ttt))
    {
        matrixProfmaile.push_back(sqrt(value));
    }

    if(nanIDX.size() > 0.1)
    {
        for(auto value:nanIDX)
        {
            timeSeries_[value-1] = NAN;
        }
    }

    RETURN_MPX result;
    result.motifIdxs = motifIdx;
    result.matrixProfile = matrixProfmaile;
    result.discordIdx = discordsIdx;
    result.matrixProfileIdx = matrixProfileIdx;

    return result;

}

