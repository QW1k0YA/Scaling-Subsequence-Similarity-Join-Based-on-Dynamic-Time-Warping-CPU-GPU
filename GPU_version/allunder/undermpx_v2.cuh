
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"

#ifndef DIST_MPX_V2O1_FAST6_UNDERMPX_V2_H
#define DIST_MPX_V2O1_FAST6_UNDERMPX_V2_H
vector<FLOAT > moving_mean(vector<FLOAT > a, int w);
std::vector<FLOAT > sortInd(const std::vector<FLOAT >& vec);
vector<FLOAT > findDiscords(const vector<FLOAT >& matrixProfile, int exclusionLen);
void findMotifs(const vector<FLOAT >& timeSeries, const vector<FLOAT >& mu, const vector<FLOAT >& invnorm,
                const vector<FLOAT >& matrixProfile, const vector<int> &profileIndex, int subseqLen, int exclusionLen, vector<vector<int>> &result);
#endif 
