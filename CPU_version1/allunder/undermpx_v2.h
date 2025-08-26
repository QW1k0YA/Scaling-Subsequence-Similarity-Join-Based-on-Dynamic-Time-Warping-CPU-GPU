
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.h"

#ifndef DIST_MPX_V2O1_FAST6_UNDERMPX_V2_H
#define DIST_MPX_V2O1_FAST6_UNDERMPX_V2_H
vector<DOUBLE> moving_mean(vector<DOUBLE> a,int w);
std::vector<DOUBLE> sortInd(const std::vector<DOUBLE>& vec);
vector<DOUBLE> findDiscords(const vector<DOUBLE>& matrixProfile,int exclusionLen);
void findMotifs(const vector<DOUBLE>& timeSeries, const vector<DOUBLE>& mu, const vector<DOUBLE>& invnorm,
                const vector<DOUBLE>& matrixProfile, const vector<int> &profileIndex, int subseqLen, int exclusionLen, RETURN_FImo& result);
#endif 
