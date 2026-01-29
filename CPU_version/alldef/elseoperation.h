
#ifndef DIST_MPX_V2O1_FAST6_ELSEOPERATION_H
#define DIST_MPX_V2O1_FAST6_ELSEOPERATION_H

#endif 
#include "vector"
#include "typedefdouble.h"

using namespace  std;

DOUBLE LB_Keogh(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, long long seqlen);
DOUBLE LB_Keogh_ea(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, long long seqlen, DOUBLE bsf);
DOUBLE lb_upd(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, DOUBLE bsf) ;
void sampling_part1(int n, int subcount, int subseqlen, vector<vector<DOUBLE>> &subs, vector<vector<DOUBLE>> &U,
                    vector<vector<DOUBLE>> &L, vector<vector<DOUBLE>> &subs2, vector<vector<DOUBLE>> &U2,
                    vector<vector<DOUBLE>> &L2);
double sampling_part2(int n, int subseqlen, vector<vector<DOUBLE>> &subs2, vector<vector<DOUBLE>> &U2,
                      vector<vector<DOUBLE>> &L2, int id, int j, vector<DOUBLE> &X2);
void normalise(vector<double>&A,vector<double> &B,double mean,double std);;