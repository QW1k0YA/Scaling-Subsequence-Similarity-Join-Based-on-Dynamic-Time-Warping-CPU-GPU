
#ifndef DIST_MPX_V2O1_FAST6_ELSEOPERATION_H
#define DIST_MPX_V2O1_FAST6_ELSEOPERATION_H

#endif 
#include "vector"
#include "typedefdouble.cuh"

using namespace  std;

FLOAT  LB_Keogh(const std::vector<FLOAT >& x, const std::vector<FLOAT >& U, const std::vector<FLOAT >& L, long long seqlen);
FLOAT  LB_Keogh_ea(const std::vector<FLOAT >& x, const std::vector<FLOAT >& U, const std::vector<FLOAT >& L, long long seqlen, FLOAT  bsf);
FLOAT  lb_upd(const std::vector<FLOAT >& x, const std::vector<FLOAT >& U, const std::vector<FLOAT >& L, FLOAT  bsf) ;
void sampling_part1(int n, int subcount, int subseqlen, vector<vector<FLOAT >> &subs, vector<vector<FLOAT >> &U,
                    vector<vector<FLOAT >> &L, vector<vector<FLOAT >> &subs2, vector<vector<FLOAT >> &U2,
                    vector<vector<FLOAT >> &L2);
FLOAT sampling_part2(int n, int subseqlen, vector<vector<FLOAT >> &subs2, vector<vector<FLOAT >> &U2,
                      vector<vector<FLOAT >> &L2, int id, int j, vector<FLOAT > &X2);
void normalise(vector<FLOAT>&A,vector<FLOAT> &B,FLOAT mean,FLOAT std);;