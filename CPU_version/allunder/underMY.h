
#ifndef DIST_MPX_V2O1_FAST6_UNDERMY_H
#define DIST_MPX_V2O1_FAST6_UNDERMY_H

#endif 
#include "vector"
using namespace std;

NEXT_RETURN next(const vector<DOUBLE>& X,const vector<DOUBLE>& f_x,const vector<DOUBLE>& bsf_X,DOUBLE start_pos_,DOUBLE end_pos_);
void dist_mpx_v2O1_selectedV6_(const vector<DOUBLE>& a, int minlag, int subseqLen, int warpmax, DOUBLE bsf,
                               double adjust_factor, RETURN_V4& result);
RETURN_V4 dist_mpx_v2O1_selectedV6_(const vector<DOUBLE>& ts,int minlag,int subseqlen,int warpmax,DOUBLE bsf,int adjust_factor);