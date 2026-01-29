
#ifndef DIST_MPX_V2O1_FAST6_UNDERDTWGUI_H
#define DIST_MPX_V2O1_FAST6_UNDERDTWGUI_H
#include "vector"
#include "../alldef/allstruct.h"
#include "../alldef/typedefdouble.h"
using namespace std;

RETURN_MPX mpx_v2(const vector<DOUBLE>& timeSeries,int minlag,int subseqlen);
RETURN_GUI dtw_motifGUI(const vector<DOUBLE>& a,int subseqLen,int maxwarp,const vector<DOUBLE>& mp_ed);
void new_dtw_motifGUI(const vector<DOUBLE> &a, int subseqLen, int maxwarp, const vector<DOUBLE> &mp_ed,
                      const vector<int> &mp_ed_index, const char *output_time_file_path);
double SS_Pruned_dtw(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb);
#endif 
