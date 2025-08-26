
#ifndef DIST_MPX_V2O1_FAST6_UNDERDTWGUI_H
#define DIST_MPX_V2O1_FAST6_UNDERDTWGUI_H
#include "vector"
#include "../alldef/allstruct.cuh"
#include "../alldef/typedefdouble.cuh"
using namespace std;

RETURN_MPX mpx_v2(const vector<FLOAT >& timeSeries, int minlag, int subseqlen);
RETURN_GUI dtw_motifGUI(const vector<FLOAT >& a, int subseqLen, int maxwarp, const vector<FLOAT >& mp_ed);
void new_dtw_motifGUI(const vector<FLOAT > &a, int subseqLen, int maxwarp, const vector<FLOAT > &mp_ed,
                      const vector<int> &mp_ed_index);
void new_dtw_motifGUI_malloc(const vector<FLOAT> &a, int subseqLen, int maxwarp, const vector<FLOAT> &mp_ed,
                             const vector<int> &mp_ed_index, const char *file);
FLOAT SS_Pruned_dtw(const vector<FLOAT> &A, const vector<FLOAT> &B, int m, int r, FLOAT threshold_2, vector<FLOAT> &cb);
#endif 
