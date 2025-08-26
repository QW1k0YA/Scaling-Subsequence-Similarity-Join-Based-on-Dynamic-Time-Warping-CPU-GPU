
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../allunder/underDTWGUI.cuh"
using namespace std;

void DTWMotifDiscoveryGUI(const vector<FLOAT> &TS, int subseqlen, int maxwarp, const char *file)
{
    auto start_time1 = chrono::high_resolution_clock ::now();

    RETURN_MPX mpx_temp = mpx_v2(TS,subseqlen,subseqlen);
    vector<FLOAT > mp = mpx_temp.matrixProfile;
    vector<int> mp_index = mpx_temp.matrixProfileIdx;

    auto end_time1 = chrono::high_resolution_clock ::now();

    auto MPlatency = std::chrono::duration_cast<std::chrono::microseconds>(end_time1 - start_time1).count()/1000000.0;
    printf("Time to finish ED MP computation: %f "
           "seconds\n",MPlatency);

    new_dtw_motifGUI_malloc(TS, subseqlen, maxwarp, mp, mp_index, file);

}