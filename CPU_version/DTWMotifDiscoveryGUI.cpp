
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
#include "allunder/underDTWGUI.h"
using namespace std;

/**
 * Main function for DTW-based motif discovery with Euclidean distance preprocessing.
 * This function first computes the Euclidean distance matrix profile and then
 * performs DTW-based motif discovery using the computed profile as initialization.
 *
 * @param TS: Input time series data
 * @param subseqlen: Length of subsequences to compare
 * @param maxwarp: Maximum warping window size for DTW
 * @param output_time_file_path: Path to output timing results
 */
void DTWMotifDiscoveryGUI(const vector<DOUBLE> &TS, int subseqlen, int maxwarp, const char *output_time_file_path)
{
    auto start_time1 = chrono::high_resolution_clock ::now();

    RETURN_MPX mpx_temp = mpx_v2(TS,subseqlen,subseqlen);
    vector<DOUBLE> mp = mpx_temp.matrixProfile;
    vector<int> mp_index = mpx_temp.matrixProfileIdx;

    auto end_time1 = chrono::high_resolution_clock ::now();

    auto MPlatency = std::chrono::duration_cast<std::chrono::microseconds>(end_time1 - start_time1).count()/1000000.0;
    printf("Time to finish ED MP computation: %f "
           "seconds\n",MPlatency);
    new_dtw_motifGUI(TS, subseqlen, maxwarp, mp, mp_index, output_time_file_path);
    cout << endl;
}