
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../allunder/underDTWGUI.cuh"
using namespace std;

/**
 * Main entry point for GPU-accelerated DTW motif discovery.
 * This function computes the Euclidean distance matrix profile first, then calls
 * the GPU-accelerated DTW motif discovery algorithm to find the most similar
 * subsequence pairs in the time series.
 *
 * @param TS: Input time series data
 * @param subseqlen: Length of subsequences to compare
 * @param maxwarp: Maximum warping window size for DTW
 * @param file: Path to output timing results
 */
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

    
    new_dtw_motifGUI_malloc_small_block(TS, subseqlen, maxwarp, mp, mp_index, file);
    

}