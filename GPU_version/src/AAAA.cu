
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"

#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "fstream"
#include "../alldef/fileoperations.cuh"
#include "../alldef/matrix.cuh"
#include "../allunder/underdtw.cuh"

using namespace std;

/**
 * Main entry point for GPU-accelerated DTW motif discovery program.
 * This program takes a time series file and parameters, then performs
 * GPU-accelerated motif discovery using Dynamic Time Warping distance.
 *
 * Command line arguments:
 * - argv[1]: Input time series file path
 * - argv[2]: Output file path for timing results
 * - argv[3]: Subsequence length (m)
 * - argv[4]: Warping window ratio (w)
 *
 * @param argc: Number of command line arguments
 * @param argv: Array of command line argument strings
 * @return int: Exit status (0 for success)
 */
int main(int argc,char* argv[])
{

    string filename1 = argv[1];
    const char* filename2 = argv[2];
    int m = atoi(argv[3]);
    float w = atof(argv[4]);

    vector<FLOAT > TS1 = loadfile(filename1);
    vector<FLOAT> TS = TS1;
    DTWMotifDiscoveryGUI(TS,m,MAX(round(m*w),1),filename2);


}

