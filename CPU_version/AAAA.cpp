
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allunder/underDIAGV3.h"
#include "alldef/allstruct.h"
#include "fstream"
#include "alldef/fileoperations.h"
#include "alldef/matrix.h"
#include "allunder/underdtw.h"

using namespace std;

/**
 * Main entry point for DTW motif discovery with various parameter configurations.
 * This function loads a time series from a file and performs DTW-based motif discovery
 * with different combinations of subsequence lengths and warping window sizes.
 *
 * @param argc: Number of command-line arguments
 * @param argv: Array of command-line arguments
 *            - argv[1]: Input time series file path
 *            - argv[2]: Output time file path
 *            - argv[3]: Subsequence length (m)
 *            - argv[4]: Warping window ratio (w)
 * @return: 0 on successful completion
 */
int main(int argc,char* argv[])
{
    string filename1 = argv[1];
    const char* filename2 = argv[2];
    int m = atoi(argv[3]);
    DOUBLE w = atof(argv[4]);

    vector<DOUBLE > TS1 = loadfile(filename1);
    vector<DOUBLE> TS = TS1;
    DTWMotifDiscoveryGUI(TS,m,MAX(round(m*w),1),filename2);

}

