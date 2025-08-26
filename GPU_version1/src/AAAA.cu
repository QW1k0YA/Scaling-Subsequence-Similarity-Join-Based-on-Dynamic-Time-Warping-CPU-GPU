
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

int main(int argc,char* argv[])
{

    string filename1 = argv[1];
    const char* filename2 = argv[2];
    int m = atoi(argv[3]);
    float w = atof(argv[4]);

    vector<FLOAT > TS1 = loadfile(filename1);
    vector<FLOAT > TS_1;
    size_t num = TS1.size()/2;
    for(int i = 0;i < num/2  ;i++)
    {
        TS_1.push_back(TS1[i]);
    }

    vector<FLOAT> TS = TS1;
    DTWMotifDiscoveryGUI(TS,m,MAX(round(m*w),1),filename2);

}

