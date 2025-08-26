
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

int main(int argc,char* argv[])
{
    string filename1 = argv[1];
    const char* filename2 = argv[2];
   
    vector<DOUBLE> TS1 = loadfile(filename1);
    vector<DOUBLE> TS_1;
    size_t num = TS1.size()/2;
    for(int i = 0;i < num/2  ;i++)
    {
        TS_1.push_back(TS1[i]);
    }

    vector<int> L_num = {32,64,128,256,512,1024,2048};

    vector<double> r_num = {0.1};

    for (int i=0;i<r_num.size();i++) {
        for (int j=0;j<L_num.size();j++) {
            DTWMotifDiscoveryGUI(TS1, L_num[j], MAX(round(L_num[j] * r_num[i]), 1), filename2);
        }
    }

    vector<DOUBLE> TS = TS1;

}

