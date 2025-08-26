
#ifndef DIST_MPX_V2O1_FAST6_FILEOPERATIONS_H
#define DIST_MPX_V2O1_FAST6_FILEOPERATIONS_H
#include "matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../allunder/underDIAGV3.h"
#include "allstruct.h"
#include "fstream"
vector<DOUBLE> loadfile(string filename);

vector<DOUBLE> loadfile_apart(string filename,int apart);
vector<DOUBLE> loadfile_2(string filename1,string filename2);
vector<DOUBLE> loadfile_(string filename);

void DTWMotifDiscoveryGUI(const vector<DOUBLE> &TS, int subseqlen, int maxwarp, const char *output_time_file_path);
#endif 
