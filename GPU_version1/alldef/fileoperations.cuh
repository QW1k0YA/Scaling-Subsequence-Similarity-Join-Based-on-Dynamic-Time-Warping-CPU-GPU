
#ifndef DIST_MPX_V2O1_FAST6_FILEOPERATIONS_H
#define DIST_MPX_V2O1_FAST6_FILEOPERATIONS_H
#include "matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allstruct.cuh"
#include "fstream"
vector<FLOAT > loadfile(string filename);

vector<FLOAT > loadfile_apart(string filename, int apart);
vector<FLOAT > loadfile_2(string filename1, string filename2);
vector<FLOAT > loadfile_(string filename);

void DTWMotifDiscoveryGUI(const vector<FLOAT> &TS, int subseqlen, int maxwarp, const char *file);
#endif 
