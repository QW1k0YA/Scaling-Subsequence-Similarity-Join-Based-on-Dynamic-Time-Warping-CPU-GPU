
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.h"
#ifndef DIST_MPX_V2O1_FAST6_UNDERFAST_H
#define DIST_MPX_V2O1_FAST6_UNDERFAST_H
using namespace std;

void
compute_shared_data_local(const vector<DOUBLE> &ts, int subseqlen, const vector<vector<DOUBLE>> &subs,
                          const vector<vector<DOUBLE>> &UU, const vector<vector<DOUBLE>> &LL, DOUBLE &real_min,
                          DOUBLE &real_max, vector<double> &pos_UU, vector<double> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, double MAX_REAL_VALUE);

#endif 

