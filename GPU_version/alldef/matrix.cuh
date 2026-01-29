
#ifndef DIST_MPX_V2O1_FAST6_MATRIX_H
#define DIST_MPX_V2O1_FAST6_MATRIX_H
#include "vector"
#include "typedefdouble.cuh"
#include <iostream>
#include <complex>
#include <cmath>
#include "allstruct.cuh"
#include "GPU_parameters.h"
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define DIST(x,y) ((x-y)*(x-y))
#define SHOW_DEBUG_STATS true
#define ALPHA_BIAS 0.1
#define FIND_POS_GLOBAL(LL,len_of_table_global) floor(LL*100 + len_of_table_global *0.5)
#define FIND_POS_LOCAL(LL,len_of_table_local) floor(LL*100 + len_of_table_local *0.5)
#define CHECK_BORDER(len_of_table_local,pos) MIN(MAX(0,pos),len_of_table-1)

using namespace std;
extern FLOAT cnt_up ,cnt_down;

std::vector<FLOAT > plusvector(const std::vector<FLOAT >&a, const std::vector<FLOAT >&b);
std::vector<FLOAT > movsum(const std::vector<FLOAT > &ts, int b);
std::vector<FLOAT > elementWiseMultiply(const std::vector<FLOAT >& vector1, const std::vector<FLOAT >& vector2);
std::vector<FLOAT > elementWiseMultiply_p(const FLOAT  vector1[], const FLOAT  vector2[], long long len);
FLOAT  norm_vector(const std::vector<FLOAT >& A, int b);
FLOAT  sum_vector(const std::vector<FLOAT >& vec);
std::vector<FLOAT > extr_vfromv(const std::vector<FLOAT >& arr, int st, int ed);
std::vector<FLOAT > addElementToFront(const std::vector<FLOAT >& arr, FLOAT  element);
vector<FLOAT > substractvector(const vector<FLOAT >&a, const vector<FLOAT >&b);
FLOAT  stddev(const std::vector<FLOAT >& data);
std::vector<FLOAT > elementWiseDivison_vv(const std::vector<FLOAT >& vector1, const std::vector<FLOAT >& vector2);

bool vectorisempty(const vector<FLOAT >& X);
std::vector<FLOAT > linspace(FLOAT  start, FLOAT  end, int num);
std::vector<FLOAT > mergeVectors(const std::vector<FLOAT >& v1, const std::vector<FLOAT >& v2, const std::vector<FLOAT >& v3);

std::vector<FLOAT > elementwiseMultiply_nv(FLOAT  scalar, const std::vector<FLOAT >& vector);
FLOAT  Mean(const std::vector<FLOAT >& values);
std::vector<bool> isinfinite(const std::vector<FLOAT >& vec);
std::vector<FLOAT > findNonZero(const std::vector<bool>& vec);
std::vector<bool> isNaN(const std::vector<FLOAT >& vec);
vector<FLOAT > convolve_valid(const std::vector<FLOAT >& vec, const std::vector<FLOAT >& core);
unsigned int nextpow2(unsigned int n);
std::vector<std::complex<FLOAT >> fft(const std::vector<FLOAT >& x, int n);
std::vector<std::complex<FLOAT >> conjugate(const std::vector<std::complex<FLOAT >>& input);
std::vector<std::complex<FLOAT >> elementWiseMultiply_complex(const std::vector<std::complex<FLOAT >>& vec1,
                                                             const std::vector<std::complex<FLOAT >>& vec2);
void fft1(std::vector<std::complex<FLOAT >>& a, bool inverse = false);
std::vector<FLOAT > ifft(const std::vector<std::complex<FLOAT >>& input, bool symmetric);
std::vector<FLOAT > min_nv_Includenan(FLOAT  x, const std::vector<FLOAT >& data);
std::vector<FLOAT > max_nv_Includenan(FLOAT  x, const std::vector<FLOAT >& data);

void mvmean(const std::vector<FLOAT > &a, int l, vector<FLOAT > &miu, vector<FLOAT > &si);

void lower_upper_lemire(const vector<FLOAT >& a, int n, int r, vector<FLOAT >& l, vector<FLOAT >& u);
void lower_upper_lemire(const FLOAT*  a, int n, int r, FLOAT*  l,FLOAT*  u);

FLOAT  elementWiseMultiply_sum(const std::vector<FLOAT >& vector1, const std::vector<FLOAT >& vector2);
__device__ FLOAT  elementWiseMultiply_p_sum(const FLOAT  vector1[], const FLOAT  vector2[], long long len);

std::vector<FLOAT > movsum_p(FLOAT *ts, int ts_size, int b);
vector<int> min_v_k(const std::vector<FLOAT >& v, int k);
#endif 
