
#ifndef DIST_MPX_V2O1_FAST6_MATRIX_H
#define DIST_MPX_V2O1_FAST6_MATRIX_H
#include "vector"
#include "typedefdouble.h"
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "allstruct.h"
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define DIST(x,y) ((x-y)*(x-y))
#define SHOW_DEBUG_STATS true

#define ALPHA_BIAS 0.1
#define FIND_POS_GLOBAL(LL,len_of_table_global) floor(LL*100 + len_of_table_global *0.5)
#define FIND_POS_LOCAL(LL,len_of_table_local) floor(LL*100 + len_of_table_local *0.5)
#define CHECK_BORDER(len_of_table_local,pos) MIN(MAX(0,pos),len_of_table-1)

using namespace std;
extern double cnt_up ,cnt_down;
std::vector<std::vector<DOUBLE>> transposeMatrix_double(const std::vector<std::vector<DOUBLE>>& matrix);
std::vector<std::vector<bool>> transposeMatrix_bool(const std::vector<std::vector<bool>>& matrix);
std::vector<std::vector<bool>> falseMatrix(std::size_t rows, std::size_t cols);
std::vector<std::vector<DOUBLE>> convertTo2DRowOrder(const std::vector<DOUBLE>& input);
std::vector<std::vector<DOUBLE>> convertTo2DColumnOrder(const std::vector<DOUBLE>& input);
std::vector<DOUBLE> movmean(const std::vector<DOUBLE>& ts, int a, int b,bool c);
std::vector<DOUBLE> movstd(const std::vector<DOUBLE>& ts, int a, int b,bool c);
std::vector<std::vector<DOUBLE>> elementWiseDivision(DOUBLE scalar, const std::vector<std::vector<DOUBLE>>& matrix);
std::vector<DOUBLE> extractVector(const std::vector<std::vector<DOUBLE>>& matrix, size_t index, bool isRow);
std::vector<DOUBLE> movmax(const std::vector<DOUBLE>& A, int kb, int kf);
std::vector<DOUBLE> movmin(const std::vector<DOUBLE>& A, int kb, int kf);
std::vector<DOUBLE> plusvector(const std::vector<DOUBLE>&a, const std::vector<DOUBLE>&b);
std::vector<DOUBLE> movsum(const std::vector<DOUBLE> &ts, int b);
std::vector<DOUBLE> elementWiseMultiply(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2);
std::vector<DOUBLE> elementWiseMultiply_from_pos(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2,int pos);
std::vector<DOUBLE> elementWiseMultiply_p(const DOUBLE vector1[], const DOUBLE vector2[], long long len);
DOUBLE norm_vector(const std::vector<DOUBLE>& A, int b);
DOUBLE sum_vector(const std::vector<DOUBLE>& vec);
std::vector<DOUBLE> extr_vfromv(const std::vector<DOUBLE>& arr, int st, int ed);
std::vector<std::vector<bool>> matrixOr(const std::vector<std::vector<bool>>& matrix1, const std::vector<std::vector<bool>>& matrix2);
int sum_bool_Matrix(const std::vector<std::vector<bool>>& matrix);
std::vector<DOUBLE> addElementToFront(const std::vector<DOUBLE>& arr, DOUBLE element);
vector<DOUBLE> substractvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b);
DOUBLE stddev(const std::vector<DOUBLE>& data);
DOUBLE maxinmatrix(const vector<vector<DOUBLE>>& ma);
DOUBLE mininmatrix(const vector<vector<DOUBLE>>& ma);
std::vector<DOUBLE> elementWiseDivison_vv(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2);
DOUBLE meanWithoutNaN(const std::vector<DOUBLE>& vec);
DOUBLE meanWithoutNaN_from_pos(const std::vector<DOUBLE>& vec,int pos);
bool vectorisempty(const vector<DOUBLE>& X);
std::vector<DOUBLE> linspace(DOUBLE start, DOUBLE end, int num);
std::vector<DOUBLE> mergeVectors(const std::vector<DOUBLE>& v1, const std::vector<DOUBLE>& v2, const std::vector<DOUBLE>& v3);
int sum_vector_bool(const std::vector<bool>& vec);
DOUBLE sum_matrix(const std::vector<std::vector<DOUBLE>>& matrix);
DOUBLE stdWithOmitnan(const std::vector<DOUBLE>& A);
DOUBLE dotProduct(const std::vector<DOUBLE>& A, const std::vector<DOUBLE>& B);
std::vector<DOUBLE> elementwiseMultiply_nv(DOUBLE scalar, const std::vector<DOUBLE>& vector);
DOUBLE Mean(const std::vector<DOUBLE>& values);
std::vector<bool> isinfinite(const std::vector<DOUBLE>& vec);
std::vector<DOUBLE> findNonZero(const std::vector<bool>& vec);
std::vector<bool> isNaN(const std::vector<DOUBLE>& vec);
vector<DOUBLE> convolve_valid(const std::vector<DOUBLE>& vec, const std::vector<DOUBLE>& core);
unsigned int nextpow2(unsigned int n);
std::vector<std::complex<DOUBLE>> fft(const std::vector<DOUBLE>& x, int n);
std::vector<std::complex<DOUBLE>> conjugate(const std::vector<std::complex<DOUBLE>>& input);
std::vector<std::complex<DOUBLE>> elementWiseMultiply_complex(const std::vector<std::complex<DOUBLE>>& vec1,
                                                            const std::vector<std::complex<DOUBLE>>& vec2);
void fft1(std::vector<std::complex<DOUBLE>>& a, bool inverse = false);
std::vector<DOUBLE> ifft(const std::vector<std::complex<DOUBLE>>& input, bool symmetric);
std::vector<DOUBLE> min_nv_Includenan(DOUBLE x, const std::vector<DOUBLE>& data);
std::vector<DOUBLE> max_nv_Includenan(DOUBLE x, const std::vector<DOUBLE>& data);
RETURN_MIN min_v(const std::vector<DOUBLE>& v);
void mvmean(const std::vector<DOUBLE> &a, int l, vector<DOUBLE> &miu, vector<DOUBLE> &si);
vector<DOUBLE> vectordividedbydouble(const vector<DOUBLE>&a, DOUBLE b);
DOUBLE* plusvector_p(const double *a, const double *b, long long len);
DOUBLE *plusvector_p_from_pos(const double *a, const double *b, long long len, long long pos);
void lower_upper_lemire(const vector<DOUBLE>& a, int n, int r, vector<DOUBLE>& l,vector<DOUBLE>& u);
void mvmean_miu(const std::vector<DOUBLE>& a,int len_a,int l,vector<DOUBLE>& miu);
void elementWiseMultiply(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE>& result);
void elementWiseMultiply_from_pos(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE>& result,int pos);
void subs_vector_p(const double *A,const double *B,double *C,long long len);
std::vector<DOUBLE> elementWiseMultiply_p_plus(const DOUBLE vector1[], const DOUBLE vector2[],const DOUBLE vector3[] ,long long len);
DOUBLE elementWiseMultiply_sum(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2);
DOUBLE elementWiseMultiply_p_plus_sum(const DOUBLE *vector1, const DOUBLE *vector2, const DOUBLE *vector3, long long len);
DOUBLE elementWiseMultiply_p_sum(const DOUBLE vector1[], const DOUBLE vector2[], long long len);
void elementWiseMultiply(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE>& result,int len);
void plusvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b, vector<DOUBLE>&result);
void elementWiseDivison_vv(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE> &result);
std::vector<DOUBLE> movsum_p(double *ts,int ts_size, int b);
void movsum(const std::vector<DOUBLE> &ts, int b,vector<DOUBLE> &result);
void movsum_p(double *ts,int ts_size, int b,vector<double> &result);
void movsum_p_from_pos(double *ts,int ts_size, int b,vector<double> &result,int pos);
void substractvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b,vector<DOUBLE> &result);
void substractvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b,vector<DOUBLE> &result,int len);
void substractvector_from_pos(const vector<DOUBLE>&a, const vector<DOUBLE>&b,vector<DOUBLE> &result,int pos);
vector<int> min_v_k(const std::vector<DOUBLE>& v,int k);
#endif 
