
#ifndef DIST_MPX_V2O1_FAST6_UNDERDTW_H
#define DIST_MPX_V2O1_FAST6_UNDERDTW_H

#include "vector"
#include "iostream"
#include "../alldef/typedefdouble.cuh"
#include "../alldef/allstruct.cuh"
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#define INITIAL_LEN_OF_TABLE 600

using namespace std;
constexpr FLOAT  INF{std::numeric_limits<FLOAT >::infinity()};
void dist_mpx_DIAGV3_fastF(const vector<FLOAT >& ts, int minlag, int subseqlen, int warpmax, FLOAT  bsf, const vector<vector<bool>>& lb_profile_previous, RETURN_FAST& result);
vector<vector<bool>> dist_mpx_v2O1_fast6(const vector<FLOAT > &ts, int minlag, int subseqlen, int warpmax, FLOAT  bsf);
void dist_mpx_v2O1_selectedV4_(const vector<FLOAT >& ts, int minlag, int subseqlen, int warpmax, FLOAT  bsf,
                               FLOAT adjust_factor, RETURN_V4& result);
void dist_mpx_v2O1_selectedV6_test_pro_max(const vector<FLOAT >& ts, int minlag,
                                           int subseqlen, int  warpmax, FLOAT  bsf, FLOAT  adjust_factor, const vector<vector<bool>>& lb_profile_previous, PROMAX_RETURN& result);
void dist_mpx_DIAGV3_fast(const vector<FLOAT >& ts, size_t minlag, size_t subseqlen,
                          size_t warpmax, FLOAT  bsf, const vector<vector<bool>>& lb_profile_previous, RETURN_FAST& result);

void dist_mpx_DIAGV3(const vector<FLOAT >& ts, int minlag, int subseqlen,
                     int warpmax, FLOAT  bsf, const vector<vector<bool>>& lb_profile_previous, RETURN_FAST& result);

FLOAT*  dtw_upd(const vector<FLOAT>& x,const vector<FLOAT> & y,long long warpmax,FLOAT best_so_far = INF);

RETURN_MY MY_dtw_mpGUIV2(const vector<FLOAT >& ts, FLOAT  current_purnning_rate, const vector<vector<bool>>& lb_profile_bool, const vector<FLOAT >& MPI, const vector<FLOAT >& MP, int subseqlen, int minlag, int warpmax, FLOAT  best_so_far);
FLOAT DTW_Motif_of_one_layer(const vector<FLOAT >& ts, long long subseqlen , int maxwarp,
                              long long & first_min, long long& sec_min, FLOAT  &p1, FLOAT  &p2, FLOAT  &p4, FLOAT  &p8,
                              FLOAT  &p16, FLOAT  &pdtw);
void
diag_fast(const vector<FLOAT > &ts, int subseqlen, int diag_ID, const vector<FLOAT > &UTS, const vector<FLOAT > &LTS,
          const vector<FLOAT > &mu, const vector<FLOAT > &sumU_sumL, const vector<FLOAT > &invsig,
          const vector<FLOAT > &norm_U_plus_norm_L_trans, const vector<FLOAT > &del, vector<bool> &lb_vector,
          const vector<FLOAT > &dr_bwdU_plus_dr_bwdL, const vector<FLOAT > &dc_bwd,
          const vector<FLOAT > &dr_fwdU_plus_dr_fwdL, const vector<FLOAT > &dc_fwd, FLOAT &cnt) ;
void
diag_fast(const FLOAT *ts, int subseqlen, int diag_ID, const FLOAT *UTS,const FLOAT *LTS,
          const FLOAT *mu,  const FLOAT *sumU_sumL, const FLOAT *invsig,
          const FLOAT *norm_U_plus_norm_L_trans,  const FLOAT *del,bool *lb_vector,
          const FLOAT *dr_bwdU_plus_dr_bwdL, const FLOAT *dc_bwd,
          const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_fwd, FLOAT &cnt,int EPOCH);
void
diag_mask_global(const vector<FLOAT > &ts, int subseqlen, int diagID, vector<bool> &lb_vector, const vector<FLOAT > &mu,
                 const vector<FLOAT > &UTS, const vector<FLOAT > &LTS, const vector<FLOAT > &MASK,
                 const vector<FLOAT > &TS2, const vector<FLOAT > &sumMASK, const vector<FLOAT > &invsig,
                 const vector<FLOAT > &sumU_sumL, const vector<FLOAT > &dr_bwdU_plus_dr_bwdL,
                 const vector<FLOAT > &dr_fwdU_plus_dr_fwdL, const vector<FLOAT > &dc_bwd, const vector<FLOAT > &dc_fwd,
                 const vector<FLOAT > &dr_bwdMASK, const vector<FLOAT > &dr_fwdMASK, const vector<FLOAT > &dc_bwdTS2,
                 const vector<FLOAT > &dc_fwdTS2, const vector<FLOAT > &DUL2, vector<FLOAT > &lb_vector_new, FLOAT  bsf,
                 const vector<vector<FLOAT >> &subs, FLOAT **special_shared_vector, FLOAT &cnt,
                 const vector<FLOAT > &norm_U_plus_norm_L_global, vector<FLOAT > &DUL);
void
diag_mask_global(const FLOAT *ts, int subseqlen, int diagID, bool *lb_vector, const FLOAT *mu, const FLOAT *UTS,
                 const FLOAT *LTS, const FLOAT *MASK, const FLOAT *TS2, const FLOAT *sumMASK,
                 const FLOAT *invsig, const FLOAT *sumU_sumL, const FLOAT *dr_bwdU_plus_dr_bwdL,
                 const FLOAT *dr_fwdU_plus_dr_fwdL, const FLOAT *dc_bwd, const FLOAT *dc_fwd,
                 const FLOAT *dr_bwdMASK, const FLOAT *dr_fwdMASK, const FLOAT *dc_bwdTS2,
                 const FLOAT *dc_fwdTS2, const FLOAT *DUL2, FLOAT *lb_vector_new, FLOAT  bsf, FLOAT **subs,
                 FLOAT **special_shared_vector, FLOAT &cnt, const FLOAT *DUL,
                 const FLOAT *norm_U_plus_norm_L_global, int start_pos, int end_pos);
void
LB_KIM(FLOAT threshold, int m, FLOAT **special_shared_vector, const vector<vector<FLOAT >> &subs, int temp_1,
       int diag, vector<bool> &lb_vector, FLOAT &cnt);
void LB_KIM_new(FLOAT threshold, int m, FLOAT **special_shared_vector, const vector<vector<FLOAT >> &subs, int temp_1,
                int diag, vector<bool> &lb_vector, FLOAT &cnt);

__device__ bool LB_KK(const FLOAT* q, const FLOAT *U, const FLOAT *L, long long seqlen, FLOAT threshold_2,
                      FLOAT *cb, FLOAT special_shared_vector[], FLOAT &lbk);
bool lbPetitjean_new(vector<FLOAT> &p, FLOAT  *q, vector<FLOAT> &up, vector<FLOAT> &lp, const FLOAT uq[], const FLOAT lq[],
                     FLOAT  *x, const FLOAT ux[], const FLOAT lx[], int w, int m, FLOAT threshold_2,
                     FLOAT &lbk2, vector<FLOAT> &cb, FLOAT miu, FLOAT si)  ;
bool lbPetitjean_new (vector<FLOAT> &p, FLOAT *q, vector<FLOAT> &up, vector<FLOAT> &lp, const FLOAT uq[], const FLOAT lq[],
                      FLOAT  *x, const FLOAT ux[], const FLOAT lx[], int w, int m, FLOAT threshold_2,
                      FLOAT &lbk2, vector<FLOAT> &cb, FLOAT miu, FLOAT si);

FLOAT MON_dtw_host(
        const FLOAT*  lines,
        const FLOAT*  cols,
        const FLOAT*  cb,
        int l,
        int w,
        FLOAT bsf
);
FLOAT SS_Pruned_dtw(const vector<FLOAT> &A, const vector<FLOAT> &B, int m, int r, FLOAT threshold_2, vector<FLOAT> &cb);
FLOAT SS_Pruned_dtw_debug(const vector<FLOAT> &A, const vector<FLOAT> &B, int m, int r, FLOAT threshold_2, vector<FLOAT> &cb,bool flag);
FLOAT*  compute_dtw(const vector<FLOAT>& A, const vector<FLOAT>& B, FLOAT *cb, int m, int r, FLOAT best_so_far);
bool LB_KK_FIRST(const vector<FLOAT> &q, const vector<FLOAT> &t, const vector<FLOAT> &U, const vector<FLOAT> &L,
                 long long seqlen, FLOAT threshold_2, vector<FLOAT> &cb);
bool LB_KK_FIRST(const FLOAT* q, const FLOAT* t, const FLOAT* U, const FLOAT* L,
                 long long seqlen, FLOAT threshold_2, vector<FLOAT> &cb);
void compute_shared_data_origin(const vector<FLOAT >& ts, int subseqlen, int warpmax, const vector<FLOAT > &mu,
                                const vector<FLOAT >& sig, const vector<vector<FLOAT >>&subs, RETURN_COM& result);
void compute_shared_data_local(const vector<FLOAT > &ts, int subseqlen, const vector<vector<FLOAT >> &subs,
                               const vector<vector<FLOAT >> &UU, const vector<vector<FLOAT >> &LL, FLOAT  &real_min,
                               FLOAT  &real_max, vector<FLOAT> &pos_UU, vector<FLOAT> &pos_LL, int &len_of_cdf,
                               vector<vector<short>> &count_table_cdf, FLOAT MAX_REAL_VALUE);
void
compute_shared_data_local(const vector<FLOAT > &ts, int subseqlen, FLOAT **subs,
                          FLOAT **UU, FLOAT **LL, FLOAT  &real_min,
                          FLOAT  &real_max, vector<FLOAT> &pos_UU, vector<FLOAT> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, FLOAT MAX_REAL_VALUE);
void compute_shared_data_global(const vector<FLOAT > &ts, int subseqLen, const vector<vector<FLOAT >> &my_subs,
                                int len_of_table_global, vector<vector<FLOAT >> &count_table_global);
void compute_shared_data_global(const vector<FLOAT > &ts, int subseqLen, FLOAT* * &my_subs,
                                long long int len_of_table_global, vector<vector<FLOAT >> &count_table_global);
void
diag_mask_local_down(const vector<FLOAT > &ts, int minlag, int subseqlen, FLOAT  bsf, int diagID,
                     vector<bool> &lb_vector, const vector<vector<SHORT>> &count_table,
                     const vector<FLOAT > &pos_UU, const vector<FLOAT > &pos_LL, const vector<FLOAT > &mu,
                     const vector<FLOAT > &UTS_, const vector<FLOAT > &LTS_, FLOAT **special_shared_vector,
                     FLOAT &cnt, vector<FLOAT > &sumU2_sumL2, vector<FLOAT > &sumU_sumL, vector<FLOAT > &sumMASK,
                     vector<FLOAT > &norm_U_plus_norm_L, vector<FLOAT > &raw_DIFF_UL, vector<FLOAT > &raw_DIFF_UL2,
                     vector<FLOAT > &raw_DIFF_UL2_temp, vector<FLOAT > &DUL2_raw, vector<FLOAT > &DUL2,
                     const vector<FLOAT > &TS2, FLOAT alpha, const vector<FLOAT> &invsig,
                     const vector<FLOAT> &invsig_2, int diag_bias);
void
diag_mask_local_down(const FLOAT* ts, int minlag, int subseqlen, FLOAT  bsf, int diagID,
                     bool *lb_vector, const vector<vector<SHORT>> &count_table,
                     const FLOAT* pos_UU, const FLOAT* pos_LL, const FLOAT* mu,
                     const FLOAT* UTS_, const FLOAT* LTS_, FLOAT **special_shared_vector,
                     FLOAT &cnt, FLOAT *sumU2_sumL2, FLOAT *sumU_sumL, FLOAT *sumMASK,
                     FLOAT *norm_U_plus_norm_L, FLOAT *raw_DIFF_UL, FLOAT *raw_DIFF_UL2,
                     FLOAT *raw_DIFF_UL2_temp, FLOAT *DUL2_raw, FLOAT *DUL2,
                     const FLOAT* TS2, FLOAT alpha, const FLOAT* invsig,
                     const FLOAT* invsig_2, int diag_bias, size_t len );
void
diag_mask_local_up(const vector<FLOAT > &ts, int minlag, int subseqlen, FLOAT  bsf, int diagID, vector<bool> &lb_vector,
                   const vector<vector<SHORT>> &count_table, const vector<FLOAT > &pos_UU, const vector<FLOAT > &pos_LL,
                   const vector<FLOAT > &mu, const vector<FLOAT > &UTS_, const vector<FLOAT > &LTS_,
                   FLOAT **special_shared_vector, FLOAT &cnt, vector<FLOAT > &sumU2_sumL2, vector<FLOAT > &sumU_sumL,
                   vector<FLOAT > &sumMASK, vector<FLOAT > &norm_U_plus_norm_L, vector<FLOAT > &raw_DIFF_UL,
                   vector<FLOAT > &raw_DIFF_UL2, vector<FLOAT > &raw_DIFF_UL2_temp, vector<FLOAT > &DUL2_raw,
                   vector<FLOAT > &DUL2, const vector<FLOAT > &TS2, FLOAT alpha, const vector<FLOAT> &invsig,
                   const vector<FLOAT> &invsig_2, int diag_bias);
__device__ void
diag_mask_local_up(const FLOAT *ts, int minlag, int subseqlen, FLOAT  bsf, int diagID, bool *lb_vector,
                   const vector<vector<SHORT>> &count_table, const FLOAT *pos_UU, const FLOAT *pos_LL,
                   const FLOAT *mu, const FLOAT *UTS_, const FLOAT *LTS_,
                   FLOAT **special_shared_vector, FLOAT &cnt, FLOAT *sumU2_sumL2, FLOAT *sumU_sumL,
                   FLOAT *sumMASK, FLOAT *norm_U_plus_norm_L, FLOAT *raw_DIFF_UL,
                   FLOAT *raw_DIFF_UL2, FLOAT *raw_DIFF_UL2_temp, FLOAT *DUL2_raw,
                   FLOAT *DUL2, const FLOAT *TS2, FLOAT alpha, FLOAT *invsig,
                   FLOAT *invsig_2, int diag_bias, int len);
void
diag_mask_local(const vector<FLOAT > &ts, int minlag, int subseqlen, FLOAT  bsf, int diagID, vector<bool> &lb_vector,
                const vector<vector<SHORT>> &count_table, const vector<FLOAT > &pos_UU, const vector<FLOAT > &pos_LL,
                const vector<FLOAT > &mu, const vector<FLOAT > &sig, int len_of_table, const vector<FLOAT > &UTS_,
                const vector<FLOAT > &LTS_, const vector<vector<FLOAT >> &subs, FLOAT **special_shared_vector,
                vector<FLOAT > &lb_vector_new, FLOAT &cnt, vector<FLOAT > &sumU, vector<FLOAT > &sumL,
                vector<FLOAT > &sumU2, vector<FLOAT > &sumL2, vector<FLOAT > &sumU2_sumL2, vector<FLOAT > &sumU_sumL,
                vector<FLOAT > &sumMASK, vector<FLOAT > &norm_U_plus_norm_L, vector<FLOAT> &del,
                vector<FLOAT > &raw_DIFF_UL, vector<FLOAT > &raw_DIFF_UL2, vector<FLOAT > &raw_DIFF_UL2_temp,
                vector<FLOAT > &DUL, vector<FLOAT > &DUL2_raw, vector<FLOAT > &DUL2, const vector<FLOAT > &TS2,
                FLOAT alpha, const vector<FLOAT> &invsig, const vector<FLOAT> &invsig_2);

void
compute_shared_data_local(const vector<FLOAT > &ts, int subseqlen, const vector<vector<FLOAT >> &subs,
                          const vector<vector<FLOAT >> &UU, const vector<vector<FLOAT >> &LL, FLOAT  &real_min,
                          FLOAT  &real_max, vector<FLOAT> &pos_UU, vector<FLOAT> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, FLOAT MAX_REAL_VALUE);
void compute_shared_data_global(const vector<FLOAT > &ts, int subseqLen, const vector<vector<FLOAT >> &my_subs,
                                int len_of_table_global, vector<vector<FLOAT >> &count_table_global);
void compute_MASK_global(const vector<FLOAT > &a, int subseqLen, int len_of_table_global,
                         vector<vector<FLOAT >> &count_table_global, const vector<vector<FLOAT >> &subs_U,
                         const vector<vector<FLOAT >> &subs_L, vector<FLOAT > &MASK_global);
void compute_MASK_global(const vector<FLOAT > &a, int subseqLen, int len_of_table_global,
                         vector<vector<FLOAT >> &count_table_global, FLOAT **subs_U,
                         FLOAT **subs_L, vector<FLOAT > &MASK_global);
#endif 

