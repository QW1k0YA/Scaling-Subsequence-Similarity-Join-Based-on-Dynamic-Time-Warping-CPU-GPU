
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.h"
#ifndef DIST_MPX_V2O1_FAST6_UNDERPROMAX_H
#define DIST_MPX_V2O1_FAST6_UNDERPROMAX_H
using namespace std;

void
compute_shared_data_local(const vector<DOUBLE> &ts, int subseqlen, const vector<vector<DOUBLE>> &subs,
                          const vector<vector<DOUBLE>> &UU, const vector<vector<DOUBLE>> &LL, DOUBLE &real_min,
                          DOUBLE &real_max, vector<double> &pos_UU, vector<double> &pos_LL, int &len_of_cdf,
                          vector<vector<short>> &count_table_cdf, double MAX_REAL_VALUE);
void compute_shared_data_global(const vector<DOUBLE> &ts, int subseqLen,const vector<vector<DOUBLE>> &my_subs,
                                int len_of_table_global, vector<vector<DOUBLE>> &count_table_global);
void compute_MASK_global(const vector<DOUBLE> &a, int subseqLen,const vector<vector<DOUBLE>> &my_subs,
                         int len_of_table_global, vector<vector<DOUBLE>> &count_table_global,
                         const vector<vector<DOUBLE>> &subs_U,const vector<vector<DOUBLE>> &subs_L,
                         vector<DOUBLE>& MASK_global);
void dist_LB_in_one_Line_O1_V2(const vector<DOUBLE>& ts, int minlag,
                               int subseqlen,  DOUBLE bsf,int diagID,const vector<vector<bool>>& lb_profile_dummy,
                               const vector<vector<SHORT>> & count_table,const vector<DOUBLE>& pos_UU,const vector<DOUBLE>& pos_LL,
                               const vector<DOUBLE>& mu,const vector<DOUBLE>& sig,int len_of_table,const vector<DOUBLE>& UTS_,const vector<DOUBLE>& LTS_,LB_RETURN& result);

void dist_LB_in_one_Line_O1_V2_onlylb(const vector<DOUBLE>& ts, int minlag,
                                      int subseqlen,  DOUBLE bsf,int diagID,vector<vector<bool>>& lb_profile_dummy,
                                      const vector<vector<SHORT>> & count_table,const vector<DOUBLE>& pos_UU,const vector<DOUBLE>& pos_LL,
                                      const vector<DOUBLE>& mu,const vector<DOUBLE>& sig,int len_of_table,const vector<DOUBLE>& UTS_,const vector<DOUBLE>& LTS_);

void
diag_mask_local_down(const vector<DOUBLE> &ts, int minlag, int subseqlen, DOUBLE bsf, int diagID,
                     vector<bool> &lb_vector, const vector<vector<SHORT>> &count_table,
                     const vector<DOUBLE> &pos_UU, const vector<DOUBLE> &pos_LL, const vector<DOUBLE> &mu,
                     const vector<DOUBLE> &UTS_, const vector<DOUBLE> &LTS_, double **special_shared_vector,
                     double &cnt, vector<DOUBLE> &sumU2_sumL2, vector<DOUBLE> &sumU_sumL, vector<DOUBLE> &sumMASK,
                     vector<DOUBLE> &norm_U_plus_norm_L, vector<DOUBLE> &raw_DIFF_UL, vector<DOUBLE> &raw_DIFF_UL2,
                     vector<DOUBLE> &raw_DIFF_UL2_temp, vector<DOUBLE> &DUL2_raw, vector<DOUBLE> &DUL2,
                     const vector<DOUBLE> &TS2, double alpha, const vector<double> &invsig,
                     const vector<double> &invsig_2, int diag_bias);
void
diag_mask_local_up(const vector<DOUBLE> &ts, int minlag, int subseqlen, DOUBLE bsf, int diagID, vector<bool> &lb_vector,
                   const vector<vector<SHORT>> &count_table, const vector<DOUBLE> &pos_UU, const vector<DOUBLE> &pos_LL,
                   const vector<DOUBLE> &mu, const vector<DOUBLE> &UTS_, const vector<DOUBLE> &LTS_,
                   double **special_shared_vector, double &cnt, vector<DOUBLE> &sumU2_sumL2, vector<DOUBLE> &sumU_sumL,
                   vector<DOUBLE> &sumMASK, vector<DOUBLE> &norm_U_plus_norm_L, vector<DOUBLE> &raw_DIFF_UL,
                   vector<DOUBLE> &raw_DIFF_UL2, vector<DOUBLE> &raw_DIFF_UL2_temp, vector<DOUBLE> &DUL2_raw,
                   vector<DOUBLE> &DUL2, const vector<DOUBLE> &TS2, double alpha, const vector<double> &invsig,
                   const vector<double> &invsig_2, int diag_bias);
void
diag_mask_local(const vector<DOUBLE> &ts, int minlag, int subseqlen, DOUBLE bsf, int diagID, vector<bool> &lb_vector,
                const vector<vector<SHORT>> &count_table, const vector<DOUBLE> &pos_UU, const vector<DOUBLE> &pos_LL,
                const vector<DOUBLE> &mu, const vector<DOUBLE> &sig, int len_of_table, const vector<DOUBLE> &UTS_,
                const vector<DOUBLE> &LTS_, const vector<vector<DOUBLE>> &subs, double **special_shared_vector,
                vector<DOUBLE> &lb_vector_new, double &cnt, vector<DOUBLE> &sumU, vector<DOUBLE> &sumL,
                vector<DOUBLE> &sumU2, vector<DOUBLE> &sumL2, vector<DOUBLE> &sumU2_sumL2, vector<DOUBLE> &sumU_sumL,
                vector<DOUBLE> &sumMASK, vector<DOUBLE> &norm_U_plus_norm_L, vector<double> &del,
                vector<DOUBLE> &raw_DIFF_UL, vector<DOUBLE> &raw_DIFF_UL2, vector<DOUBLE> &raw_DIFF_UL2_temp,
                vector<DOUBLE> &DUL, vector<DOUBLE> &DUL2_raw, vector<DOUBLE> &DUL2, const vector<DOUBLE> &TS2,
                double alpha, const vector<double> &invsig, const vector<double> &invsig_2);

#endif 
