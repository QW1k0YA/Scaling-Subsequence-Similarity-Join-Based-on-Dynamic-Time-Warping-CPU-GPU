
#ifndef DIST_MPX_V2O1_FAST6_UNDERDTW_H
#define DIST_MPX_V2O1_FAST6_UNDERDTW_H

#include "vector"
#include "iostream"
#include "../alldef/typedefdouble.h"
#include "../alldef/allstruct.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#define INITIAL_LEN_OF_TABLE 600

using namespace std;
constexpr DOUBLE INF{std::numeric_limits<DOUBLE>::infinity()};
void dist_mpx_DIAGV3_fastF(const vector<DOUBLE>& ts,int minlag,int subseqlen,int warpmax,DOUBLE bsf,const vector<vector<bool>>& lb_profile_previous,RETURN_FAST& result);
vector<vector<bool>> dist_mpx_v2O1_fast6(const vector<DOUBLE> &ts, int minlag, int subseqlen, int warpmax, DOUBLE bsf);
void dist_mpx_v2O1_selectedV4_(const vector<DOUBLE>& ts, int minlag, int subseqlen, int warpmax, DOUBLE bsf,
                               double adjust_factor, RETURN_V4& result);
void dist_mpx_v2O1_selectedV6_test_pro_max(const vector<DOUBLE>& ts, int minlag,
                                           int subseqlen, int  warpmax, DOUBLE bsf,DOUBLE adjust_factor,const vector<vector<bool>>& lb_profile_previous,PROMAX_RETURN& result);
void dist_mpx_DIAGV3_fast(const vector<DOUBLE>& ts,size_t minlag,size_t subseqlen,
                          size_t warpmax,DOUBLE bsf,const vector<vector<bool>>& lb_profile_previous,    RETURN_FAST& result);

void dist_mpx_DIAGV3(const vector<DOUBLE>& ts,int minlag,int subseqlen,
                     int warpmax,DOUBLE bsf,const vector<vector<bool>>& lb_profile_previous,RETURN_FAST& result);

double* dtw_upd(const vector<double>& x,const vector<double> & y,long long warpmax,double best_so_far = INF);

RETURN_MY MY_dtw_mpGUIV2(const vector<DOUBLE>& ts,DOUBLE current_purnning_rate,const vector<vector<bool>>& lb_profile_bool,const vector<DOUBLE>& MPI,const vector<DOUBLE>& MP,int subseqlen,int minlag,int warpmax,DOUBLE best_so_far);
double DTW_Motif_of_one_layer(const vector<DOUBLE>& ts, long long subseqlen , int maxwarp,
                              long long & first_min, long long& sec_min, DOUBLE &p1, DOUBLE &p2, DOUBLE &p4, DOUBLE &p8,
                              DOUBLE &p16, DOUBLE &pdtw);
void
diag_fast(const vector<DOUBLE> &ts, int subseqlen, int diag_ID, const vector<DOUBLE> &UTS, const vector<DOUBLE> &LTS,
          const vector<DOUBLE> &mu, const vector<DOUBLE> &sumU_sumL, const vector<DOUBLE> &invsig,
          const vector<DOUBLE> &norm_U_plus_norm_L_trans, const vector<DOUBLE> &del, vector<bool> &lb_vector,
          const vector<DOUBLE> &dr_bwdU_plus_dr_bwdL, const vector<DOUBLE> &dc_bwd,
          const vector<DOUBLE> &dr_fwdU_plus_dr_fwdL, const vector<DOUBLE> &dc_fwd, double &cnt,
          const vector<DOUBLE> &DUL, const vector<DOUBLE> &DUL2, double bsf);
void
diag_mask_global(const vector<DOUBLE> &ts, int subseqlen, int diagID, vector<bool> &lb_vector, const vector<DOUBLE> &mu,
                 const vector<DOUBLE> &UTS, const vector<DOUBLE> &LTS, const vector<DOUBLE> &MASK,
                 const vector<DOUBLE> &TS2, const vector<DOUBLE> &sumMASK, const vector<DOUBLE> &invsig,
                 const vector<DOUBLE> &sumU_sumL, const vector<DOUBLE> &dr_bwdU_plus_dr_bwdL,
                 const vector<DOUBLE> &dr_fwdU_plus_dr_fwdL, const vector<DOUBLE> &dc_bwd, const vector<DOUBLE> &dc_fwd,
                 const vector<DOUBLE> &dr_bwdMASK, const vector<DOUBLE> &dr_fwdMASK, const vector<DOUBLE> &dc_bwdTS2,
                 const vector<DOUBLE> &dc_fwdTS2, const vector<DOUBLE> &DUL2, vector<DOUBLE> &lb_vector_new, DOUBLE bsf,
                 const vector<vector<DOUBLE>> &subs, double **special_shared_vector, double &cnt,
                 const vector<DOUBLE> &norm_U_plus_norm_L_global, vector<DOUBLE> &DUL);
void
LB_KIM(double threshold, int m, double **special_shared_vector, const vector<vector<DOUBLE>> &subs, int temp_1,
       int diag, vector<bool> &lb_vector, double &cnt);
void LB_KIM_new(double threshold, int m, double **special_shared_vector, const vector<vector<DOUBLE>> &subs, int temp_1,
                int diag, vector<bool> &lb_vector, double &cnt);
bool LB_KK(const vector<double> &q, const double *U, const double *L, long long seqlen, double threshold_2,
           vector<double> &cb, double special_shared_vector[], double &lbk);
bool lbPetitjean_new(vector<double> &p, const vector<double>& q, vector<double> &up, vector<double> &lp, const double uq[], const double lq[],
                     const vector<double> &x, const double ux[], const double lx[], int w, int m, double threshold_2,
                     double &lbk2, vector<double> &cb, double miu, double si)  ;
double dtw(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb,
           int &cb_prune, int low_index, int high_index);
double MON_dtw(
        const vector<double>& lines,
        const vector<double>& cols,
        const vector<double>& cb,
        int l,
        int w,
        double bsf
);
double SS_Pruned_dtw(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb);
double SS_Pruned_dtw_debug(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb,bool flag);
double* compute_dtw(const vector<double>& A, const vector<double>& B, double *cb, int m, int r, double best_so_far);
bool LB_KK_FIRST(const vector<double> &q, const vector<double> &t, const vector<double> &U, const vector<double> &L,
                 long long seqlen, double threshold_2, vector<double> &cb);
void compute_shared_data_origin(const vector<DOUBLE>& ts,int subseqlen,int warpmax,const vector<DOUBLE> &mu,
                                const vector<DOUBLE>& sig,const vector<vector<DOUBLE>>&subs,RETURN_COM& result);
#endif 

