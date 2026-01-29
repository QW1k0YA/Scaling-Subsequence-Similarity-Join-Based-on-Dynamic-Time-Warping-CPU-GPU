
#include "vector"
#ifndef DIST_MPX_V2O1_FAST6_ALLSTRUCT_H
#define DIST_MPX_V2O1_FAST6_ALLSTRUCT_H
using namespace std;
struct RETURN_COM{
    vector<vector<SHORT>> count_table;
    vector<DOUBLE> pos_UU;
    vector<DOUBLE> pos_LL;
    vector<DOUBLE> UTS;
    vector<DOUBLE> LTS;
    vector<bool> attached_mask;
    DOUBLE len_of_table;
};
struct RETURN_V4{
    vector<vector<bool>> lb_profile_dummy;
    vector<DOUBLE> MPI;
    vector<DOUBLE> MP;
    DOUBLE Pratio;
};
struct RETURN_FAST{
    vector<vector<bool>> lb_profile_dummy;
    vector<DOUBLE> MPI;
    vector<DOUBLE> MP;
    DOUBLE P__ratio;
};
struct PROMAX_RETURN
{
    vector<vector<bool>> lb_profile_dummy;
    vector<DOUBLE> MPI;
    vector<DOUBLE> MP;
    DOUBLE Pratio;
};
struct LB_RETURN{
    vector<vector<bool>> lb_profile_dummy;
    vector<DOUBLE> MPI;
    vector<DOUBLE> MP;
    DOUBLE Pratio;
};
struct LB_RETURN_diag{
    vector<bool> lb_profile_dummy;
    vector<DOUBLE> MPI;
    vector<DOUBLE> MP;
    DOUBLE Pratio;
};
struct NEXT_RETURN
{
    DOUBLE pos;
    DOUBLE start_pos;
    DOUBLE end_pos;
};
struct RETURN_MPX{
    vector<DOUBLE> matrixProfile;
    vector<int> matrixProfileIdx;
    vector<DOUBLE> discordIdx;
    vector<vector<int>> motifIdxs;
};
struct RETURN_FImo{
    vector<vector<int>> motifIdxs;
    vector<DOUBLE> matrixProfile;
};

struct RETURN_MIN{
    DOUBLE minval;
    DOUBLE minidx;
};
struct RETURN_GUI
{
    DOUBLE sp_rates;
    DOUBLE pr_rate;
    DOUBLE prrate;
    DOUBLE best_so_far;
    DOUBLE first_min;
    DOUBLE sec_min;
    DOUBLE DTW_time;
};
struct RETURN_MY{
    DOUBLE best_so_far;
    int motiffirst;
    int motifsec;
};

#endif 
