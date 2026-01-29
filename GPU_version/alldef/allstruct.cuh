
#include "vector"
#ifndef DIST_MPX_V2O1_FAST6_ALLSTRUCT_H
#define DIST_MPX_V2O1_FAST6_ALLSTRUCT_H
using namespace std;
struct RETURN_COM{
    vector<vector<SHORT>> count_table;
    vector<FLOAT > pos_UU;
    vector<FLOAT > pos_LL;
    vector<FLOAT > UTS;
    vector<FLOAT > LTS;
    vector<bool> attached_mask;
    FLOAT  len_of_table;
};
struct RETURN_V4{
    vector<vector<bool>> lb_profile_dummy;
    vector<FLOAT > MPI;
    vector<FLOAT > MP;
    FLOAT  Pratio;
};
struct RETURN_FAST{
    vector<vector<bool>> lb_profile_dummy;
    vector<FLOAT > MPI;
    vector<FLOAT > MP;
    FLOAT  P__ratio;
};
struct PROMAX_RETURN
{
    vector<vector<bool>> lb_profile_dummy;
    vector<FLOAT > MPI;
    vector<FLOAT > MP;
    FLOAT  Pratio;
};
struct LB_RETURN{
    vector<vector<bool>> lb_profile_dummy;
    vector<FLOAT > MPI;
    vector<FLOAT > MP;
    FLOAT  Pratio;
};
struct LB_RETURN_diag{
    vector<bool> lb_profile_dummy;
    vector<FLOAT > MPI;
    vector<FLOAT > MP;
    FLOAT  Pratio;
};
struct NEXT_RETURN
{
    FLOAT  pos;
    FLOAT  start_pos;
    FLOAT  end_pos;
};
struct RETURN_MPX{
    vector<FLOAT > matrixProfile;
    vector<int> matrixProfileIdx;
    vector<FLOAT > discordIdx;
    vector<vector<int>> motifIdxs;
};
struct RETURN_FImo{
    vector<vector<int>> motifIdxs;
    vector<FLOAT > matrixProfile;
};

struct RETURN_MIN{
    FLOAT  minval;
    FLOAT  minidx;
};
struct RETURN_GUI
{
    FLOAT  sp_rates;
    FLOAT  pr_rate;
    FLOAT  prrate;
    FLOAT  best_so_far;
    FLOAT  first_min;
    FLOAT  sec_min;
    FLOAT  DTW_time;
};
struct RETURN_MY{
    FLOAT  best_so_far;
    int motiffirst;
    int motifsec;
};

#endif 
