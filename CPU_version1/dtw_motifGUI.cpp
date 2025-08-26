
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
#include "alldef/typedefdouble.h"
#include "allunder/underdtw.h"

using namespace  std;

RETURN_GUI dtw_motifGUI(const vector<DOUBLE>& a,int subseqLen,int maxwarp,const vector<DOUBLE>& mp_ed)
{
    int subcount = a.size() - subseqLen + 1;
    auto my_total_time = std::chrono::high_resolution_clock::now();
    printf("\n\n\nlength: %d, maxwarp: %d\n",subseqLen,maxwarp);

    int minspacing = subseqLen;

    RETURN_MIN min01 = min_v(mp_ed);
    DOUBLE best_so_far = min01.minval;
    long long first_min = min01.minidx;

    vector<DOUBLE> mp_ed_ = mp_ed;

    for(DOUBLE i = max(1.0,first_min*1.0-subseqLen+1);i <= first_min + subseqLen -1;i++)
    {
        mp_ed_[i-1] = INFINITY;
    }
    RETURN_MIN min02 = min_v(mp_ed_);
    long long sec_min = min02.minidx;

    vector<DOUBLE> aa = extr_vfromv(a, first_min, first_min + subseqLen - 1);
    vector<DOUBLE> bb = extr_vfromv(a, sec_min, sec_min + subseqLen - 1);

    vector<DOUBLE> aa_temp;
    DOUBLE meanaa = meanWithoutNaN(aa);
    DOUBLE stdaa = stddev(aa);
    for(auto value:aa)
    {
        aa_temp.push_back((value - meanaa)/stdaa);
    }
    aa = aa_temp;

    vector<DOUBLE> bb_temp;
    DOUBLE meanbb = meanWithoutNaN(bb);
    DOUBLE stdbb = stddev(bb);
    for(auto value:bb)
    {
        bb_temp.push_back((value - meanbb)/stdbb);
    }
    bb = bb_temp;

    DOUBLE *res = dtw_upd(aa,bb,maxwarp);
    if(best_so_far > res[0])
    {
        best_so_far = res[0];
    }

    printf("################Start calculating lower bound#################  \n");

    int debug_sum = 0;
    auto t1 = chrono::high_resolution_clock ::now();

    vector<vector<bool>> lb_profile = dist_mpx_v2O1_fast6(a, subseqLen, subseqLen, maxwarp, best_so_far);

    RETURN_V4 temp_v4;
    dist_mpx_v2O1_selectedV4_(a, subseqLen, subseqLen,maxwarp,best_so_far,0,temp_v4);

    vector<vector<bool>> lb_profile2 = temp_v4.lb_profile_dummy;
    vector<DOUBLE> my_MPI = temp_v4.MPI;
    vector<DOUBLE> my_MP = temp_v4.MP;
    DOUBLE Pratio = temp_v4.Pratio;

    lb_profile = matrixOr(lb_profile,lb_profile2);

    auto endtime1 = chrono::high_resolution_clock ::now();
    auto O1_lb_time = std::chrono::duration_cast<std::chrono::microseconds>(endtime1 - t1).count()/1000000.0;

    int total = lb_profile[0].size();
    vector<DOUBLE> by_column;

    DOUBLE the_sereis_fully_purned = 0;
    
    if (SHOW_DEBUG_STATS){
        size_t lb_profile_size = lb_profile.size();
        for(int i = 1;i <= lb_profile_size;i++)
        {
            by_column.push_back(sum_vector_bool(lb_profile[i-1]));
        }

        for(auto value:by_column)
        {
            if(value > total - 0.5)
            {
                the_sereis_fully_purned+= value;
            }
        }
    }

    DOUBLE SSS = sum_bool_Matrix(lb_profile);
    DOUBLE purned_pieces_dived_the_total=SSS/DOUBLE (lb_profile.size()*lb_profile[0].size());

    printf("over all purnning rate %f  in series :%5.3f  \n",purned_pieces_dived_the_total,the_sereis_fully_purned/(a.size() - subseqLen +1));
    printf("Lower bound calculation time  %fs  bsf:%5.3f  \n",O1_lb_time,best_so_far);
    printf("####################start calculate the specific distance############## \n");

    vector<vector<bool>> empty_lb_profile(a.size() - subseqLen +1,vector<bool>(a.size() - subseqLen +1,0));
    vector<vector<bool>> lb_profile3;

    PROMAX_RETURN PRO_temp;
    dist_mpx_v2O1_selectedV6_test_pro_max(a, subseqLen, subseqLen,maxwarp,best_so_far, 0 , lb_profile,PRO_temp);
    lb_profile3 = PRO_temp.lb_profile_dummy;
    my_MPI = PRO_temp.MPI;
    my_MP = PRO_temp.MP;
    Pratio = PRO_temp.Pratio;
    lb_profile = matrixOr(lb_profile,lb_profile3);
    SSS = sum_bool_Matrix(lb_profile);
    purned_pieces_dived_the_total = SSS/DOUBLE(lb_profile.size()*lb_profile[0].size());

    auto t11 = chrono::high_resolution_clock ::now();

    RETURN_MY GUI_temp = MY_dtw_mpGUIV2(a,purned_pieces_dived_the_total,lb_profile, my_MPI, my_MP, subseqLen, minspacing, maxwarp , best_so_far);

    DOUBLE my_best_dtw = GUI_temp.best_so_far;
    DOUBLE my_motiffirst = GUI_temp.motiffirst;
    DOUBLE my_motifsec = GUI_temp.motifsec;

    auto endtime2 = chrono::high_resolution_clock ::now();
    auto our_compute_DTW_time = chrono::duration_cast<chrono::seconds>(endtime2 - t11).count();
    printf("DTW distance calculation time is %fs   input of bsf:%5.3f  ,final distance :%f  \n",DOUBLE(our_compute_DTW_time),best_so_far,my_best_dtw);
    fflush(stdout);
    if(my_motiffirst == -1)
    {
        printf("all series are cutted by ED motif");
        fflush(stdout);
        my_motiffirst = first_min;
        my_motifsec = sec_min;
    }

    aa = extr_vfromv(a, my_motiffirst, my_motiffirst + subseqLen - 1);
    bb = extr_vfromv(a, my_motifsec, my_motifsec + subseqLen - 1);

    vector<DOUBLE> aa_temp_;
    aa_temp = aa_temp_;
    meanaa = meanWithoutNaN(aa);
    stdaa = stddev(aa);
    for(auto value:aa)
    {
        aa_temp.push_back((value - meanaa)/stdaa);
    }
    aa = aa_temp;

    vector<DOUBLE> bb_temp_;
    bb_temp = bb_temp_;
    meanbb = meanWithoutNaN(bb);
    stdbb = stddev(bb);
    for(auto value:bb)
    {
        bb_temp.push_back((value - meanbb)/stdbb);
    }
    bb = bb_temp;

    DOUBLE dist_dtw = dtw_upd(aa,bb,maxwarp,best_so_far)[0];
    printf("\ntrue distance %f\n",dist_dtw);
    fflush(stdout);

    auto endtime3 = chrono::high_resolution_clock ::now();
    auto my_final_time = DOUBLE(chrono::duration_cast<chrono::microseconds>(endtime3 - my_total_time).count());

    printf("total time :%5.3fs\n",my_final_time/1000000);
    fflush(stdout);

    RETURN_GUI result;
    result.sp_rates = 0;
    result.pr_rate = 0;
    result.prrate = 0;
    result.DTW_time= 0;
    result.sec_min = sec_min;
    result.first_min = first_min;
    result.best_so_far = best_so_far;

    return  result;
}

