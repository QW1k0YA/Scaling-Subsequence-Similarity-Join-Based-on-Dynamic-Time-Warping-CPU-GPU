
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/elseoperation.h"
#include "allunder/underfastF.h"
#include "allunder/underdtw.h"
#include "alldef/allstruct.h"
#include "allunder/underMY.h"
using namespace std;

RETURN_MY MY_dtw_mpGUIV2(const vector<DOUBLE>& ts,DOUBLE current_purnning_rate,const vector<vector<bool>>& lb_profile_bool,const vector<DOUBLE>& MPI,const vector<DOUBLE>& MP,int subseqlen,int minlag,int warpmax,DOUBLE best_so_far)
{
    DOUBLE  old_best_so_far=best_so_far;
    int motiffirst=-1;
    int motifsec=-1;

    vector<DOUBLE> mu = movmean(ts,0,subseqlen-1,1);
    vector<DOUBLE> sig = movstd(ts,0,subseqlen-1,1);

    size_t subcount = ts.size() - subseqlen  + 1;
    vector<vector<DOUBLE>> subs(subseqlen,vector<DOUBLE>(subcount,0.0));

    for (int i = 1; i <= subcount; i++) {
        for (int j = i; j <= i + subseqlen - 1; j++) {
            subs[j - i][i - 1] = (ts[j - 1] - mu[i - 1]) / sig[i - 1];
        }
    }
    vector<vector<DOUBLE>> subs_trans = transposeMatrix_double(subs);

    vector<vector<DOUBLE>> U(subseqlen,vector<DOUBLE>(subcount,NAN));
    vector<vector<DOUBLE>> L(subseqlen,vector<DOUBLE>(subcount,NAN));

    if(warpmax > 0)
    {
        for(int i = 1;i <= subcount;i++)
        {
            vector<DOUBLE> U_temp = movmax(extractVector(subs,i-1,0),warpmax,warpmax);
            vector<DOUBLE> L_temp = movmin(extractVector(subs,i-1,0),warpmax,warpmax);
            size_t U_temp_size = U_temp.size();
            for(int j = 0;j < U_temp_size;j++)
            {
                U[j][i-1] = U_temp[j];
                L[j][i-1] = L_temp[j];
            }
        }
    }
    else
    {
        for(int i = 1;i <= subcount; i++)
        {
            vector<DOUBLE> U_temp = extractVector(subs,i-1,0);
            size_t U_temp_size = U_temp.size();
            for(int j = 0;j < U_temp_size;j++)
            {
                U[j][i-1] = U_temp[j];
                L[j][i-1] = U_temp[j];
            }
        }
    }

    vector<vector<double>> U_trans = transposeMatrix_double(U);
    vector<vector<double>> L_trans = transposeMatrix_double(L);

    DOUBLE lb_kim_time = 0,lb_keogh_time = 0,dtw_time = 0,lb_kim_iter = 0,lb_keogh_iter = 0,dtw_iter = 0,ea_cnt = 0;
    DOUBLE best_so_far2=best_so_far * best_so_far;

    DOUBLE cost_estimate_for_LB_KEOGH = (1-current_purnning_rate)* subseqlen;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    int cnt = 1;
    vector<DOUBLE> f_x;
    vector<DOUBLE> X;
    vector<DOUBLE> bsf_X;
    DOUBLE start_pos = -1.5;
    DOUBLE end_pos = 1.5;
    vector<vector<bool>> lb_profile_bool_ = lb_profile_bool;
    size_t lb_profile_bool_size = lb_profile_bool_.size();
    size_t lb_profile_bool0_size = lb_profile_bool_[0].size();
    vector<DOUBLE> MPI_ = MPI;
    vector<DOUBLE> MP_ = MP;

    while((cnt == 1) || (cost_estimate_for_LB_KEOGH >= 10)) {
        auto start_time = std::chrono::high_resolution_clock::now();

        for (int row = 1; row <= subcount; row++) {
            if (MPI_[row - 1] == INFINITY) {
                continue;
            }
            if (lb_profile_bool_[row - 1][static_cast<int>(MPI_[row - 1])-1]) {
                continue;
            }
            lb_profile_bool_[int(MIN(DOUBLE(row), MPI_[row - 1])) - 1][int(MAX(DOUBLE(row), MPI_[row - 1])) - 1] = true;

            const vector<DOUBLE>& a = subs_trans[row-1];
            const vector<DOUBLE>& b = subs_trans[int(MPI_[row - 1]) - 1];

            DOUBLE lb_Kim = MAX(abs(a[0] - b[0]), abs(a[a.size() - 1] - b[b.size() - 1]));
            if (lb_Kim >= best_so_far) {
                continue;
            }

            const vector<DOUBLE>& Ua = U_trans[row - 1];
            const vector<DOUBLE>& La = L_trans[row - 1];
            const vector<DOUBLE>& Ub = U_trans[int(MPI_[row - 1] - 1)];
            const vector<DOUBLE>& Lb = L_trans[int(MPI_[row - 1] - 1)];

            DOUBLE LB_Keogh = lb_upd(b, Ua, La, best_so_far2);

            if (LB_Keogh >= best_so_far2) {
                lb_keogh_iter = lb_keogh_iter + 1;
                continue;
            } else {
                LB_Keogh = lb_upd(a, Ub, Lb, best_so_far2);
                if (LB_Keogh >= best_so_far2) {
                    lb_keogh_iter = lb_keogh_iter + 1;
                    continue;
                }
            }

            auto resul2dtw = dtw_upd(a, b, warpmax,best_so_far);
            DOUBLE dist = resul2dtw[0];
            DOUBLE ea = resul2dtw[1];

            ea_cnt = ea + ea_cnt;
            dtw_iter = dtw_iter + 1;

            if (dist < best_so_far) {
                best_so_far = dist;
                best_so_far2 = best_so_far * best_so_far;
                motiffirst = row;
                motifsec = int(MPI_[row - 1]);

            }

        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        printf("round %d bsf update to %f, cost %5.3fs  \n", cnt, best_so_far, DOUBLE(duration.count())/1000000);

        DOUBLE del_1st = cost_estimate_for_LB_KEOGH * old_best_so_far / best_so_far;
        DOUBLE SSS;

        vector<vector<bool>> lb_profile2;

        DOUBLE Pratio;
        DOUBLE pos;

        if (cnt == 1) {
            
            if ((del_1st >= 10) || (old_best_so_far * 0.95 > best_so_far)) {
                auto ttt = std::chrono::high_resolution_clock::now();
                RETURN_V4 result1;
                dist_mpx_v2O1_selectedV4_(ts, subseqlen, subseqlen, warpmax, best_so_far, 0,result1);
                lb_profile2 = result1.lb_profile_dummy;
                MPI_ = result1.MPI;
                MP_ = result1.MP;

                Pratio = result1.Pratio;

                lb_profile_bool_ = matrixOr(lb_profile_bool_, lb_profile2);

                SSS = sum_bool_Matrix(lb_profile_bool_);
                current_purnning_rate = SSS / (lb_profile_bool_size * lb_profile_bool0_size);

                auto update_time = std::chrono::high_resolution_clock::now();
                auto duration_ = std::chrono::duration_cast<std::chrono::microseconds>(update_time - ttt);

                f_x.push_back(Pratio);
                X.push_back(0.0);
                bsf_X.push_back(best_so_far);
                printf("cnt:%d del= %5.2f, after update bsf , purnning rate %f ,cost %5.3fs  \n", cnt, del_1st,
                       current_purnning_rate, DOUBLE(duration_.count())/1000000);

            } else {
                printf("cnt:%d del= %5.2f,  Skip the first filter \n", cnt, del_1st);
                SSS = sum_bool_Matrix(lb_profile_bool_);
                Pratio = SSS / (lb_profile_bool_size * lb_profile_bool0_size);
                f_x.push_back(Pratio);
                X.push_back(0.0);
                bsf_X.push_back(best_so_far);
                cnt += 1;
                continue;
            }

        } else {
            if ((current_purnning_rate < 0.95) && cost_estimate_for_LB_KEOGH>10*(cnt-1)) {
                
                auto ttt = std::chrono::high_resolution_clock::now();
                NEXT_RETURN result1 = next(X, f_x, bsf_X, start_pos, end_pos);
                
                pos = result1.pos;
                start_pos = result1.start_pos;
                end_pos = result1.end_pos;

                fflush(stdout);
                RETURN_V4 result2;
                dist_mpx_v2O1_selectedV4_(ts, subseqlen, subseqlen, warpmax, best_so_far, pos,result2);
                fflush(stdout);
                lb_profile2 = result2.lb_profile_dummy;

                MPI_ = result2.MPI;
                MP_ = result2.MP;
                Pratio = result2.Pratio;
                lb_profile_bool_ = matrixOr(lb_profile_bool_, result2.lb_profile_dummy);
                SSS = sum_bool_Matrix(lb_profile_bool_);
                current_purnning_rate = SSS / (lb_profile_bool_size * lb_profile_bool0_size);

                auto update_time = std::chrono::high_resolution_clock::now();
                auto duration_ = std::chrono::duration_cast<std::chrono::microseconds>(update_time - ttt);

                f_x.push_back(Pratio);
                X.push_back(pos);
                bsf_X.push_back(best_so_far);
                fflush(stdout);
                printf("cnt:%d RK= %5.2f, after update bsf , purnning rate %f ,cost %5.3fs  \n", cnt,
                       cost_estimate_for_LB_KEOGH, current_purnning_rate, DOUBLE(duration_.count()/1000000));
                fflush(stdout);
            } else {
                RETURN_V4 result3;
                dist_mpx_v2O1_selectedV6_(ts, subseqlen, subseqlen, warpmax, best_so_far, -100,result3);
                lb_profile_bool_ = matrixOr(lb_profile_bool_, result3.lb_profile_dummy);
                printf("cnt:%d RK= %5.2f,  No further filtering \n", cnt, cost_estimate_for_LB_KEOGH);
                SSS = sum_bool_Matrix(lb_profile_bool_);
                current_purnning_rate = SSS / (lb_profile_bool_size * lb_profile_bool0_size);
                break;
            }
        }
        cost_estimate_for_LB_KEOGH = (1 - current_purnning_rate) * subseqlen;
        cnt = cnt + 1;

        if( end_pos - start_pos < 0.1)
        {
            RETURN_V4 result2;
            dist_mpx_v2O1_selectedV6_(ts,subseqlen,subseqlen,warpmax,best_so_far,-100,result2);
            lb_profile2 = result2.lb_profile_dummy;
            MPI_ = result2.MPI;
            MP_ = result2.MP;
            Pratio = result2.Pratio;
            lb_profile_bool_ = matrixOr(lb_profile_bool_, lb_profile2);
            printf("cnt:%d RK= %5.2f,  No further filtering \n", cnt, cost_estimate_for_LB_KEOGH);
            SSS = sum_bool_Matrix(lb_profile_bool_);
            current_purnning_rate = SSS / (lb_profile_bool_size * lb_profile_bool0_size);
            printf("The range is too small, we won¡¯t do it anymore \n ,%d,%f",cnt,cost_estimate_for_LB_KEOGH);
        }

    }

    auto end_time2 = std::chrono::high_resolution_clock::now();
    auto tmp_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time2 - start_time);
    printf("Two-factor authentication overhead %5.3fs  \n",DOUBLE(tmp_time.count())/1000000);

    auto start_time2 = std::chrono::high_resolution_clock::now();
    DOUBLE ea;

    for(int i = 1; i <=subcount; i++)
    {
        int id = i;

        const vector<DOUBLE>& a =subs_trans[id - 1];
        const vector<DOUBLE>& Ua = U_trans[id - 1];
        const vector<DOUBLE>& La = L_trans[id - 1];
        size_t a_size = a.size();
        size_t b_size = subs_trans[0].size();

        for(int idp = id + minlag;idp <= subcount;idp++)
        {
            if(lb_profile_bool_[id-1][idp-1])
            {
                continue;
            }

            const vector<DOUBLE>& b = subs_trans[idp - 1];
            DOUBLE lb_Kim = MAX(abs(a[0] - b[0]), abs(a[a_size - 1] - b[b_size - 1]));

            if(lb_Kim >= best_so_far)
            {
                lb_kim_iter = lb_kim_iter + 1;
                continue;
            }
            else
            {
                DOUBLE LB_Keogh = lb_upd(b,Ua,La,best_so_far2);

                if(LB_Keogh >= best_so_far2)
                {
                    lb_keogh_iter = lb_keogh_iter + 1;
                    continue;
                }
                else
                {

                    const vector<DOUBLE>&Ub = U_trans[idp - 1];
                    const vector<DOUBLE>&Lb = L_trans[idp - 1];
                    LB_Keogh=lb_upd(a,Ub,Lb,best_so_far2);
                    if(LB_Keogh >= best_so_far2)
                    {
                        lb_keogh_iter = lb_keogh_iter + 1;
                        continue;
                    }
                }
            }

            auto start_time3 = std::chrono::high_resolution_clock::now();

            double* temp = dtw_upd(a,b,warpmax,best_so_far);
            DOUBLE dist = temp[0];
            ea = temp[1];

            auto end_time3 = std::chrono::high_resolution_clock::now();
            auto tmp_time2 = std::chrono::duration_cast<std::chrono::seconds>(end_time3 - start_time3);
            ea_cnt += ea;
            dtw_iter+=1;
            dtw_time+= DOUBLE(tmp_time.count());

            if(dist < best_so_far)
            {
                best_so_far = dist;
                best_so_far2 = best_so_far*best_so_far;
                motiffirst = id;
                motifsec = idp;
            }
        }
    }

    auto end_time4 = std::chrono::high_resolution_clock::now();
    auto phase_time = std::chrono::duration_cast<std::chrono::seconds>(end_time4 - start_time2);
    ea_cnt += ea;
    DOUBLE pr_rate =  dtw_iter;
    printf("DTW is computed by %f times, costs: %5.3f \n",dtw_iter,dtw_time/1000000);

    RETURN_MY result;
    result.motifsec = motifsec;
    result.motiffirst = motiffirst;
    result.best_so_far = best_so_far;

    return result;

}

