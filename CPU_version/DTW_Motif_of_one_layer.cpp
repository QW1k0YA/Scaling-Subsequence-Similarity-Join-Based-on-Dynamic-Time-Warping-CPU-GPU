
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "time.h"
#include "alldef/elseoperation.h"
#include "allunder/underdtw.h"
using namespace std;

double DTW_Motif_of_one_layer(const vector<DOUBLE>& ts, long long subseqlen , int maxwarp,
                              long long & first_min, long long& sec_min, DOUBLE &p1, DOUBLE &p2, DOUBLE &p4, DOUBLE &p8,
                              DOUBLE &p16, DOUBLE &pdtw)

{
    
    DOUBLE my_total_time = clock();
    
    long long minlag = subseqlen;
    long long ts_size = ts.size();
    long long subcount = ts_size - subseqlen + 1;

    vector<DOUBLE> mu(subcount);
    vector<DOUBLE> sig(subcount);
    mvmean(ts, subseqlen, mu, sig);

    DOUBLE p64=0;
    DOUBLE p32=0;
    p16=0;
    p8=0;
    p4=0;
    p2=0;
    p1=0;
    DOUBLE p_my=0;
    DOUBLE dtw_cnt=0;

    vector<DOUBLE> aa = extr_vfromv(ts,first_min,first_min+subseqlen-1);
    vector<DOUBLE> bb = extr_vfromv(ts,sec_min,sec_min+subseqlen-1);

    DOUBLE meanaa = meanWithoutNaN(aa);
    DOUBLE stdaa = stddev(aa);
    for(int i = 0; i < subseqlen;i++)
    {
        aa[i] = (aa[i] - meanaa)/stdaa;
    }
    DOUBLE meanbb = meanWithoutNaN(bb);
    DOUBLE stdbb = stddev(bb);
    for(int i = 0; i < subseqlen;i++)
    {
        bb[i] = (bb[i] - meanbb)/stdbb;
    }

    double best_so_far = dtw_upd(aa,bb,maxwarp,INF) [0];
    double best_so_far2 = best_so_far*best_so_far;

    vector<vector<DOUBLE>> subs(subcount,vector<DOUBLE>(subseqlen,0));
    vector<vector<DOUBLE>> U(subcount,vector<DOUBLE>(subseqlen,NAN));
    vector<vector<DOUBLE>> L(subcount,vector<DOUBLE>(subseqlen,NAN));

    int real_len = floor(subseqlen / 2);
    vector<vector<DOUBLE>> U2(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> L2(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> subs2(subcount, vector<DOUBLE>(real_len, 0));

    real_len = floor(subseqlen / 4);
    vector<vector<DOUBLE>> U4(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> L4(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> subs4(subcount, vector<DOUBLE>(real_len, 0));

    real_len = floor(subseqlen / 8);
    vector<vector<DOUBLE>> U8(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> L8(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> subs8(subcount, vector<DOUBLE>(real_len, 0));

    real_len = floor(subseqlen / 16);
    vector<vector<DOUBLE>> U16(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> L16(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> subs16(subcount, vector<DOUBLE>(real_len, 0));

    real_len = floor(subseqlen / 32);
    vector<vector<DOUBLE>> U32(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> L32(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> subs32(subcount, vector<DOUBLE>(real_len, 0));

    real_len = floor(subseqlen / 64);
    vector<vector<DOUBLE>> U64(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> L64(subcount,vector<DOUBLE>(real_len,NAN));
    vector<vector<DOUBLE>> subs64(subcount, vector<DOUBLE>(real_len, 0));

    for(int i = 0;i < subcount;i++)
    {
        for(int j = 0;j < subseqlen;j++)
        {
            subs[i][j] = (ts[i + j] - mu[i])/sig[i];
        }
    }

    if(maxwarp > 0)
    {
       for(int i = 0;i < subcount;i++)
       {
           lower_upper_lemire(subs[i],subseqlen,maxwarp,L[i],U[i]);
       }
    }
    else
    {
        U = subs;
        L = subs;
    }

    sampling_part1(2, subcount, subseqlen, subs, U, L, subs2, U2, L2);
    sampling_part1(4, subcount, subseqlen, subs, U, L, subs4, U4, L4);
    sampling_part1(8, subcount, subseqlen, subs, U, L, subs8, U8, L8);
    sampling_part1(16, subcount, subseqlen, subs, U, L, subs16, U16, L16);
    sampling_part1(32, subcount, subseqlen, subs, U, L, subs32, U32, L32);
    sampling_part1(64, subcount, subseqlen, subs, U, L, subs64, U64, L64);

    double D = 10;
    double distA = 0,distB = 0;
    
    vector<DOUBLE> mean_t(ts_size - subseqlen + D + 1),std_t(ts_size - subseqlen + D + 1);

    mvmean(ts, subseqlen - D, mean_t, std_t);

    for(int id = 0;id < subcount;id++ )
    {
        vector<DOUBLE>& X2 = subs2[id];
        vector<DOUBLE>& X4 = subs4[id];
        vector<DOUBLE>& X8 = subs8[id];
        vector<DOUBLE>& X16 = subs16[id];
        vector<DOUBLE>& X32 = subs32[id];
        vector<DOUBLE>& X64 = subs64[id];

        for(int j = id + subseqlen;j < subcount;j ++)
        {
            vector<DOUBLE>& Q = subs[j];
            vector<DOUBLE>& Uq = U[j];
            vector<DOUBLE>& Lq = L[j];

            DOUBLE dist2 = sampling_part2(2, subseqlen, subs2, U2, L2, id, j, X2);
            if (dist2 > best_so_far2) {
                p2 += 1;
                continue;
            }
            DOUBLE dist4 = sampling_part2(4, subseqlen, subs4, U4, L4, id, j, X4);
            if (dist4 > best_so_far2) {
                p4 += 1;
                continue;
            }
            DOUBLE dist8 = sampling_part2(8, subseqlen, subs8, U8, L8, id, j, X8);
            if (dist8 > best_so_far2) {
                p8 += 1;
                continue;
            }
            DOUBLE dist16 = sampling_part2(16, subseqlen, subs16, U16, L16, id, j, X16);
            if (dist16 > best_so_far2) {
                p16 += 1;
                continue;
            }
            DOUBLE dist32 = sampling_part2(32, subseqlen, subs32, U32, L32, id, j, X32);
            if (dist32 > best_so_far2) {
                p32 += 1;
                continue;
            }
            DOUBLE dist64 = sampling_part2(64, subseqlen, subs64, U64, L64, id, j, X64);
            if (dist64 > best_so_far2) {
                p64 += 1;
                continue;
            }

            vector<DOUBLE> t_previous_id = extr_vfromv(ts,id + 1,id + subseqlen - D);
            normalise(t_previous_id,t_previous_id,mean_t[id],std_t[id]);

            vector<DOUBLE> t_previous_j = extr_vfromv(ts,j + 1,j + subseqlen - D);
            normalise(t_previous_j,t_previous_j,mean_t[j],std_t[j]);

            size_t previous_len = t_previous_j.size();

            vector<DOUBLE> U_previous_id = movmax(t_previous_id,maxwarp,maxwarp);
            vector<DOUBLE> L_previous_id = movmin(t_previous_id,maxwarp,maxwarp);
            vector<DOUBLE> U_previous_j = movmax(t_previous_j,maxwarp,maxwarp);
            vector<DOUBLE> L_previous_j = movmin(t_previous_j,maxwarp,maxwarp);

            vector<DOUBLE> dist_previous_id_v(previous_len);
            vector<DOUBLE> dist_previous_j_v(previous_len);

            double dist_previous_id = 0;
            double dist_previous_j = 0;

            for(int i = 0;i < previous_len;i++)
            {
                
                if(U_previous_j[i] < t_previous_id[i])
                {
                    dist_previous_id += (U_previous_j[i] - t_previous_id[i])*(U_previous_j[i] - t_previous_id[i]);
                }
                else if (L_previous_j[i] > t_previous_id[i])
                {
                    dist_previous_id += (L_previous_j[i] - t_previous_id[i])*(L_previous_j[i] - t_previous_id[i]);
                }

                if(U_previous_id[i] < t_previous_j[i])
                {
                    dist_previous_j += (U_previous_id[i] - t_previous_j[i])*(U_previous_id[i] - t_previous_j[i]);
                }
                else if (L_previous_id[i] > t_previous_j[i])
                {
                    dist_previous_j += (L_previous_id[i] - t_previous_j[i])*(L_previous_id[i] - t_previous_j[i]);
                }
            }

            vector<DOUBLE> changeU_id(previous_len);
            vector<DOUBLE> changeL_id(previous_len);
            vector<DOUBLE> changeU_j(previous_len);
            vector<DOUBLE> changeL_j(previous_len);
            vector<DOUBLE> change_ts_j(previous_len);
            vector<DOUBLE> change_ts_id(previous_len);

            for(int i = 0;i < previous_len;i++)
            {
                changeU_id[i] = U[id][i] - U_previous_id[i];
                changeL_id[i] = L[id][i] - L_previous_id[i];
                changeU_j[i] = Uq[i] - U_previous_j[i];
                changeL_j[i] = Lq[i] - L_previous_j[i];
                change_ts_id[i] = subs[id][i] - t_previous_id[i];
                change_ts_j[i] = Q[i] - t_previous_j[i];
            }

            vector<DOUBLE> mask_id(previous_len);
            vector<DOUBLE> mask_j(previous_len);

            for(int i = 0;i < previous_len;i++)
            {
                if((U_previous_j[i] < t_previous_id[i]) || (L_previous_j[i] > t_previous_id[i]))
                {
                    mask_id[i] = true;
                }
                if((U_previous_id[i] < t_previous_j[i]) || (L_previous_id[i] > t_previous_j[i]))
                {
                    mask_j[i] = true;
                }
            }

            double change1 = 0;
            vector<DOUBLE> changeUj2(changeU_j.size()),changeLj2(changeL_j.size());

            int len_change = change_ts_id.size();
            for(int i = 0;i < len_change; i++)
            {
                change1 += (change_ts_id[i]*change_ts_id[i]*mask_id[i]);
            }

            elementWiseMultiply(changeU_j, changeU_j, changeUj2);
            elementWiseMultiply(changeL_j, changeL_j, changeLj2);

            vector<DOUBLE> ttt(previous_len,0);
            for(int i =0;i < previous_len;i++)
            {
                if((abs(mask_id[i] - 1) < 0.001) && (changeUj2[i] >= changeLj2[i]) )
                {
                    ttt[i] += changeUj2[i];
                }
                else if((abs(mask_id[i] - 1) < 0.001) && (changeUj2[i] < changeLj2[i]) )
                {
                    ttt[i] += changeLj2[i];
                }
            }

            double tttt = sum_vector(ttt);
            change1 += 2*tttt;
            double dist_my_a= pow((MAX(sqrt(dist_previous_id) - sqrt(change1), 0.0) ), 2);

            DOUBLE change2 = 0;
            int len_change2 = change_ts_j.size();
            for(int i = 0;i < len_change2; i++)
            {
                change2 += (change_ts_j[i]*change_ts_j[i]);
            }

            vector<DOUBLE> changeUid2 = std::move(elementWiseMultiply(changeU_id,changeU_id));
            
            vector<DOUBLE> changeLid2 = std::move(elementWiseMultiply(changeL_id,changeL_id));

            ttt = std::move(plusvector(changeUid2,changeLid2));

            change2 += 2* sum_vector(ttt);
            double dist_my_b = pow((MAX(sqrt(dist_previous_j) - sqrt(change2), 0.0) ), 2);
            double dist_my = MAX(dist_my_a, dist_my_b);

            t_previous_id = extr_vfromv(ts,id + 1,id + subseqlen - D);
            double sigma1 = std_t[id];
            double mu1= mean_t[id];

            double sigma2 =sig[id];
            double mu2 = mu[id];

            double SUM_X = sum_vector(t_previous_id);

            double SUM_X2 = elementWiseMultiply_sum(t_previous_id,t_previous_id);

            DOUBLE SUM_DELT2= ( pow((sigma2-sigma1),2) * SUM_X2 + 2*(sigma1 * mu2 - sigma2 * mu1) * (sigma2-sigma1) * SUM_X +
                                (subseqlen-D)*pow((sigma1 *mu2- sigma2 * mu1),2) )/((sigma1 * sigma2)*(sigma1 * sigma2));

            vector<DOUBLE> normalized_t_previous_id(t_previous_id.size());
            normalise(t_previous_id,normalized_t_previous_id,mean_t[id],std_t[id]);

            subs_vector_p(subs[id].data(),normalized_t_previous_id.data(),change_ts_id.data(),previous_len);

            double del = SUM_DELT2 - elementWiseMultiply_sum(change_ts_id,change_ts_id);
            if(del >= 0.00001)
            {
                cerr << "wrong in samping";
            }

            if(dist_my > best_so_far2)
            {
                p_my += 1;
            }

            distA = 0;
            distB = 0;
            
            for(int i = 0;i < subseqlen;i++)
            {
                if(U[j][i] < subs[id][i])
                {
                    distA += (U[j][i] - subs[id][i])*(U[j][i] - subs[id][i]);
                }
                else if(L[j][i] > subs[id][i])
                {
                    distA += (L[j][i] - subs[id][i])*(L[j][i] - subs[id][i]);
                }
            }

            for(int i = 0;i < subseqlen;i++)
            {
                if(U[id][i] < subs[j][i])
                {
                    distB += (U[id][i] - subs[j][i])*(U[id][i] - subs[j][i]);
                }
                else if(L[id][i] > subs[j][i])
                {
                    distB += (L[id][i] - subs[j][i])*(L[id][i] - subs[j][i]);
                }
            }

            double dist = MAX(distA, distB);
            double ea = 0;
            if(dist_my > dist)
            {
                cerr << "error in one layer";
            }

            if(dist < dist2)
            {
                cerr << "error in samling2";
            }
            if(dist2 < dist4)
            {
                cerr << "error in samling2";
            }
            if(dist4 < dist8)
            {
                cerr << "error in samling2";
            }
            if(dist8 < dist16)
            {
                cerr << "error in samling2";
            }
            if(dist16 < dist32)
            {
                cerr << "error in samling2";
            }
            if(dist32 < dist64)
            {
                cerr << "error in samling2";
            }

            if(dist > best_so_far2)
            {
                p1 += 1;
                continue;
            }
            double*res = dtw_upd(subs[id],subs[j],maxwarp,best_so_far);
            dist = res[0];
            ea = res[1];

            dtw_cnt += 1;

            if(dist*dist < best_so_far2)
            {
                best_so_far = dist;
                best_so_far2 = best_so_far*best_so_far;
                first_min = id + 1;
                sec_min = j + 1;

            }

        }
    }

    aa = extr_vfromv(ts,first_min,first_min+subseqlen-1);
    bb = extr_vfromv(ts,sec_min,sec_min+subseqlen-1);

    meanaa = meanWithoutNaN(aa);
    stdaa = stddev(aa);
    for(int i = 0; i < subseqlen;i++)
    {
        aa[i] = (aa[i] - meanaa)/stdaa;
    }
    meanbb = meanWithoutNaN(bb);
    stdbb = stddev(bb);
    for(int i = 0; i < subseqlen;i++)
    {
        bb[i] = (bb[i] - meanbb)/stdbb;
    }

    double dist_dtw = dtw_upd(aa,bb,maxwarp,INF) [0];

    double total_time = clock();
    cout << "complete computing" << dist_dtw <<"," << (total_time - my_total_time)/CLOCKS_PER_SEC <<endl;

    double total_cnt = p16 + p8 + p4 + p2 + p1 + dtw_cnt;
    p_my /= total_cnt;
    p1 = p1/ total_cnt;
    p2 = p2/ total_cnt;
    p4 /= total_cnt;
    p8 /= total_cnt;
    p16 /= total_cnt;
    p32 /= total_cnt;
    p64 /= total_cnt;
    pdtw = dtw_cnt/total_cnt;
    printf("allocate  %5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f,%5.3f, %5.3f, %5.3f \n", p64,p32,p16,p8,p4,p2,p_my,p1,pdtw);

    return best_so_far;
}

