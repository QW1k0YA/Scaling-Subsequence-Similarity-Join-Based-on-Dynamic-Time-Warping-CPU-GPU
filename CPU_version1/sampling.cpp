
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "alldef/allstruct.h"
#include "alldef//elseoperation.h"
using namespace std;

void sampling_part1(int n, int subcount, int subseqlen, vector<vector<DOUBLE>> &subs, vector<vector<DOUBLE>> &U,
                      vector<vector<DOUBLE>> &L, vector<vector<DOUBLE>> &subs2, vector<vector<DOUBLE>> &U2,
                      vector<vector<DOUBLE>> &L2)
{
    double p_my = 0;
    int real_len = floor(subseqlen / n);
    
    long long kml = subseqlen - n + 1;
    vector<DOUBLE> a2_temp(subseqlen);

    vector<DOUBLE> Ua2_temp(subseqlen);

    vector<DOUBLE> La2_temp(subseqlen);

    for(int id = 0;id < subcount;id++)
    {
        vector<DOUBLE>&Ua = U[id];
        vector<DOUBLE>&La = L[id];
        vector<DOUBLE>&a = subs[id];
        vector<DOUBLE> a2;
        vector<DOUBLE> Ua2;
        vector<DOUBLE> La2;
        mvmean_miu(a, subseqlen, n, a2_temp);
        mvmean_miu(Ua, subseqlen, n, Ua2_temp);
        mvmean_miu(La, subseqlen, n, La2_temp);

        for(int i = 0;i < kml;i++)
        {
            if(i%n == 0)
            {
                a2.push_back(a2_temp[i]);
                Ua2.push_back(Ua2_temp[i]);
                La2.push_back(La2_temp[i]);
            }
        }
        U2[id] = Ua2;
        L2[id] = La2;
        subs2[id] = a2;
    }
}

double sampling_part2(int n, int subseqlen, vector<vector<DOUBLE>> &subs2, vector<vector<DOUBLE>> &U2,
                      vector<vector<DOUBLE>> &L2, int id, int j, vector<DOUBLE> &X2)
{
    double p_my = 0;
    int real_len = floor(subseqlen / n);
    
    double distA = 0,distB = 0,dist2;

    vector<DOUBLE> &Ux2 = U2[id];
    vector<DOUBLE> &Lx2 = L2[id];
    vector<DOUBLE> &Q2 = subs2[j];
    vector<DOUBLE> &Uq2 = U2[j];
    vector<DOUBLE> &Lq2 = L2[j];

    for(int k = 0;k < real_len;k++) {
        if (Uq2[k] < X2[k]) {
            distA += (Uq2[k] - X2[k]) * (Uq2[k] - X2[k]);
        } else if (Lq2[k] > X2[k]) {
            distA += (Lq2[k] - X2[k]) * (Lq2[k] - X2[k]);
        }

        if (Ux2[k] < Q2[k]) {
            distB += (Ux2[k] - Q2[k]) * (Ux2[k] - Q2[k]);
        } else if (Lx2[k] > Q2[k]) {
            distB += (Lx2[k] - Q2[k]) * (Lx2[k] - Q2[k]);
        }
    }
    dist2 = n * max(distA, distB);
    return dist2;
}