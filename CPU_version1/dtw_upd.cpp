#include <vector>
#include <cmath>
#include <iostream>
#include <limits>
#include "alldef/elseoperation.h"
#include "alldef/typedefdouble.h"
#include "alldef/matrix.h"
constexpr DOUBLE INF{std::numeric_limits<DOUBLE>::infinity()};

#include<vector>
#include<limits>
#include<algorithm>
#include<cmath>

#include <iostream>
using namespace std;

double* compute_dtw(const vector<double>& A, const vector<double>& B, double *cb, int m, int r, double best_so_far) {
    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i, j, k;
    double x, y, z, min_cost;
    int ea = 0;
    double* res;
    
    res = (double*)calloc(2, sizeof(double));
    cost = (double*) calloc(2 * r + 1, sizeof(double));
    cost_prev = (double*) calloc(2 * r + 1, sizeof(double));
    for (k = 0; k < 2 * r + 1; k++) {
        cost[k] = INF;
        cost_prev[k] = INF;
    }

    for (i = 0; i < m; i++) {
        k = max(0, r - i);
        min_cost = INF;

        for (j = max(0, i - r); j <= MIN(m - 1, i + r); j++, k++) {
            
            if ((i == 0) && (j == 0)) {
                double c = (A[0] - B[0]);
                cost[k] = c * c;
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0)) {
                y = INF;
            } else {
                y = cost[k - 1];
            }
            if ((i < 1) || (k > 2 * r - 1)) {
                x = INF;
            } else {
                x = cost_prev[k + 1];
            }
            if ((i < 1) || (j < 1)) {
                z = INF;
            } else {
                z = cost_prev[k];
            }

            double d = A[i] - B[j];
            cost[k] = MIN(MIN(x, y) , z) + d * d;

            if (cost[k] < min_cost) {
                min_cost = cost[k];
            }
        }

        if (double((i+r+1))/m<=0.5 && min_cost + cb[i + r + 1] >= best_so_far) {
            free(cost);
            free(cost_prev);
            
            ea = ea + 1;

            res[0] = min_cost + cb[i + r + 1];
            res[1] = ea;
            return res;
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    res[0] = final_dtw;
    res[1] = ea;
    
    return res;
}

double* dtw(const vector<double>& A,const vector<double>& B, double *cb, int m, int r, double bsf = INF)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;
    int ea = 0;
    double*res = (double*)calloc(2, sizeof(double));

    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j <= MIN(m - 1, i + r); j++, k++)
        {
            
            if ((i==0)&&(j==0))
            {
                cost[k]=DIST(A[0], B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            cost[k] = MIN(MIN(x, y) , z) + DIST(A[i], B[j]);

            if (cost[k] < min_cost)
            {   min_cost = cost[k];
            }
        }

        if (i+r < m-1 && min_cost + cb[i+r+1] >= bsf)
        {   free(cost);
            free(cost_prev);
            ea = ea + 1;

            res[0] = min_cost + cb[i + r + 1];
            res[1] = ea;
            return res;
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    double final_dtw = cost_prev[k];
    res[0] = final_dtw;
    res[1] = ea;
    
    return res;
}

double* dtw_upd(const vector<double>& x,const vector<double> & y,long long warpmax,double best_so_far = INF){

    long long seqlen = x.size();
    std::vector<double> cb(seqlen, 0);
    double max_threshold;
    if(abs(best_so_far - INF) > 0.1)
    {
        max_threshold =  best_so_far * best_so_far;
    }
    else
    {
        max_threshold = INF;
    }

    double *res = dtw(x, y, cb.data(), seqlen, warpmax, max_threshold);
    double dist = res[0];
    res[0] = sqrt(dist);
    return res;

}
