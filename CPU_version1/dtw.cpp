
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define INF 1E20

double dtw(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb,
           int &cb_prune, int low_index, int high_index)
{

    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;
    cb_prune = 0;
    
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        k = MAX(0,r-i);
        
        min_cost = INF;
        int index;

        for(j=max(0,i-r); j<=MIN(m-1,i+r); j++, k++)
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

            cost[k] = MIN( MIN( x, y) , z) + DIST(A[i], B[j]);
            double i_ = i;
            double j_ = j;

            DOUBLE d_ = cost[k];
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
                
                index = k;
            }
        }

        if(i == m/2)
        {
            if (i+r < m-1 && min_cost + cb[i+r+1] >= threshold_2)
            {   free(cost);
                free(cost_prev);
                cb_prune ++;
                return INFINITY;
            }
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return sqrt(final_dtw);
}