
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define INF 1E20
__device__ FLOAT dtw(const FLOAT *A, const FLOAT *B, int m, int r, FLOAT threshold_2,   FLOAT *cost, FLOAT *cost_prev)
{

    FLOAT *cost_tmp;
    int i,j,k;
    FLOAT x,y,z,min_cost;

    for(k=0; k<2*r+1; k++)    cost[k]=INF;

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
            FLOAT i_ = i;
            FLOAT j_ = j;

            FLOAT  d_ = cost[k];
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
                
                index = k;
            }
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    FLOAT final_dtw = cost_prev[k];
    return sqrt(final_dtw);
}
FLOAT dtw(const vector<FLOAT> &A, const vector<FLOAT> &B, int m, int r, FLOAT threshold_2, vector<FLOAT> &cb)
{

    FLOAT *cost;
    FLOAT *cost_prev;
    FLOAT *cost_tmp;
    int i,j,k;
    FLOAT x,y,z,min_cost;

    cost = (FLOAT* )malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (FLOAT* )malloc(sizeof(double)*(2*r+1));
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
            FLOAT i_ = i;
            FLOAT j_ = j;

            FLOAT  d_ = cost[k];
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
                
                index = k;
            }
        }

        if (i+r < m-1 && min_cost + cb[i+r+1] >= threshold_2)
        {   free(cost);
            free(cost_prev);
            return sqrt(min_cost + cb[i+r+1]);
        }

        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    FLOAT final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return sqrt(final_dtw);
}

