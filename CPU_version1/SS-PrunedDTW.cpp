
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
#define INF 1E20

double SS_Pruned_dtw(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb)
{

    int sc = 0;
    int ec = 0;
    int lp = 0;
    bool smaller_found;
    bool pruned_ec;
    int ec_next;
    int sc_last = 0;
    double ub;

    double *D;
    double *D_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    D = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    D[k]=INF;

    D_prev = (double*)malloc(sizeof(double)*(2*r+1));
    
    for(k=0; k<2*r+1; k++)    D_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        smaller_found = false;
        pruned_ec = false;
        ec_next = i;
        ub = (threshold_2 - cb[i+r+1]);

        k = MAX(0, MAX(r-i,sc-(i - r)));
        
        min_cost = INF;

        for(j=MAX(MAX(0,sc),i-r); j<=MIN(m-1,i+r); j++, k++)
        {

            if ((i==0)&&(j==0))
            {
                D[k]=DIST(A[0], B[0]);
                min_cost = D[k];
                if(D[k] <= ub)
                {
                    smaller_found = true;
                    ec_next = 1;
                }
                continue;
                
            }

            if ((j-1<0)||(k-1<0) ||  j - 1 < sc )     y = INF; 
            else                      y = D[k-1]; 
            if ((i-1<0)||(k+1>2*r) ||(j >= lp))   x = INF; 
            else                      x = D_prev[k+1]; 
            if ((i-1<0)||(j-1<0 ) || (j > lp) || j - 1 < sc )     z = INF; 
            else                      z = D_prev[k]; 

            D[k] = MIN(MIN(x, y) , z) + DIST(A[i], B[j]);
            double i_ = i;
            double j_ = j;

            DOUBLE d_ = D[k];
            if (D[k] < min_cost)
            {
                min_cost = D[k];
            }

            if(D[k] > ub)
            {
                if(j >= ec)
                {
                    lp = j;
                    pruned_ec = true;
                    
                    break;
                }
            }
            else
            {
                if(!smaller_found) {
                    
                    sc = j;
                    
                    smaller_found = true;
                }
                ec_next = j + 1;
                
            }
            
        }

        if (i+r < m-1 && min_cost + cb[i+r+1] >= threshold_2)
        {   free(D);
            free(D_prev);
            return sqrt(min_cost + cb[i+r+1]);
        }

        cost_tmp = D;
        D = D_prev;
        D_prev = cost_tmp;

        if(!pruned_ec)
        {
            lp = i + 1 + r;
        }
        ec = ec_next;
    }

    if(k == 0)
    {
        free(D_prev);
        free(D);
        
        return INF;
    }
    
    k--;
    if(pruned_ec)
    {
        D_prev[k] = INF;
        
    }
    
    double final_dtw = D_prev[k];
    free(D_prev);
    free(D);

    return sqrt(final_dtw);
}

double SS_Pruned_dtw_debug(const vector<double> &A, const vector<double> &B, int m, int r, double threshold_2, vector<double> &cb,bool flag)
{

    int sc = 0;
    int ec = 0;
    int lp = 0;
    bool smaller_found;
    bool pruned_ec;
    int ec_next;
    int sc_last = 0;
    double ub;

    double *D;
    double *D_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    D = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    D[k]=INF;

    D_prev = (double*)malloc(sizeof(double)*(2*r+1));
    
    for(k=0; k<2*r+1; k++)    D_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        smaller_found = false;
        pruned_ec = false;
        ec_next = i;
        ub = threshold_2 - cb[i+r+1];

        k = MAX(0, MAX(r-i,sc-(i - r)));
        
        min_cost = INF;

        for(j=MAX(MAX(0,sc),i-r); j<=MIN(m-1,i+r); j++, k++)
        {

            if ((i==0)&&(j==0))
            {
                D[k]=DIST(A[0], B[0]);
                min_cost = D[k];
                if(D[k] <= ub)
                {
                    smaller_found = true;
                    ec_next = 1;
                }
                continue;
                
            }

            if ((j-1<0)||(k-1<0) ||  j - 1 < sc )     y = INF; 
            else                      y = D[k-1]; 
            if ((i-1<0)||(k+1>2*r) ||(j >= lp))   x = INF; 
            else                      x = D_prev[k+1]; 
            if ((i-1<0)||(j-1<0 ) || (j > lp) || j - 1 < sc )     z = INF; 
            else                      z = D_prev[k]; 

            D[k] = MIN(MIN(x, y) , z) + DIST(A[i], B[j]);
            double i_ = i;
            double j_ = j;

            DOUBLE d_ = D[k];
            if(flag) cout << "(" << i << "," << j << ") " << d_ << "  ";
            if (D[k] < min_cost)
            {
                min_cost = D[k];
                
            }

            if(D[k] > ub)
            {
                if(j >= ec)
                {
                    lp = j;
                    pruned_ec = true;
                    
                    break;
                }
            }
            else
            {
                if(!smaller_found) {
                    sc_last = sc;
                    sc = j;
                    
                    smaller_found = true;
                }
                ec_next = j + 1;
                
            }
            
        }
        if(flag)
        {
            cout << endl;
        }

        if (i+r < m-1 && min_cost + cb[i+r+1] >= threshold_2)
        {   free(D);
            free(D_prev);

            return sqrt(min_cost + cb[i+r+1]);
        }

        cost_tmp = D;
        D = D_prev;
        D_prev = cost_tmp;

        if(!pruned_ec)
        {
            lp = i + 1 + r;
        }
        ec = ec_next;
    }

    if(k == 0)
    {
        free(D_prev);
        free(D);
        
        return INF;
    }
    
    k--;
    if(pruned_ec)
    {
        D_prev[k] = INF;
        cout << "pruned out  ub = "<< ub << endl;
        
    }
    
    double final_dtw = D_prev[k];
    free(D_prev);
    free(D);

    return sqrt(final_dtw);
}