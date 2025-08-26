
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
using namespace std;
bool LB_KK(const vector<double> &q, const double *U, const double *L, long long seqlen, double threshold_2,
           vector<double> &cb, double special_shared_vector[], double &lbk) {

    DOUBLE dist = 0.0;
    double temp_u;
    double temp_l;
    
    dist=special_shared_vector[0]+special_shared_vector[1];

    double temp_dist = 0;
    for (size_t i = 0; i < 3 ; i++) {

        double d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        temp_dist += d;
    }

    for (size_t i = seqlen - 3; i < seqlen ; i++) {
        double d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        temp_dist += d;
    }

    cb[0]=special_shared_vector[0];
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=special_shared_vector[1];
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;
    dist = MAX(temp_dist,dist);

        for (size_t i = 3; i < seqlen - 3 ; i++)
        {
        double d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }

        dist += d;
        if(threshold_2 < dist)
        {
            return 0;
        }
        cb[i] = d;
    }
    lbk = dist;
    return 1;
}

bool LB_KK_FIRST(const vector<double> &q, const vector<double> &t, const vector<double> &U, const vector<double> &L,
                 long long seqlen, double threshold_2, vector<double> &cb) {

    DOUBLE dist = 0.0;
    double temp_u;
    double temp_l;

    int m = seqlen;
    double d = 0;
    double x0 = t[0] ;
    double y0 = t[(m - 1 )] ;
    double dleft=DIST(x0, q[0]);
    double dright=DIST(y0, q[m - 1]);

    double x1 = (t[( 1)] );
    d = min(DIST(x1, q[0]), DIST(x0, q[1]));
    d = min(d, DIST(x1, q[1]));
    dleft+=d;

    double y1 = (t[(m - 2 )]);
    d = min(DIST(y1, q[m - 1]), DIST(y0, q[m - 2]) );
    d = min(d, DIST(y1, q[m - 2]));
    dright+=d;

    if (dleft+dright>=threshold_2){
        return false;
    }
    else{
        double x2 = (t[(2)]);
        d = min(DIST(x0, q[2]), DIST(x1, q[2]));
        d = min(d, DIST(x2, q[2]));
        d = min(d, DIST(x2, q[1]));
        d = min(d, DIST(x2, q[0]));
        dleft += d;

        double y2 = (t[(m - 3 )]);
        d = min(DIST(y0, q[m - 3]), DIST(y1, q[m - 3]));
        d = min(d, DIST(y2, q[m - 3]));
        d = min(d, DIST(y2, q[m - 2]));
        d = min(d, DIST(y2, q[m - 1]));
        dright += d;

        if (dleft+dright> threshold_2){
            return false;
        }
    }

    cb[0]=dleft;
    cb[1]=0;
    cb[2]=0;
    cb[seqlen-3]=dright;
    cb[seqlen-2]=0;
    cb[seqlen-1]=0;

    dist = dleft + dright;
    for (size_t i = 3; i < seqlen - 3; i++) {
        d = 0;
        temp_u = U[i];
        temp_l = L[i];
        if (q[i] > temp_u) {
            d = (q[i] - temp_u) * (q[i] - temp_u);
        } else if (q[i] < temp_l) {
            d = (q[i] - temp_l) * (q[i] - temp_l);
        }
        dist += d;
        if(threshold_2 < dist)
        {
            return false;
        }
        cb[i] = d;
    }
    return true;
}