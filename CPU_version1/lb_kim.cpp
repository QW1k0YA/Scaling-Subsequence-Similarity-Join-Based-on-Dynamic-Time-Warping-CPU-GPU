
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"

using namespace std;
void
LB_KIM(double threshold, int m, double **special_shared_vector, const vector<vector<DOUBLE>> &subs, int temp_1,
       int diag, vector<bool> &lb_vector, double &cnt) {

    for (int col = 0; col < temp_1; col++) {

        int row = diag + col - 1;
        if(lb_vector[col])
        {
            
            continue;
        }
        const vector<DOUBLE> &t_ = subs[row];
        const vector<DOUBLE> &q = subs[col];

        double  d;
        double threshold2=threshold*threshold;

        double x0 = t_[0] ;
        double y0 = t_[(m - 1 )] ;

        const double dleft_orgin=DIST(x0, q[0]);
        const double dright_orgin=DIST(y0, q[m - 1]);
        double dleft = dleft_orgin;
        double dright = dright_orgin;

        double x1 = (t_[( 1)] );
        const double d_left_weak = min(DIST(x1, q[0]), DIST(x1, q[1]));
        d = min(d_left_weak, DIST(x0, q[1]));
        dleft+=d;

        double y1 = (t_[(m - 2 )]);
        const double d_right_weak = min(DIST(y1, q[m - 1]),  DIST(y1, q[m - 2]));
        d = min(d_right_weak,DIST(y0, q[m - 2]));
        dright+=d;

        if (dleft+dright>=threshold2){
            cnt++;
            lb_vector[col] = true;
            continue;
        }
        else{

            d = MIN(DIST(x1,q[0]) + DIST(t_[2],q[0]),d_left_weak + DIST(t_[2],q[1]));
            d = MIN(d,DIST(q[1],t_[1]) + DIST(q[2],t_[2]));
            d = MIN(d,dleft + DIST(q[2],t_[1]));
            d = MIN(d,DIST(q[1],t_[0]) + DIST(q[2],t_[0]));
            dleft = d + dleft_orgin;

            d = MIN(DIST(t_[m-2],q[m-1]) + DIST(t_[m-3],q[m-1]),d_right_weak + DIST(t_[m-3],q[m-2]));
            d = MIN(d,DIST(q[m-2],t_[m-2]) + DIST(q[m-3],t_[m-3]));
            d = MIN(d,dright + DIST(q[m-3],t_[m-2]));
            d = MIN(d,DIST(q[m-2],t_[m-1]) + DIST(q[m-3],t_[m-1]));
            dright = d + dright_orgin;

            if (dleft+dright >=threshold2){
                lb_vector[col] = true;
                cnt++;
                continue;
            }
            special_shared_vector[col][0]=dleft; special_shared_vector[col][1]=dright;
        }
    }
}

void LB_KIM_new(double threshold, int m, double **special_shared_vector, const vector<vector<DOUBLE>> &subs, int temp_1,
       int diag, vector<bool> &lb_vector, double &cnt) {

    for (int col = 1; col <= temp_1; col++) {
        int row = diag + col - 1;

        if (lb_vector[col - 1]) {
            continue;
        }
        const vector<DOUBLE> &t_ = subs[row - 1];
        const vector<DOUBLE> &q = subs[col - 1];

        double d;
        double threshold2 = threshold * threshold;

        double x0 = t_[0];
        double y0 = t_[(m - 1)];

        const double dleft_orgin = DIST(x0, q[0]);
        const double dright_orgin = DIST(y0, q[m - 1]);
        double dleft = dleft_orgin;
        double dright = dright_orgin;

        double x1 = (t_[(1)]);
        const double d_left_weak = min(DIST(x1, q[0]), DIST(x1, q[1]));
        d = min(d_left_weak, DIST(x0, q[1]));
        dleft += d;

        double y1 = (t_[(m - 2)]);
        const double d_right_weak = min(DIST(y1, q[m - 1]), DIST(y1, q[m - 2]));
        d = min(d_right_weak, DIST(y0, q[m - 2]));
        dright += d;

        if (dleft + dright  >= threshold2) {
            cnt++;
            lb_vector[col-1] = true;
            continue;
        } else {

            d = MIN(DIST(x1, q[0]) + DIST(t_[2], q[0]), d_left_weak + DIST(t_[2], q[1]));
            d = MIN(d, DIST(q[1], t_[1]) + DIST(q[2], t_[2]));
            d = MIN(d, dleft + DIST(q[2], t_[1]));
            d = MIN(d, DIST(q[1], t_[0]) + DIST(q[2], t_[0]));
            dleft = d + dleft_orgin;

            d = MIN(DIST(t_[m - 2], q[m - 1]) + DIST(t_[m - 3], q[m - 1]), d_right_weak + DIST(t_[m - 3], q[m - 2]));
            d = MIN(d, DIST(q[m - 2], t_[m - 2]) + DIST(q[m - 3], t_[m - 3]));
            d = MIN(d, dright + DIST(q[m - 3], t_[m - 2]));
            d = MIN(d, DIST(q[m - 2], t_[m - 1]) + DIST(q[m - 3], t_[m - 1]));
            dright = d + dright_orgin;

            if (dleft + dright>= threshold2) {
                lb_vector[col-1] = true;
                cnt++;
                continue;
            }
            special_shared_vector[col-1][0] = dleft;
            special_shared_vector[col-1][1] = dright;
        }
    }
}
