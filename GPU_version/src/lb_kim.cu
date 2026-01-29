
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"

using namespace std;
void
LB_KIM(FLOAT threshold, int m, FLOAT **special_shared_vector, const vector<vector<FLOAT >> &subs, int temp_1,
       int diag, vector<bool> &lb_vector, FLOAT &cnt) {

    for (int col = 1; col <= temp_1; col++) {

        int row = diag + col - 1;

        if(lb_vector[col - 1])
        {
            
            continue;
        }
        const vector<FLOAT > &t_ = subs[row - 1];
        const vector<FLOAT > &q = subs[col - 1];

        FLOAT  d;
        FLOAT threshold2=threshold*threshold;

        FLOAT x0 = t_[0] ;
        FLOAT y0 = t_[(m - 1 )] ;

        const FLOAT dleft_orgin=DIST(x0, q[0]);
        const FLOAT dright_orgin=DIST(y0, q[m - 1]);
        FLOAT dleft = dleft_orgin;
        FLOAT dright = dright_orgin;

        FLOAT x1 = (t_[( 1)] );
        const FLOAT d_left_weak = min(DIST(x1, q[0]), DIST(x1, q[1]));
        d = min(d_left_weak, DIST(x0, q[1]));
        dleft+=d;

        FLOAT y1 = (t_[(m - 2 )]);
        const FLOAT d_right_weak = min(DIST(y1, q[m - 1]),  DIST(y1, q[m - 2]));
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

void LB_KIM_new(FLOAT threshold, int m, FLOAT **special_shared_vector, const vector<vector<FLOAT >> &subs, int temp_1,
                int diag, vector<bool> &lb_vector, FLOAT &cnt) {

    for (int col = 1; col <= temp_1; col++) {
        int row = diag + col - 1;

        if (lb_vector[col - 1]) {
            continue;
        }
        const vector<FLOAT > &t_ = subs[row - 1];
        const vector<FLOAT > &q = subs[col - 1];

        FLOAT d;
        FLOAT threshold2 = threshold * threshold;

        FLOAT x0 = t_[0];
        FLOAT y0 = t_[(m - 1)];

        const FLOAT dleft_orgin = DIST(x0, q[0]);
        const FLOAT dright_orgin = DIST(y0, q[m - 1]);
        FLOAT dleft = dleft_orgin;
        FLOAT dright = dright_orgin;

        FLOAT x1 = (t_[(1)]);
        const FLOAT d_left_weak = min(DIST(x1, q[0]), DIST(x1, q[1]));
        d = min(d_left_weak, DIST(x0, q[1]));
        dleft += d;

        FLOAT y1 = (t_[(m - 2)]);
        const FLOAT d_right_weak = min(DIST(y1, q[m - 1]), DIST(y1, q[m - 2]));
        d = min(d_right_weak, DIST(y0, q[m - 2]));
        dright += d;

        if (dleft + dright  >= threshold2) {
            cnt++;
            lb_vector[col] = true;
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
                lb_vector[col] = true;
                cnt++;
                continue;
            }
            special_shared_vector[col][0] = dleft;
            special_shared_vector[col][1] = dright;
        }
    }
}
