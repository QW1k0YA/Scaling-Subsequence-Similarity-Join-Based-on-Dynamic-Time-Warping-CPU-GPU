
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/allstruct.h"
#include "alldef/allclass.h"

using namespace std;

NEXT_RETURN alpha_generater::next(const vector<DOUBLE> &X, const vector<DOUBLE> &f_x, const vector<DOUBLE> &bsf_X,
                                  DOUBLE start_pos_, DOUBLE end_pos_) {
    NEXT_RETURN result;
    DOUBLE pos;
    DOUBLE start_pos = start_pos_;
    DOUBLE end_pos = end_pos_;

    pos = -1000;

    int len = static_cast<int>(X.size());
    if(vectorisempty(X))
    {
        pos = 0;
    }
    else if(fabs(X.size() - 1) < 0.001)
    {
        pos = -0.25;
    }
    else
    {
        for(int i = 1;i <= len;i++)
        {
            if((X[i-1] <= start_pos)||(X[i-1] >= end_pos) )
            {
                continue;
            }
            for(int j = 1;j <= len;j++)
            {
                if((X[j-1] <= start_pos) || (X[j-1] >= end_pos))
                {
                    continue;
                }
                if(bsf_X[i-1] > bsf_X[j-1]*1.05)
                {
                    continue;
                }

                if(f_x[j-1] < f_x[i-1])
                {
                    if((X[i-1]<=start_pos) || X[i-1] >= end_pos)
                    {
                        continue;
                    }

                    if(X[j-1] < X[i-1])
                    {
                        pos = X[j-1] + 0.5*(end_pos - X[j-1]);
                        start_pos = X[j-1];
                    }

                    if(X[j-1] > X[i-1])
                    {
                        pos = start_pos + 0.5*(X[i-1] - start_pos);
                        end_pos = X[j-1];
                    }

                }
            }
        }

        DOUBLE m = start_pos + end_pos;
        DOUBLE tttt;
        int fact;
        if(m <= 0)
        {
            fact = -1;
        }
        if(m > 0)
        {
            fact = 1;
        }

        if(fabs(pos+1000)<0.0001)
        {
            for(int i = 1;i <= 10;i++)
            {
                pos = i *0.25*fact;
                if((pos < start_pos) || (pos > end_pos))
                {
                    pos = -1000;
                    break;
                }

                tttt = 0.0;
                for(auto value:X)
                {
                    if(fabs(value - pos) < 0.0001)
                    {
                        tttt++;
                    }
                }

                if(fabs(tttt - 0) <0.0001)
                {
                    break;
                }
            }
        }

        if(fabs(pos + 1000) < 0.0001)
        {
            for(int i = 1;i <= 10;i ++)
            {
                pos = i*(-0.25)*fact;
                if((pos < start_pos) || (pos > end_pos))
                {
                    pos = -1000;
                    break;
                }

                tttt = 0.0;
                for(auto value:X)
                {
                    if(fabs(value - pos)< 0.0001)
                    {
                        tttt++;
                    }
                }

                if(fabs(tttt-0)<0.0001)
                {
                    break;
                }
            }
        }
    }

    result.end_pos = end_pos;
    result.start_pos = start_pos;
    result.pos = pos;

    return result;
}

