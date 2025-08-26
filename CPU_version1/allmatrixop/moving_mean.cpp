
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.h"
#include "../alldef/typedefdouble.h"
using namespace std;

vector<DOUBLE> moving_mean(vector<DOUBLE> a,int w)
{
    vector<DOUBLE> res(a.size() - w +1,0.0);
    DOUBLE p = a[0];
    DOUBLE s = 0;
    DOUBLE x,z;
    for(int i = 2;i <= w ;i++)
    {
        x = p +a[i-1];
        z = x -p;
        s+= (p-x+z) + a[i-1] -z;
        p =x;
    }

    res[0] = p + s;
    size_t a_size = a.size();
    for(int i = w+1;i <= a_size;i++)
    {
        x = p - a[i-w-1];
        z = x - p;
        s += p - x + z - a[i-w-1]-z;
        p = x;

        x = p + a[i-1];
        z = x - p;
        s += p - x +z + a[i-1] - z;
        p = x;

        res[i -w] = p +s;
    }

    size_t res_size = res.size();
    for(int i = 1;i <= res_size;i++)
    {
        res[i-1] = res[i-1]/w;
    }

    return res;
}