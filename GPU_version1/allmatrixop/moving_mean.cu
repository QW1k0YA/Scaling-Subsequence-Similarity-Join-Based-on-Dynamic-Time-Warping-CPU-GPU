
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../alldef/typedefdouble.cuh"
using namespace std;

vector<FLOAT > moving_mean(vector<FLOAT > a, int w)
{
    vector<FLOAT > res(a.size() - w + 1, 0.0);
    FLOAT  p = a[0];
    FLOAT  s = 0;
    FLOAT  x,z;
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