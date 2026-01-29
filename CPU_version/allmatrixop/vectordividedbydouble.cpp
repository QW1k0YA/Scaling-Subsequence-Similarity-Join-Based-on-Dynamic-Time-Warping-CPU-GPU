
#include<iostream>
#include<vector>
#include "../alldef/matrix.h"
using namespace std;

vector<DOUBLE> vectordividedbydouble(const vector<DOUBLE>&a, DOUBLE b)
{
    {
        vector<DOUBLE> result(a.size());
        size_t a_size = a.size();
        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i]/b;
        }
        return result;
    }
}