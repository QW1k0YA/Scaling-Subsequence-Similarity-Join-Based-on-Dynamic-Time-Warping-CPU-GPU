
#include<iostream>
#include<vector>
#include "../alldef/matrix.h"
using namespace std;

vector<DOUBLE> plusvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in plusvector";
    }
    else
    {
        vector<DOUBLE> result(a.size());
        size_t a_size = a.size();
        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }
}

void plusvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b, vector<DOUBLE>&result)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in plusvector";
    }
    else
    {
        size_t a_size = a.size();
        result.reserve(a_size);

        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i] + b[i];
        }
    }
}

DOUBLE* plusvector_p(const double *a, const double *b, long long len)
{
    double* result = new double[len];
    for (int i = 0; i < len;i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;

}

DOUBLE *plusvector_p_from_pos(const double *a, const double *b, long long len, long long pos)
{
    double* result = new double[len];
    for (int i = pos; i < len;i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;

}