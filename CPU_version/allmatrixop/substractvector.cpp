
#include "../alldef/matrix.h"
#include<iostream>
#include<vector>

using namespace std;

vector<DOUBLE> substractvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in subvector";
    }
    else
    {
        vector<DOUBLE> result(a.size());
        size_t a_size = a.size();
        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }
}

void substractvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b,vector<DOUBLE> &result)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in subvector";
    }
    else
    {
        result.reserve(a.size());
        size_t a_size = a.size();
        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i] - b[i];
        }
    }
}

void substractvector_from_pos(const vector<DOUBLE>&a, const vector<DOUBLE>&b,vector<DOUBLE> &result,int pos)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in subvector";
    }
    else
    {
        result.reserve(a.size());
        size_t a_size = a.size();
        for (int i = pos; i < a_size;i++)
        {
            result[i] = a[i] - b[i];
        }
    }
}
void substractvector(const vector<DOUBLE>&a, const vector<DOUBLE>&b,vector<DOUBLE> &result,int len)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in subvector";
    }
    else
    {
        result.reserve(len);

        for (int i = 0; i < len;i++)
        {
            result[i] = a[i] - b[i];
        }
    }
}

void subs_vector_p(const double *A,const double *B,double *C,long long len)
{
    for(int i = 0;i < len;i++)
    {
        C[i] = A[i] - B[i];
    }
}
