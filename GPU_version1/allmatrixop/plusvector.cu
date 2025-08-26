
#include<iostream>
#include<vector>
#include "../alldef/matrix.cuh"
using namespace std;

vector<FLOAT > plusvector(const vector<FLOAT >&a, const vector<FLOAT >&b)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in plusvector";
    }
    else
    {
        vector<FLOAT > result(a.size());
        size_t a_size = a.size();
        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i] + b[i];
        }
        return result;
    }
}