
#include "../alldef/matrix.cuh"
#include<iostream>
#include<vector>

using namespace std;

vector<FLOAT > substractvector(const vector<FLOAT >&a, const vector<FLOAT >&b)
{
    if (a.size()!= b.size())
    {
        cerr << "a.size() is not equal to b.size() in subvector";
    }
    else
    {
        vector<FLOAT > result(a.size());
        size_t a_size = a.size();
        for (int i = 0; i < a_size;i++)
        {
            result[i] = a[i] - b[i];
        }
        return result;
    }
}
