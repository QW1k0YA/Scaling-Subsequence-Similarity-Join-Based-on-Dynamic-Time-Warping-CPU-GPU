
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"

using namespace std;

void normalise(vector<double>&A,vector<double> &B,double mean,double std)
{
    int len = A.size();
    for(int i = 0;i < len;i++)
    {
        B[i] = (A[i] - mean)/std;
    }
}