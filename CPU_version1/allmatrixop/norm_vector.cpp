
#include "iostream"
#include "vector"
#include "../alldef/matrix.h"
#include "cmath"

using namespace std;

DOUBLE norm_vector(const vector<DOUBLE>& A, int b)
{
    if (b == 2)
    {
        DOUBLE sum = 0;
        size_t A_Size = A.size();
        for (int i = 0; i < A_Size; i++)
        {
            sum += pow(A[i], 2);
        }
        sum = pow(sum, 0.5);
        return sum;
    }
}