
#include "iostream"
#include "vector"
#include "../alldef/matrix.h"
#include "cmath"

using namespace std;

/**
 * Compute the L2 norm (Euclidean norm) of a vector.
 * This function calculates the square root of the sum of squares of vector elements.
 *
 * @param A: Input vector of DOUBLE values
 * @param b: Norm type (only supports b=2 for L2 norm)
 * @return: L2 norm of the input vector
 */
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