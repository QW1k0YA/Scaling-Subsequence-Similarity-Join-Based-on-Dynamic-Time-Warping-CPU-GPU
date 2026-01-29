
#include "iostream"
#include "vector"
#include "../alldef/matrix.cuh"
#include "cmath"

using namespace std;

/**
 * Computes the L2 (Euclidean) norm of a vector.
 * This function calculates the square root of the sum of squares of vector elements.
 *
 * @param A: Input vector
 * @param b: Norm type (currently only supports b=2 for L2 norm)
 * @return FLOAT: L2 norm of the vector
 */
FLOAT  norm_vector(const vector<FLOAT >& A, int b)
{
    if (b == 2)
    {
        FLOAT  sum = 0;
        size_t A_Size = A.size();
        for (int i = 0; i < A_Size; i++)
        {
            sum += pow(A[i], 2);
        }
        sum = pow(sum, 0.5);
        return sum;
    }
}