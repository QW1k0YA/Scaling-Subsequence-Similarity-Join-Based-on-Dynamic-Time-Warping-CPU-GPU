
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"

using namespace std;

/**
 * Check if a vector is empty.
 * This function determines whether the input vector has zero elements.
 *
 * @param X: Input vector of DOUBLE values
 * @return: True if the vector is empty (size is 0), false otherwise
 */
bool vectorisempty(const vector<DOUBLE>& X)
{

    size_t X_size = X.size();
    if(X_size<0.1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}