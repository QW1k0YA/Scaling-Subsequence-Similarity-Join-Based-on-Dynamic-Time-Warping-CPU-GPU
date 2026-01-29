
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../alldef/typedefdouble.h"

using namespace std;
/**
 * Get the indices that would sort a vector in ascending order.
 * This function returns the indices that would sort the input vector, maintaining
 * the original values' relationships.
 *
 * @param vec: Input vector of DOUBLE values
 * @return: Vector containing 1-indexed positions that would sort the input vector
 */
std::vector<DOUBLE> sortInd(const std::vector<DOUBLE>& vec) {
    vector<DOUBLE> result(vec.size());

    iota(result.begin(),result.end(),1);

    sort(result.begin(),result.end(),[&vec](int i,int j){return vec[i-1] < vec[j-1];});

    return result;
}

