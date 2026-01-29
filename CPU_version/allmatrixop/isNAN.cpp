
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/typedefdouble.h"
/**
 * Check for NaN values in a vector.
 * This function creates a boolean vector indicating which elements of the input vector are NaN.
 *
 * @param vec: Input vector of DOUBLE values
 * @return: Boolean vector where true indicates the corresponding element is NaN
 */
std::vector<bool> isNaN(const std::vector<DOUBLE>& vec) {
    std::vector<bool> result;

    for (DOUBLE value : vec) {
        result.push_back(std::isnan(value));
    }

    return result;
}

