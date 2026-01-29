
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/typedefdouble.cuh"
/**
 * Checks each element in a vector for NaN (Not a Number) values.
 * This function returns a boolean vector indicating which elements are NaN.
 *
 * @param vec: Input vector of floating-point values
 * @return std::vector<bool>: Boolean vector where true indicates NaN values
 */
std::vector<bool> isNaN(const std::vector<FLOAT >& vec) {
    std::vector<bool> result;

    for (FLOAT  value : vec) {
        result.push_back(std::isnan(value));
    }

    return result;
}

