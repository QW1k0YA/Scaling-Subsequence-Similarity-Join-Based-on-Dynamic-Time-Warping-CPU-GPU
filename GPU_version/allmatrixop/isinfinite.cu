
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/typedefdouble.cuh"
/**
 * Checks each element in a vector for infinite values.
 * This function returns a boolean vector indicating which elements are infinite (+inf or -inf).
 *
 * @param vec: Input vector of floating-point values
 * @return std::vector<bool>: Boolean vector where true indicates infinite values
 */
std::vector<bool> isinfinite(const std::vector<FLOAT >& vec) {
    std::vector<bool> result;

    for (FLOAT  value : vec) {
        result.push_back(!std::isfinite(value));
    }

    return result;
}