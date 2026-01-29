
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/typedefdouble.h"
/**
 * Check for infinite values in a vector.
 * This function creates a boolean vector indicating which elements of the input vector are infinite.
 *
 * @param vec: Input vector of DOUBLE values
 * @return: Boolean vector where true indicates the corresponding element is infinite (positive or negative)
 */
std::vector<bool> isinfinite(const std::vector<DOUBLE>& vec) {
    std::vector<bool> result;

    for (DOUBLE value : vec) {
        result.push_back(!std::isfinite(value));
    }

    return result;
}