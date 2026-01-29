
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.cuh"
/**
 * Finds the indices of non-zero (true) elements in a boolean vector.
 * This function returns the 1-indexed positions of all true values in the input vector.
 *
 * @param vec: Input boolean vector
 * @return std::vector<FLOAT>: Vector containing 1-indexed positions of true values
 */
std::vector<FLOAT > findNonZero(const std::vector<bool>& vec) {
    std::vector<FLOAT > result;

    size_t vec_size = vec.size();
    for (size_t i = 0; i < vec_size; ++i) {
        if (vec[i] != 0) {
            result.push_back(static_cast<FLOAT >(i) + 1); 
        }
    }

    return result;
}