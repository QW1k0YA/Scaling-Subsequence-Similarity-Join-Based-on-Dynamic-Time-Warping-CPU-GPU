
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.cuh"
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