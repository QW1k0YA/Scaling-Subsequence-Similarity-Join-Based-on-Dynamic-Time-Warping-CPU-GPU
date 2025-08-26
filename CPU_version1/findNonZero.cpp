
#include <iostream>
#include <vector>
#include "alldef/typedefdouble.h"
std::vector<DOUBLE> findNonZero(const std::vector<bool>& vec) {
    std::vector<DOUBLE> result;

    size_t vec_size = vec.size();
    for (size_t i = 0; i < vec_size; ++i) {
        if (vec[i] != 0) {
            result.push_back(static_cast<DOUBLE>(i)+1); 
        }
    }

    return result;
}