
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/typedefdouble.cuh"
std::vector<bool> isinfinite(const std::vector<FLOAT >& vec) {
    std::vector<bool> result;

    for (FLOAT  value : vec) {
        result.push_back(!std::isfinite(value));
    }

    return result;
}