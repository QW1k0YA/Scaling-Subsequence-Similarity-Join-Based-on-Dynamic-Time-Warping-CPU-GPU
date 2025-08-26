
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/typedefdouble.h"
std::vector<bool> isNaN(const std::vector<DOUBLE>& vec) {
    std::vector<bool> result;

    for (DOUBLE value : vec) {
        result.push_back(std::isnan(value));
    }

    return result;
}

