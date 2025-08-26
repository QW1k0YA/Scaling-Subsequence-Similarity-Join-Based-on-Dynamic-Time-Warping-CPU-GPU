
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.h"
std::vector<DOUBLE> mergeVectors(const std::vector<DOUBLE>& v1, const std::vector<DOUBLE>& v2, const std::vector<DOUBLE>& v3) {
    std::vector<DOUBLE> result;

    result.insert(result.end(), v1.begin(), v1.end());
    result.insert(result.end(), v2.begin(), v2.end());
    result.insert(result.end(), v3.begin(), v3.end());

    return result;
}