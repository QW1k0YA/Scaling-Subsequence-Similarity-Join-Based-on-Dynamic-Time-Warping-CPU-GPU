
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.cuh"
std::vector<FLOAT > mergeVectors(const std::vector<FLOAT >& v1, const std::vector<FLOAT >& v2, const std::vector<FLOAT >& v3) {
    std::vector<FLOAT > result;

    result.insert(result.end(), v1.begin(), v1.end());
    result.insert(result.end(), v2.begin(), v2.end());
    result.insert(result.end(), v3.begin(), v3.end());

    return result;
}