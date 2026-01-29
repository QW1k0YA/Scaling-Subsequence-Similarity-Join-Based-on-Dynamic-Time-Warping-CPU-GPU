
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.cuh"
std::vector<FLOAT > elementwiseMultiply_nv(FLOAT  scalar, const std::vector<FLOAT >& vector) {
    std::vector<FLOAT > result;
    result.reserve(vector.size());

    for (const auto& element : vector) {
        result.push_back(scalar * element);
    }

    return result;
}