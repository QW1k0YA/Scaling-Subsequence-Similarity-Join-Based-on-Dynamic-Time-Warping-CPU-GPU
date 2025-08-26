
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.h"
std::vector<DOUBLE> elementwiseMultiply_nv(DOUBLE scalar, const std::vector<DOUBLE>& vector) {
    std::vector<DOUBLE> result;
    result.reserve(vector.size());

    for (const auto& element : vector) {
        result.push_back(scalar * element);
    }

    return result;
}