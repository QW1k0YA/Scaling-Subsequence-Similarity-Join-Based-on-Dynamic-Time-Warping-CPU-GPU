
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.h"

DOUBLE dotProduct(const std::vector<DOUBLE>& A, const std::vector<DOUBLE>& B) {

    if (A.size() != B.size()) {
        std::cerr << "Error: Vectors must have the same size. in dotproduct" << std::endl;
        return 0.0;
    }

    DOUBLE result = 0.0;
    size_t temp = A.size();
    for (std::size_t i = 0; i < temp; ++i) {
        result += A[i] * B[i];
    }

    return result;
}