
#include <iostream>
#include <vector>
#include <complex>
#include "../alldef/typedefdouble.h"
std::vector<std::complex<DOUBLE>> elementWiseMultiply_complex(const std::vector<std::complex<DOUBLE>>& vec1,
                                                            const std::vector<std::complex<DOUBLE>>& vec2) {
    std::vector<std::complex<DOUBLE>> result;
    result.reserve(vec1.size());

    size_t vec1_size = vec1.size();
    for (size_t i = 0; i < vec1_size; ++i) {
        result.push_back(vec1[i] * vec2[i]);
    }
    return result;
}

