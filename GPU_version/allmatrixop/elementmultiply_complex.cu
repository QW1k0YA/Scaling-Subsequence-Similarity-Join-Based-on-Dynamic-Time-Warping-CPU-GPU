
#include <iostream>
#include <vector>
#include <complex>
#include "../alldef/typedefdouble.cuh"
std::vector<std::complex<FLOAT >> elementWiseMultiply_complex(const std::vector<std::complex<FLOAT >>& vec1,
                                                             const std::vector<std::complex<FLOAT >>& vec2) {
    std::vector<std::complex<FLOAT >> result;
    result.reserve(vec1.size());

    size_t vec1_size = vec1.size();
    for (size_t i = 0; i < vec1_size; ++i) {
        result.push_back(vec1[i] * vec2[i]);
    }
    return result;
}

