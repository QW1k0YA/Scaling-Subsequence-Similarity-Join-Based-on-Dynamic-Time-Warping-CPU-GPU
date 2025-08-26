
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include "../alldef/typedefdouble.cuh"
std::vector<std::complex<FLOAT >> conjugate(const std::vector<std::complex<FLOAT >>& input) {
    std::vector<std::complex<FLOAT >> result;
    result.reserve(input.size());
    std::transform(input.begin(), input.end(), std::back_inserter(result),
                   [](const std::complex<FLOAT >& c) { return std::conj(c); });
    return result;
}