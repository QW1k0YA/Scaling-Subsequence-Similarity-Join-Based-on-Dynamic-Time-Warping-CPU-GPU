
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include "../alldef/typedefdouble.h"
std::vector<std::complex<DOUBLE>> conjugate(const std::vector<std::complex<DOUBLE>>& input) {
    std::vector<std::complex<DOUBLE>> result;
    result.reserve(input.size());
    std::transform(input.begin(), input.end(), std::back_inserter(result),
                   [](const std::complex<DOUBLE>& c) { return std::conj(c); });
    return result;
}