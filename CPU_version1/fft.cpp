
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "alldef/typedefdouble.h"
std::vector<std::complex<DOUBLE>> fft(const std::vector<DOUBLE>& x, int n) {
    size_t N = x.size();

    std::vector<std::complex<DOUBLE>> X(n, std::complex<DOUBLE>(0, 0));
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < N; ++j) {
            X[k] += x[j] * std::polar(1.0, -2 * 3.1415926 * k * j / N);
        }
    }
    return X;
}
