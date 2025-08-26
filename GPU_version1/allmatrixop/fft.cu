
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "../alldef/typedefdouble.cuh"
std::vector<std::complex<FLOAT >> fft(const std::vector<FLOAT >& x, int n) {
    size_t N = x.size();

    std::vector<std::complex<FLOAT >> X(n, std::complex<FLOAT >(0, 0));
    for (int k = 0; k < n; ++k) {

        for (int j = 0; j < N; ++j) {
            X[k] += x[j] * std::polar(1.0f, -2 * 3.14159265f * k * j / N);
        }
    }
    return X;
}
