
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include "alldef/typedefdouble.h"

void fft1(std::vector<std::complex<DOUBLE>>& a, bool inverse = false) {
    int n = a.size();
    if (n <= 1) return;
    std::vector<std::complex<DOUBLE>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    fft1(even, inverse);
    fft1(odd, inverse);

    DOUBLE angle = (inverse ? 1.0 : -1.0) * 2.0 * 3.1415926 / n;
    std::complex<DOUBLE> w(1.0, 0.0), wn(cos(angle), sin(angle));

    for (int i = 0; i < n / 2; ++i) {
        std::complex<DOUBLE> t = w * odd[i];
        a[i] = even[i] + t;
        a[i + n / 2] = even[i] - t;
        w *= wn;
    }
}

std::vector<DOUBLE> ifft(const std::vector<std::complex<DOUBLE>>& input, bool symmetric) {
    std::vector<std::complex<DOUBLE>> a = input;
    int n = a.size();

    if (symmetric) {
        fft1(a, true);
        for (int i = 0; i < n; ++i) {
            a[i] /= n;
        }
    }
    else {
        fft1(a, true);
    }

    std::vector<DOUBLE> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = a[i].real();
    }

    return result;
}