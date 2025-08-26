
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include "../alldef/matrix.h"

DOUBLE meanWithoutNaN(const std::vector<DOUBLE>& vec) {
    DOUBLE sum = 0.0;
    int count = 0;

    for (const auto& value : vec) {

        if (!std::isnan(value)) {

            sum += value;
            count++;
        }
    }

    return (count == 0) ? std::nan("") : sum / count;
}

DOUBLE meanWithoutNaN_from_pos(const std::vector<DOUBLE>& vec,int pos) {
    DOUBLE sum = 0.0;
    int count = 0;

    for (const auto& value : vec) {

        if (!std::isnan(value)) {

            sum += value;
            count++;
        }
    }

    return (count == 0) ? std::nan("") : sum / count;
}