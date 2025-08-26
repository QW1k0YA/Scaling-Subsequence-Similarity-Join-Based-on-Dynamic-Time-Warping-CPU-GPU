
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "../alldef/typedefdouble.h"
DOUBLE stdWithOmitnan(const std::vector<DOUBLE>& A) {
    
    std::vector<DOUBLE> filteredA;
    std::copy_if(A.begin(), A.end(), std::back_inserter(filteredA), [](DOUBLE val) {
        return !std::isnan(val);
    });

    DOUBLE mean = std::accumulate(filteredA.begin(), filteredA.end(), 0.0) / filteredA.size();

    DOUBLE squaredSum = std::accumulate(filteredA.begin(), filteredA.end(), 0.0, [mean](DOUBLE sum, DOUBLE val) {
        return sum + (val - mean) * (val - mean);
    });

    DOUBLE stdDev = std::sqrt(squaredSum / (filteredA.size()-1));

    return stdDev;
}