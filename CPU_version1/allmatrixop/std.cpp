
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/matrix.h"
DOUBLE stddev(const std::vector<DOUBLE>& data) {

    DOUBLE mean = 0.0;
    for (DOUBLE value : data) {
        mean += value;
    }
    mean /= data.size();

    DOUBLE sumSquaredDifferences = 0.0;
    for (DOUBLE value : data) {
        sumSquaredDifferences += (value - mean) * (value - mean);
    }

    DOUBLE variance = sumSquaredDifferences / (data.size());
    
    DOUBLE standardDeviation = std::sqrt(variance);

    return standardDeviation;
}