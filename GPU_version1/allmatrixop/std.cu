
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/matrix.cuh"
FLOAT  stddev(const std::vector<FLOAT >& data) {

    FLOAT  mean = 0.0;
    for (FLOAT  value : data) {
        mean += value;
    }
    mean /= data.size();

    FLOAT  sumSquaredDifferences = 0.0;
    for (FLOAT  value : data) {
        sumSquaredDifferences += (value - mean) * (value - mean);
    }

    FLOAT  variance = sumSquaredDifferences / (data.size());
    
    FLOAT  standardDeviation = std::sqrt(variance);

    return standardDeviation;
}