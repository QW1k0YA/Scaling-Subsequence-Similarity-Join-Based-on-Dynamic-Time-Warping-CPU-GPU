
#include <vector>
#include "../alldef/matrix.cuh"
#include "iostream"
FLOAT  sum_vector(const std::vector<FLOAT >& vec) {

    FLOAT  sum = 0.0;
    int size = vec.size();
    for (int i = 0;i < size;i++) {
        sum += vec[i];
    }
    return sum;
}
