
#include <vector>
#include "../alldef/matrix.h"
#include "iostream"
DOUBLE sum_vector(const std::vector<DOUBLE>& vec) {

    DOUBLE sum = 0.0;
    int size = vec.size();
    for (int i = 0;i < size;i++) {
        sum += vec[i];
    }
    return sum;
}
