
#include <vector>
#include "../alldef/matrix.h"
#include "iostream"
int sum_vector_bool(const std::vector<bool>& vec) {
    int sum = 0.0;
    for (const auto& element : vec) {
        sum += element;
    }

    return sum;
}
