
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"
DOUBLE sum_matrix(const std::vector<std::vector<DOUBLE>>& matrix) {
    
    if (matrix.empty() || matrix[0].empty()) {
        
        return 0;
    }

    DOUBLE sum = 0;

    size_t matr_si = matrix.size();
    size_t matr0_si = matrix[0].size();
    for (size_t i = 0; i < matr_si; ++i) {
        for (size_t j = 0; j < matr0_si; ++j) {
            sum += static_cast<DOUBLE>(matrix[i][j]);
        }
    }

    return sum;
}