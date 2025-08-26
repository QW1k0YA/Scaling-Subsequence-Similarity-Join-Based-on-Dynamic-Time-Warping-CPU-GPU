
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"

std::vector<std::vector<DOUBLE>> transposeMatrix_double(const std::vector<std::vector<DOUBLE>>& matrix) {
    
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    std::vector<std::vector<DOUBLE>> result(cols, std::vector<DOUBLE>(rows, 0.0));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}

std::vector<std::vector<bool>> transposeMatrix_bool(const std::vector<std::vector<bool>>& matrix) {
    
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();

    std::vector<std::vector<bool>> result(cols, std::vector<bool>(rows, false));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }

    return result;
}

