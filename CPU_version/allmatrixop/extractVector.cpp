
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"

std::vector<DOUBLE> extractVector(const std::vector<std::vector<DOUBLE>>& matrix, size_t index, bool isRow) {
    std::vector<DOUBLE> result;

    size_t matrix_size = matrix.size();
    size_t matrix0_size = matrix[0].size();
    if (isRow) {
        if (index < matrix_size) {
            result = matrix[index];
        } else {
            std::cerr << "Error: Invalid row index. in extracVector" << std::endl;
        }
    } else {
        if (!matrix.empty() && index < matrix0_size) {
            for (const auto& row : matrix) {
                result.push_back(row[index]);
            }
        } else {
            std::cerr << "Error: Invalid column index.in extracVector" << std::endl;
        }
    }
    return result;
}