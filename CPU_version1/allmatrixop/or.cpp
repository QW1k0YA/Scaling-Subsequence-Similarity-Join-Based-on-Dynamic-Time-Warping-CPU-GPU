
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"
std::vector<std::vector<bool>> matrixOr(const std::vector<std::vector<bool>>& matrix1, const std::vector<std::vector<bool>>& matrix2) {
    
    if (matrix1.empty() || matrix2.empty() || matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
        
        std::cerr<<"size false in matrixor"<<std::endl;
        exit(0);
    }

    std::vector<std::vector<bool>> result(matrix1.size(), std::vector<bool>(matrix1[0].size(), false));

    size_t matrix1si = matrix1.size();
    size_t matrix10_si = matrix1.size();
    for (size_t i = 0; i < matrix1si; ++i) {
        for (size_t j = 0; j < matrix10_si; ++j) {
            result[i][j] = matrix1[i][j] || matrix2[i][j];
        }
    }

    return result;
}