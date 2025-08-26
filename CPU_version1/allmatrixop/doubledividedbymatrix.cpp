
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"

std::vector<std::vector<DOUBLE>> elementWiseDivision(DOUBLE scalar, const std::vector<std::vector<DOUBLE>>& matrix) {

    std::vector<std::vector<DOUBLE>> result(matrix.size(), std::vector<DOUBLE>(matrix[0].size(), 0.0));

   size_t temp_1 = matrix.size();
   size_t temp_2 = matrix[0].size();

    for (size_t i = 0; i < temp_1; ++i) {
        for (size_t j = 0; j < temp_2; ++j) {

            if (matrix[i][j] != 0.0) {
                result[i][j] = scalar / matrix[i][j];
            }
            else {
                std::cerr << "Error: Division by zero encountered. in doubledividedbymatrix" << std::endl;
                return {};
            }
        }
    }

    return result;
}
