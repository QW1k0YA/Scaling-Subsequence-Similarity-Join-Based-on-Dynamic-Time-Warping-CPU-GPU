
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"

std::vector<std::vector<bool>> falseMatrix(size_t rows, size_t cols) {
    return std::vector<std::vector<bool>>(rows, std::vector<bool>(cols, false));
}