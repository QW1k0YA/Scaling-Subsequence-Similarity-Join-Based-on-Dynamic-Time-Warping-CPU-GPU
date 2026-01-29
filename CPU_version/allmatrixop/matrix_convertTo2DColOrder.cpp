
#include "../alldef/matrix.h"
#include <iostream>
#include <vector>

std::vector<std::vector<DOUBLE>> convertTo2DColumnOrder(const std::vector<DOUBLE>& input) {
    int numRows = input.size();
    int numCols = 1;

    std::vector<std::vector<DOUBLE>> result( numRows,std::vector<DOUBLE>(1));

    for (int i = 0; i < numRows; i++)
    {
        result[i][0] = input[i];
    }

    return result;
}