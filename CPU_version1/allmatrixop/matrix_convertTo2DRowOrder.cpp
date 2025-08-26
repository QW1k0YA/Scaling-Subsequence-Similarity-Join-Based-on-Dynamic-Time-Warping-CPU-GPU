
#include "../alldef/matrix.h"
#include <iostream>
#include <vector>

std::vector<std::vector<DOUBLE>> convertTo2DRowOrder(const std::vector<DOUBLE>& input)
{
    int numCols = input.size();
    int numRows = 1;

    std::vector<std::vector<DOUBLE>> result(1, std::vector<DOUBLE>(numCols));

    for (int i = 0; i < numCols; i++)
    {
        result[0][i] = input[i];
    }

    return result;
}
