
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "../alldef/typedefdouble.cuh"
#include "../alldef/elseoperation.cuh"
#include "../alldef/matrix.cuh"
using namespace std;
std::vector<FLOAT > min_nv_Includenan(FLOAT  x, const std::vector<FLOAT >& data) {

    std::vector<FLOAT > result;

    for (auto value : data)
    {
        if (isnan(value))
        {
            result.push_back(NAN);
        }
        else
        {
            result.push_back(MIN(x, value));
        }
    }

    return result;
}

std::vector<FLOAT > max_nv_Includenan(FLOAT  x, const std::vector<FLOAT >& data) {

    std::vector<FLOAT > result;

    for (auto value : data)
    {
        if (isnan(value))
        {
            result.push_back(NAN);
        }
        else
        {
            result.push_back(MAX(x, value));
        }
    }

    return result;
}
