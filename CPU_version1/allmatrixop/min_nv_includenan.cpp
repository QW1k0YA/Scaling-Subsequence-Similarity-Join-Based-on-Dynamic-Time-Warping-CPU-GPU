
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "../alldef/typedefdouble.h"
#include "../alldef/elseoperation.h"
#include "../alldef/matrix.h"
using namespace std;
std::vector<DOUBLE> min_nv_Includenan(DOUBLE x, const std::vector<DOUBLE>& data) {

    std::vector<DOUBLE> result;

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

std::vector<DOUBLE> max_nv_Includenan(DOUBLE x, const std::vector<DOUBLE>& data) {

    std::vector<DOUBLE> result;

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
