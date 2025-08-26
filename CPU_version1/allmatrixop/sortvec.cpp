
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../alldef/typedefdouble.h"

using namespace std;
std::vector<DOUBLE> sortInd(const std::vector<DOUBLE>& vec) {
    vector<DOUBLE> result(vec.size());

    iota(result.begin(),result.end(),1);

    sort(result.begin(),result.end(),[&vec](int i,int j){return vec[i-1] < vec[j-1];});

    return result;
}

