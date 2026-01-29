
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../alldef/typedefdouble.cuh"

using namespace std;
std::vector<FLOAT > sortInd(const std::vector<FLOAT >& vec) {
    vector<FLOAT > result(vec.size());

    iota(result.begin(),result.end(),1);

    sort(result.begin(),result.end(),[&vec](int i,int j){return vec[i-1] < vec[j-1];});

    return result;
}

