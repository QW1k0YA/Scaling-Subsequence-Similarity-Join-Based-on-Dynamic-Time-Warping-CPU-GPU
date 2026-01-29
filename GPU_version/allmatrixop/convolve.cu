
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../alldef/typedefdouble.cuh"
#include "../allunder/undermpx_v2.cuh"
#include "numeric"
using namespace std;

vector<FLOAT > convolve_valid(const std::vector<FLOAT >& vec, const std::vector<FLOAT >& core)
{
    vector<FLOAT > result(vec.size() - core.size() + 1, 0.0);
    vector<FLOAT > core_;
    unsigned long long core_size = core.size();
    for (int i = core_size; i >= 1; i--)
    {
        core_.push_back(core[i - 1]);
    }

    int j = 0;
    unsigned long long core_size_ = core_.size();
    auto temp_1 = vec.end() - core_size_+1;
    for (auto st_v = vec.begin(); st_v != temp_1; st_v++)
    {
        for (int i = 1; i <= core_size_; i++)
        {
            FLOAT  temp = core_[i - 1] * (*(st_v + i - 1));
            result[j] += temp;
        }
        j++;
    }

    return result;

}

