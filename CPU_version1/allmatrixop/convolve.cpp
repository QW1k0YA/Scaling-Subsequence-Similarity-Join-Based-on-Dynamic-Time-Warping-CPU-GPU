
#include <vector>
#include <algorithm>  
#include <numeric>    
#include <functional> 
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.h"
#include "../alldef/typedefdouble.h"
#include "../allunder/undermpx_v2.h"
#include "numeric"
using namespace std;

vector<DOUBLE> convolve_valid(const std::vector<DOUBLE>& vec, const std::vector<DOUBLE>& core)
{
    vector<DOUBLE> result(vec.size() - core.size() + 1, 0.0);
    vector<DOUBLE> core_;
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
            DOUBLE temp = core_[i - 1] * (*(st_v + i - 1));
            result[j] += temp;
        }
        j++;
    }

    return result;

}

