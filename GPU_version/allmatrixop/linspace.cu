
#include <vector>
#include "../alldef/typedefdouble.cuh"
std::vector<FLOAT > linspace(FLOAT  start, FLOAT  end, int num) {
    std::vector<FLOAT > result;
    if (num <= 1) {
        result.push_back(start);
        return result;
    }

    FLOAT  step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
}
