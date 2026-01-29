
#include <vector>
#include "alldef/typedefdouble.h"
std::vector<DOUBLE> linspace(DOUBLE start, DOUBLE end, int num) {
    std::vector<DOUBLE> result;
    if (num <= 1) {
        result.push_back(start);
        return result;
    }

    DOUBLE step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result.push_back(start + i * step);
    }

    return result;
}
