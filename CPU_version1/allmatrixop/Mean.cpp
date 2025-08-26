
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.h"
DOUBLE Mean(const std::vector<DOUBLE>& values) {
    if (values.empty()) {
        std::cerr << "Error: Vector is empty. Cannot calculate mean. in Mean" << std::endl;
        return 0.0; 
    }

    DOUBLE sum = 0.0;
    for (const DOUBLE& value : values) {
        sum += value;
    }

    return sum / static_cast<DOUBLE>(values.size());
}