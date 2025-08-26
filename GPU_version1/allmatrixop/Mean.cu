
#include <iostream>
#include <vector>
#include "../alldef/typedefdouble.cuh"
FLOAT  Mean(const std::vector<FLOAT >& values) {
    if (values.empty()) {
        std::cerr << "Error: Vector is empty. Cannot calculate mean. in Mean" << std::endl;
        return 0.0; 
    }

    FLOAT  sum = 0.0;
    for (const FLOAT & value : values) {
        sum += value;
    }

    return sum / static_cast<FLOAT >(values.size());
}