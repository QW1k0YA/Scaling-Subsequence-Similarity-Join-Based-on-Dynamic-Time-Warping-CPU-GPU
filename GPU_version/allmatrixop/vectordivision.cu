
#include <iostream>
#include <vector>
#include "../alldef/matrix.cuh"
std::vector<FLOAT > elementWiseDivison_vv(const std::vector<FLOAT >& vector1, const std::vector<FLOAT >& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match. in vectordivision" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    std::vector<FLOAT > result(vector1_size);
    for (int i = 0; i < vector1_size; i++) {
        result[i] = (vector1[i] / vector2[i]);
    }

    return result;
}
