
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"
std::vector<DOUBLE> elementWiseDivison_vv(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match. in vectordivision" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    std::vector<DOUBLE> result(vector1_size);
    for (int i = 0; i < vector1_size; i++) {
        result[i] = (vector1[i] / vector2[i]);
    }

    return result;
}

void elementWiseDivison_vv(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE> &result) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match. in vectordivision" << std::endl;
    }

    size_t vector1_size = vector1.size();
    result.reserve(vector1_size);
    for (int i = 0; i < vector1_size; i++) {
        result[i] = (vector1[i] / vector2[i]);
    }

}