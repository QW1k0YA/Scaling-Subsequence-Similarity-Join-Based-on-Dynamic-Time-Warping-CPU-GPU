
#include <iostream>
#include <vector>
#include "../alldef/matrix.cuh"
std::vector<FLOAT > elementWiseMultiply(const std::vector<FLOAT >& vector1, const std::vector<FLOAT >& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    std::vector<FLOAT > result(vector1_size);
    for (int i = 0; i < vector1_size; i++) {
        result[i] = (vector1[i] * vector2[i]);
    }

    return result;
}

FLOAT  elementWiseMultiply_sum(const std::vector<FLOAT >& vector1, const std::vector<FLOAT >& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply_sum" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    FLOAT  sum = 0;
    for (int i = 0; i < vector1_size; i++) {
        sum += (vector1[i] * vector2[i]);
    }

    return sum;
}

std::vector<FLOAT > elementWiseMultiply_p(const FLOAT  vector1[], const FLOAT  vector2[], long long len) {

    std::vector<FLOAT > result(len);

    for (int i = 0; i < len; i++) {
        result[i]=(vector1[i] * vector2[i]);
    }
    return result;
}

__device__ FLOAT  elementWiseMultiply_p_sum(const FLOAT  vector1[], const FLOAT  vector2[], long long len) {

    FLOAT  result = 0;

    for (int i = 0; i < len; i++) {
        result+=(vector1[i] * vector2[i]);
    }
    return result;
}

__device__ FLOAT  elementWiseMultiply_p_plus_sum(const FLOAT  *vector1, const FLOAT  *vector2, const FLOAT  *vector3, long long len) {

    FLOAT  result = 0;

    for (int i = 0; i < len; i++) {
        result += (vector1[i] * (vector2[i] + vector3[i]));
    }
    return result;
}
