
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"
std::vector<DOUBLE> elementWiseMultiply(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    std::vector<DOUBLE> result(vector1_size);
    for (int i = 0; i < vector1_size; i++) {
        result[i] = (vector1[i] * vector2[i]);
    }

    return result;
}

std::vector<DOUBLE> elementWiseMultiply_from_pos(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2,int pos) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    std::vector<DOUBLE> result(vector1_size);
    for (int i = pos; i < vector1_size; i++) {
        result[i] = (vector1[i] * vector2[i]);
    }

    return result;
}

DOUBLE elementWiseMultiply_sum(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2) {
    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply_sum" << std::endl;
        return {};
    }

    size_t vector1_size = vector1.size();
    DOUBLE sum = 0;
    for (int i = 0; i < vector1_size; i++) {
        sum += (vector1[i] * vector2[i]);
    }

    return sum;
}

void elementWiseMultiply(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE>& result) {

    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply" << std::endl;
    }

    size_t vector1_size = vector1.size();
    result.reserve(vector1_size);
    for (int i = 0; i < vector1_size; i++) {
        result[i] = (vector1[i] * vector2[i]);
    }

}
void elementWiseMultiply_from_pos(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE>& result,int pos) {

    if (vector1.size() != vector2.size()) {
        std::cerr << "Error: Vector sizes do not match.in vectormultiply" << std::endl;
    }

    size_t vector1_size = vector1.size();
    result.reserve(vector1_size);
    for (int i = pos; i < vector1_size; i++) {
        result[i] = (vector1[i] * vector2[i]);
    }

}

void elementWiseMultiply(const std::vector<DOUBLE>& vector1, const std::vector<DOUBLE>& vector2, vector<DOUBLE>& result,int len) {

    len = MIN(vector1.size(), MIN(vector2.size(),len));
    result.reserve(len);
    for (int i = 0; i < len; i++) {
        result[i] = (vector1[i] * vector2[i]);
    }

}
std::vector<DOUBLE> elementWiseMultiply_p(const DOUBLE vector1[], const DOUBLE vector2[], long long len) {

    std::vector<DOUBLE> result(len);

    for (int i = 0; i < len; i++) {
        result[i]=(vector1[i] * vector2[i]);
    }
    return result;
}

DOUBLE elementWiseMultiply_p_sum(const DOUBLE vector1[], const DOUBLE vector2[], long long len) {

    DOUBLE result = 0;

    for (int i = 0; i < len; i++) {
        result+=(vector1[i] * vector2[i]);
    }
    return result;
}

std::vector<DOUBLE> elementWiseMultiply_p_plus(const DOUBLE vector1[], const DOUBLE vector2[],const DOUBLE vector3[] ,long long len) {

    std::vector<DOUBLE> result(len);

    for (int i = 0; i < len; i++) {
        result[i] = (vector1[i] * (vector2[i] + vector3[i]));
    }
    return result;
}

DOUBLE elementWiseMultiply_p_plus_sum(const DOUBLE *vector1, const DOUBLE *vector2, const DOUBLE *vector3, long long len) {

    DOUBLE result = 0;

    for (int i = 0; i < len; i++) {
        result += (vector1[i] * (vector2[i] + vector3[i]));
    }
    return result;
}