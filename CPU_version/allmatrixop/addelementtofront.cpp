
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"
std::vector<DOUBLE> addElementToFront(const std::vector<DOUBLE>& arr, DOUBLE element) {

    std::vector<DOUBLE> newArr = { element };
    newArr.insert(newArr.end(), arr.begin(), arr.end());

    return newArr;
}