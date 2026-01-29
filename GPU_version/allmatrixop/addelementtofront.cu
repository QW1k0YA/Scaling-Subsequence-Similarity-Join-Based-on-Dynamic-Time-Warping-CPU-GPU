
#include <iostream>
#include <vector>
#include "../alldef/matrix.cuh"
std::vector<FLOAT > addElementToFront(const std::vector<FLOAT >& arr, FLOAT  element) {

    std::vector<FLOAT > newArr = {element };
    newArr.insert(newArr.end(), arr.begin(), arr.end());

    return newArr;
}