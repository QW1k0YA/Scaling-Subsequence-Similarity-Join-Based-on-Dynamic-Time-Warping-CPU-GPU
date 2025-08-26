
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"

using namespace std;

vector<FLOAT > extr_vfromv(const std::vector<FLOAT > &arr, int st, int ed) {
    if (st > ed)
    {
        cerr << "extr_vfromv has a error in extravectorfromvector" << endl;
        exit(0);
    }

    vector<FLOAT > result((ed - st + 1), 0.0);
    memcpy(&result[0], &arr[0] + (st - 1), sizeof(FLOAT ) * (ed - st + 1));

    return result;
}