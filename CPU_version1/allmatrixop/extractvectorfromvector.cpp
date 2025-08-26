
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"

using namespace std;

vector<DOUBLE> extr_vfromv(const std::vector<DOUBLE> &arr, int st, int ed) {
    if (st > ed)
    {
        cerr << "extr_vfromv has a error in extravectorfromvector" << endl;
        exit(0);
    }

    vector<DOUBLE> result( (ed - st + 1), 0.0);
    memcpy(&result[0], &arr[0] + (st - 1), sizeof(DOUBLE) * (ed - st + 1));

    return result;
}