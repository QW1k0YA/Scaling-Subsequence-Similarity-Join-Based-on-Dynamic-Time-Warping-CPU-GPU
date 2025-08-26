
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"

std::vector<DOUBLE> movmin(const std::vector<DOUBLE>& A, int kb, int kf) {

    vector<DOUBLE> result;
    vector<DOUBLE> a(kb + kf + A.size(), INFINITY);
    memcpy(&a[0] + kb , &A[0], sizeof(DOUBLE) * A.size());

    vector<DOUBLE> q(a.size(),0);

    int k = kb + kf + 1;
    int n = a.size();

    int hh = 0, tt = -1;

    for (int i = 0; i < n; i++) {

        if (hh <= tt && i - k + 1 > q[hh]) hh++;
        while (hh <= tt && a[q[tt]] >= a[i]) tt--;
        q[++tt] = i;
        if (i >= k - 1)
        {
            result.push_back(a[q[hh]]);
        }
    }
    return result;

}