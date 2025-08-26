
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/typedefdouble.h"

DOUBLE LB_Keogh_ea(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, long long seqlen, DOUBLE bsf) {
    DOUBLE dist = 0.0;
    for (long long i = 0; i < seqlen; i++) {
        if (x[i] > U[i]) {
            dist += (x[i] - U[i]) * (x[i] - U[i]);
        } else if (x[i] < L[i]) {
            dist += (x[i] - L[i]) * (x[i] - L[i]);
        }
        if (dist >= bsf) {
            return dist;
        }
    }
    return dist;
}

DOUBLE LB_Keogh(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, long long seqlen) {
    DOUBLE dist = 0.0;
    for (long long i = 0; i < seqlen; i++) {
        if (x[i] > U[i]) {
            dist += (x[i] - U[i]) * (x[i] - U[i]);
        } else if (x[i] < L[i]) {
            dist += (x[i] - L[i]) * (x[i] - L[i]);
        }
    }
    return dist;
}

DOUBLE lb_upd(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, DOUBLE bsf) {
    return (bsf != 0.0) ? LB_Keogh_ea(x, U, L, static_cast<long long>(x.size()), bsf) : LB_Keogh(x, U, L, static_cast<long long>(x.size()));
}

