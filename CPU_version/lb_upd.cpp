
#include "alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "alldef/typedefdouble.h"

/**
 * LB_Keogh lower bound calculation with early abandoning for DTW distance estimation.
 * This function computes the LB_Keogh lower bound between a query sequence and an envelope,
 * with early abandoning when the accumulated distance exceeds the threshold.
 *
 * @param x: Query sequence
 * @param U: Upper envelope sequence
 * @param L: Lower envelope sequence
 * @param seqlen: Length of the sequences
 * @param bsf: Best-so-far distance threshold for early abandoning
 * @return: LB_Keogh lower bound distance, or a value >= bsf if early abandoned
 */
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

/**
 * Standard LB_Keogh lower bound calculation for DTW distance estimation.
 * This function computes the LB_Keogh lower bound between a query sequence and an envelope
 * without early abandoning optimization.
 *
 * @param x: Query sequence
 * @param U: Upper envelope sequence
 * @param L: Lower envelope sequence
 * @param seqlen: Length of the sequences
 * @return: LB_Keogh lower bound distance
 */
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

/**
 * Updated version of the LB_Keogh lower bound calculation for DTW distance estimation.
 * This function computes the LB_Keogh lower bound between a query sequence and an envelope,
 * with early abandoning when the accumulated distance exceeds the threshold.
 *
 * @param x: Query sequence
 * @param U: Upper envelope sequence
 * @param L: Lower envelope sequence
 * @param bsf: Best-so-far distance threshold for early abandoning (if 0, no early abandon)
 * @return: LB_Keogh lower bound distance, or a value >= bsf if early abandoned
 */
DOUBLE lb_upd(const std::vector<DOUBLE>& x, const std::vector<DOUBLE>& U, const std::vector<DOUBLE>& L, DOUBLE bsf) {
    return (bsf != 0.0) ? LB_Keogh_ea(x, U, L, static_cast<long long>(x.size()), bsf) : LB_Keogh(x, U, L, static_cast<long long>(x.size()));
}

