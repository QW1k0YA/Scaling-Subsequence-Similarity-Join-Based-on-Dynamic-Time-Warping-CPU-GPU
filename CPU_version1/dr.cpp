
#include "vector"
#include "iostream"
#include "alldef/elseoperation.h"
#include "alldef/typedefdouble.h"
using namespace std;

vector<DOUBLE> dr(DOUBLE a, vector<DOUBLE>& ts, int st, int en)
{
    vector<DOUBLE> result;

    auto start = ts.begin() + st - 1;
    auto end = ts.begin() + en - 1;

    while(start != end + 1)
    {
        result.push_back(*start);
        start++;
    }

    auto first = result.begin();
    result.insert(first, a);

    return result;
}
