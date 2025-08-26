
#include "vector"
#include "iostream"
#include "alldef/elseoperation.h"
#include "alldef/typedefdouble.h"
using namespace std;

vector<DOUBLE> dr_f(const std::vector<DOUBLE>& ts, int st)
{
    vector<DOUBLE> result;

    auto start = ts.begin() + st - 1;
    auto end = ts.end();

    while (start != end )
    {
        result.push_back(*start);
        start++;
    }

    return result;
}
