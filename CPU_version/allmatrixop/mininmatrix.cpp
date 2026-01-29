
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"

using namespace std;

DOUBLE mininmatrix(const vector<vector<DOUBLE>>& ma)
{
    DOUBLE result = ma[0][0];
    for(const vector<DOUBLE>& value1:ma)
    {
        for(DOUBLE value2:value1)
        {
            if(value2 < result)
            {
                result = value2;
            }
        }
    }

    return result;
}