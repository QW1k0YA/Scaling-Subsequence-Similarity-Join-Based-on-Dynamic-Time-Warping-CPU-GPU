
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"

using namespace std;

bool vectorisempty(const vector<DOUBLE>& X)
{

    size_t X_size = X.size();
    if(X_size<0.1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}