
#include <queue>
#include "../alldef/matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.h"
#include "../alldef/typedefdouble.h"
#include "../allunder/undermpx_v2.h"

using namespace std;

RETURN_MIN min_v(const std::vector<DOUBLE>& v) {
    RETURN_MIN result;
    if (v.empty()) {
        result.minval = 0.0;
        result.minidx = 0.0;
        return result;
    }

    DOUBLE value = v[0];
    DOUBLE idx = 0.0;

    size_t v_size = v.size();
    for (int i = 1; i < v_size; ++i) {
        if (v[i] < value) {
            value = v[i];
            idx = static_cast<DOUBLE>(i); 
        }
    }

    result.minval = value;
    result.minidx = idx + 1;

    return result;
}

vector<int> min_v_k(const std::vector<DOUBLE>& v,int k) {

    priority_queue <double> topk;
    int len = v.size();
    for(int i = 0;i < k;i++)
    {
        topk.push(INFINITY);
    }

    for(int i = 0;i < len;i++)
    {
        if(v[i] < topk.top())
        {
            topk.pop();
            topk.push(v[i]);
        }
    }

    DOUBLE topK = topk.top();
    
    vector<int> result(k + 10);
    int index_result = 0;
    for(int i = 0;i < len;i++)
    {
        if(v[i] < topK)
        {
            result[index_result++] = i;
        }
    }

    return result;
}