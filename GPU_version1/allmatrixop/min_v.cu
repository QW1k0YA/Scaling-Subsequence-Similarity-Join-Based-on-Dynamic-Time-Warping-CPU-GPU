
#include <queue>
#include "../alldef/matrix.cuh"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "../alldef/allstruct.cuh"
#include "../alldef/typedefdouble.cuh"
#include "../allunder/undermpx_v2.cuh"

using namespace std;

vector<int> min_v_k(const std::vector<FLOAT >& v, int k) {

    priority_queue <FLOAT> topk;
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

    FLOAT  topK = topk.top();
    
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