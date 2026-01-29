
#include <iostream>
#include <vector>
#include "../alldef/matrix.h"

std::vector<DOUBLE> movsum(const std::vector<DOUBLE> &ts, int b) {

    long long ts_size = ts.size();
    long long len = ts_size - b;

    DOUBLE sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }
    std::vector<DOUBLE> result(len);

    for(long long i = 0;i < len - 1;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
    result[len-1] = sum;
    return result;

}
void movsum(const std::vector<DOUBLE> &ts, int b,vector<DOUBLE> &result) {

    long long ts_size = ts.size();
    long long len = ts_size - b;
    DOUBLE sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }
    for(long long i = 0;i < len-1;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
    result[len-1] = sum;
}
std::vector<DOUBLE> movsum_p(double *ts,int ts_size, int b) {

    long long len = ts_size - b;

    DOUBLE sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }
    std::vector<DOUBLE> result(len);

    for(long long i = 0;i < len-1;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
    result[len-1] = sum;
    return result;

}

void movsum_p(double *ts,int ts_size, int b,vector<double> &result) {

    long long len = ts_size - b;

    DOUBLE sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }

    for(long long i = 0;i < len-1;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
    result[len-1] = sum;

}

void movsum_p_from_pos(double *ts,int ts_size, int b,vector<double> &result,int pos) {

    long long len = ts_size - b;

    DOUBLE sum = 0;
    for(int i = pos;i <= b + pos;i++)
    {
        sum += ts[i];
    }

    for(long long i = pos;i < len-1;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
    result[len-1] = sum;

}