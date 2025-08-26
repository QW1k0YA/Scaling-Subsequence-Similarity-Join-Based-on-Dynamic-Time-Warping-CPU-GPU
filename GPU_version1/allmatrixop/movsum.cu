
#include <iostream>
#include <vector>
#include "../alldef/matrix.cuh"

std::vector<FLOAT > movsum(const std::vector<FLOAT > &ts, int b) {

    long long ts_size = ts.size();
    long long len = ts_size - b;

    FLOAT  sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }
    std::vector<FLOAT > result(len);

    for(long long i = 0;i < len - 1;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
    result[len - 1] = sum;

    return result;

}
void movsum(const std::vector<FLOAT > &ts, int b, vector<FLOAT > &result) {

    long long ts_size = ts.size();
    long long len = ts_size - b;
    FLOAT  sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }
    for(long long i = 0;i < len;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
}
std::vector<FLOAT > movsum_p(FLOAT *ts, int ts_size, int b) {

    long long len = ts_size - b;

    FLOAT  sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }
    std::vector<FLOAT > result(len);

    for(long long i = 0;i < len;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }

    return result;

}

void movsum_p(FLOAT *ts,int ts_size, int b,vector<FLOAT> &result) {

    long long len = ts_size - b;

    FLOAT  sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }

    for(long long i = 0;i < len;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }

}

void movsum_p(const FLOAT *ts,int ts_size, int b,FLOAT *result) {

    long long len = ts_size - b;

    FLOAT  sum = 0;
    for(int i = 0;i <= b;i++)
    {
        sum += ts[i];
    }

    for(long long i = 0;i < len;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }
}
void movsum_p_from_pos(FLOAT *ts,int ts_size, int b,vector<FLOAT> &result,int pos) {

    long long len = ts_size - b;

    FLOAT  sum = 0;
    for(int i = pos;i <= b + pos;i++)
    {
        sum += ts[i];
    }

    for(long long i = pos;i < len;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }

}

void movsum_p_from_pos(const FLOAT *ts,int ts_size, int b,FLOAT *result,int pos) {

    long long len = ts_size - b;

    FLOAT  sum = 0;
    for(int i = pos;i <= b + pos;i++)
    {
        sum += ts[i];
    }

    for(long long i = pos;i < len;i++)
    {
        result[i] = sum;
        sum = sum - ts[i] + ts[i+b+1];
    }

}