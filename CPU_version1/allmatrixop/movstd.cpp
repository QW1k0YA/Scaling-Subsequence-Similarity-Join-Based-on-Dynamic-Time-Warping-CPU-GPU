
#include <iostream>
#include <vector>
#include <cmath>
#include "../alldef/matrix.h"

std::vector<DOUBLE> movstd(const std::vector<DOUBLE>& ts, int a, int b,bool c) {
    size_t ts_size = ts.size();
    if (c) {
        std::vector<DOUBLE> result;

        if (ts.empty() || a < 0 || b < 0 || a + b == 0) {
            std::cerr << "Invalid input parameters. in movstd" << std::endl;
            return result;
        }

        for (size_t i = 0; i < ts_size; ++i) {
            DOUBLE sum = 0.0;
            DOUBLE sumSquared = 0.0;
            int count = 0;

            for (int j = -a; j <= b; ++j) {
                int index = i + j;
                if (index >= 0 && index < ts_size) {
                    DOUBLE value = ts[index];
                    sum += value;
                    sumSquared += value * value;
                    count++;
                }
            }

            if (count == a + b + 1) {
                
                DOUBLE mean = sum / count;
                DOUBLE variance = (sumSquared - sum * mean) / count;
                DOUBLE stdDev = std::sqrt(variance);
                result.push_back(stdDev);
            }
        }

        return result;
    } else {
        std::vector<DOUBLE> result;

        if (ts.empty() || a < 0 || b < 0 || a + b == 0) {
            std::cerr << "Invalid input parameters.in movstd" << std::endl;
            return result;
        }

        for (size_t i = 0; i < ts_size; ++i) {
            DOUBLE sum = 0.0;
            DOUBLE sumSquared = 0.0;
            int count = 0;

            for (int j = -a; j <= b; ++j) {
                if (i + j >= 0 && i + j < ts_size) {
                    DOUBLE value = ts[i + j];
                    sum += value;
                    sumSquared += value * value;
                    count++;
                }
            }

            if (count > 0) {
                DOUBLE mean = sum / count;
                DOUBLE variance = (sumSquared - count * mean * mean) / count;
                DOUBLE stdDev = std::sqrt(variance);
                result.push_back(stdDev);
            } else {
                
                result.push_back(ts[i]);
            }
        }

        return result;
    }
}