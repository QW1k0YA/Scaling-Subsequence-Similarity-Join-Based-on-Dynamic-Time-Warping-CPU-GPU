
#ifndef PRO_MAX_ALLCLASS_H
#define PRO_MAX_ALLCLASS_H

#include "matrix.h"
#include "iostream"
#include "vector"
#include "algorithm"
#include "cmath"
#include "chrono"
#include "allstruct.h"
using namespace std;

class alpha_generater
{
public:
    NEXT_RETURN next(const vector<DOUBLE>& X,const vector<DOUBLE>& f_x,const vector<DOUBLE>& bsf_X,DOUBLE start_pos_,DOUBLE end_pos_);
private:
    double alpha_value;
    double diag;
};
#endif 
