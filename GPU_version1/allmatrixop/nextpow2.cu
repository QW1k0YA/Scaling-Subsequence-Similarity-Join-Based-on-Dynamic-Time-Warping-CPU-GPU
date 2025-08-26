
#include <iostream>
#include<cmath>
#include<algorithm>
using namespace std;

unsigned int nextpow2(unsigned int n) {

    int a = 0;
    while (pow(2, a) < n)
    {
        a++;
    }
    return a;
}