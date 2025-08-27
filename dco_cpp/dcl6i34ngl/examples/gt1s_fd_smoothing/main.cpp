/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <cstdlib>
#include <iostream>
using namespace std;

#include "driver.hpp"

int main() {
    double x = 0.1,xt1=1;
    driver(x,xt1);
    cout.precision(15);
    cout << "x=" << x << endl;
    cout << "x^{(1)}=" << xt1 << endl;
    return 0;
}

