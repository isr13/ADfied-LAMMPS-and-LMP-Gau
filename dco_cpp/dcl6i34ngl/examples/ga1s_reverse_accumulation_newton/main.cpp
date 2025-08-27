/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <iostream>
using namespace std;

#include "driver.hpp"

int main() {
  double x=100,p=5,xa1=1,pa1=0;
  driver(p,pa1,x,xa1);
  cout.precision(15);
  cout << "x=" << x << endl;
  cout << "p_{(1)}=" << pa1 << endl;
  return 0;
}

