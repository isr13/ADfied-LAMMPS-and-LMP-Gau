/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include "g_gap.hpp"

template<typename DCO_TYPE>
void f(int n, int m, const DCO_TYPE& x, DCO_TYPE& y) {
  double r;
  DCO_TYPE sum;
  for (int i=0;i<n;i++) {
    r=rand();
    g_make_gap(m,x,r,y);
    sum+=y;
  }
  y=sum/n;
}

