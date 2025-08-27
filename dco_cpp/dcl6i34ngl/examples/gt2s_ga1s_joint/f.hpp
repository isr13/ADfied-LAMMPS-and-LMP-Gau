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
void f(int n, DCO_TYPE& x) {
  g(n/3,x);
  g_make_gap(n/3,x);
  g(n-n/3*2,x);
}

