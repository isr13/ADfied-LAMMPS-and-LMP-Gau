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
void f(int n, int m, DCO_TYPE& x) {
  for (int i=0;i<n;i+=m)
    g_make_gap(std::min(m,n-i),x);
}


