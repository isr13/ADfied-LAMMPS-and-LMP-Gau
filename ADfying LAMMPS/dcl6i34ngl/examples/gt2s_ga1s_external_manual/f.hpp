/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include<vector>
#include "g_gap.hpp"

template<typename DCO_TYPE>
void f(std::vector<DCO_TYPE>& x, DCO_TYPE& y) {
  for (unsigned i=0;i<x.size();i++) x[i]*=x[i];
  g(x,y);
  // g_make_gap(x,y);
  y*=y;
}

