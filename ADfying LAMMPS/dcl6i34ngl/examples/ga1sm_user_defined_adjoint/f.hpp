/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include "g_gap.hpp"

template<typename AD_MODE, class ATYPE>
void f(ATYPE p, ATYPE& x) {
  p=exp(p);
  g_make_gap<AD_MODE>(p,x);
  x=x*p;
}


