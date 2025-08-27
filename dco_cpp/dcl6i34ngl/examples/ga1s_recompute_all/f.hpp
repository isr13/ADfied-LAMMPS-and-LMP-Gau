/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include "g_gap.hpp"

template<typename AD_MODE>
void f(
    int n,     
    typename AD_MODE::type& x
) {
  x=sin(x);
  // g(n-2,x); // for split reversal
  g_make_gap<AD_MODE>(n-2,x);
  x=sin(x);
}


