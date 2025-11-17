/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include "g.hpp"

template<typename DCO_TYPE>
void f(DCO_TYPE p, DCO_TYPE& x) {
  p=exp(p);
  g(p,x);
  x=x*p;
}


