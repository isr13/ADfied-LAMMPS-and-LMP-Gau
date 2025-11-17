/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include "dco.hpp"
#include "f.hpp"

void driver (double& x, double& xt1) {
  dco::gt1s<double>::type t1s_x=x;
  dco::derivative(t1s_x) = xt1;
  f(t1s_x); 
  x = dco::value(t1s_x);
  xt1 = dco::derivative(t1s_x);
}

