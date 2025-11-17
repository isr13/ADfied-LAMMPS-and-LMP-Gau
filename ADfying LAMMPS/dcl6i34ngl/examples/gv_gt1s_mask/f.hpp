/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

template<typename T>
void f(const T& x, T& y) {
  T v=sin(x*x);
  T w=-cos(v); 
  y=v*w;
#ifdef DCO_CPP11
  y += dco::vmask(y < -0.5, sin(exp(x)));
#else
  if (y < -0.5) y += sin(exp(x));
#endif
}

