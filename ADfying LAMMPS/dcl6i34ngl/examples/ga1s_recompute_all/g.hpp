/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

template<typename AD_TYPE>
void g(
    int n,   
    AD_TYPE& x
) {
  for (int i=0;i<n;i++) x=sin(x);
}

