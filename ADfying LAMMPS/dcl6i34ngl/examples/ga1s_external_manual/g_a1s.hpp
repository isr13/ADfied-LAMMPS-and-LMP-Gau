/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include<vector>

template<typename DCO_TYPE>
void g_a1s(std::vector<DCO_TYPE>& xa1, const DCO_TYPE& ya1) {
  typename std::vector<DCO_TYPE>::iterator i;
  for (i=xa1.begin();i!=xa1.end();i++) *i+=ya1;
}

