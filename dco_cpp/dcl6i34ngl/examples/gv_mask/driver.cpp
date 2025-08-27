/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <vector> 
using namespace std;

#include "dco.hpp" 
using namespace dco;

#ifdef DCO_CPP11
#include "f.hpp"
typedef double DCO_BASE_TYPE;
static const unsigned long VECTOR_SIZE = 8;
typedef gv<DCO_BASE_TYPE,VECTOR_SIZE> DCO_MODE;
typedef DCO_MODE::type DCO_TYPE;

void driver(
            const vector<double>& x,
            vector<double>& y
) {
  const size_t n=x.size();
  DCO_TYPE xv, yv;
  for (auto range : dco::subranges<VECTOR_SIZE>(n)) {
    for (auto i : range) { xv[i.sub_index]=x[i.global_index]; }
    f(xv, yv); 
    for (auto i : range) { y[i.global_index]=yv[i.sub_index]; }
  }
}
#else
#include "f_branch.hpp"
void driver(
            const vector<double>& x,
            vector<double>& y
) {
  for (size_t i=0;i<x.size();i++) { f(x[i],y[i]); }
}
#endif

