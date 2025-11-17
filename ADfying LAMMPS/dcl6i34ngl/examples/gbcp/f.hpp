/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

template<typename T, typename U>
void f(const vector<T>& x, vector<U>& y) {
  T v=tan(x[2]*x[3]); T w=x[1]-v; 
  y[0]=x[0]*v/w; 
  y[1]=y[0]*x[1];
}


