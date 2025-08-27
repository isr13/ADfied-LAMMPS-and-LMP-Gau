/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

template<class DCO_TYPE>
void g(
    const DCO_TYPE& p, 
    DCO_TYPE& x
) {
  const DCO_TYPE eps=1e-12;
  DCO_TYPE xp=x+1;
  while (fabs(x-xp)>eps) {
    xp=x;
    x=xp-(xp*xp-p)/(2*xp);
  }
}

