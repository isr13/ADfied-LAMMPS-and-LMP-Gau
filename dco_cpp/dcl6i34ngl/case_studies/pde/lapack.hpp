/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#ifndef __LAPACK_INCLUDED_
#define __LAPACK_INCLUDED_

template<typename ATYPE>
void dgttrf(int n, ATYPE *dl, ATYPE *d, ATYPE *du, ATYPE *du2, int *ipiv, int *info);

template<typename ATYPE>
void dgttrs(const char *trans, int n, int nrhs,
            ATYPE *dl, ATYPE *d, ATYPE *du, ATYPE *du2,
            int *ipiv, ATYPE *b, int ldb, int *info);

#include "lapack.cpp"

#endif

