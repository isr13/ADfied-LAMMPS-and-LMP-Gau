/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include<iostream>
#include<vector>
using namespace std;

#include "driver.hpp" 

int main() {
  const int n=4; cout.precision(15);

  {
    vector<double> xv(n); double yv;
    vector<vector<double> > hessian(n, vector<double>(n));

    for (int i=0;i<n;i++) { xv[i]=1; }

    t2s_a1s_Hessian(xv,yv,hessian);

    cout << "y=" << yv << endl;

    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) {
        cout << "hessian[" << i << "][" << j << "]=" << hessian[i][j] << endl;
      }
    }
  }

  {
    vector<double> xv(n); double yv;
    vector<vector<double> > hessian(n, vector<double>(n));

    for (int i=0;i<n;i++) { xv[i]=1; }

    a2s_a1s_Hessian(xv,yv,hessian);

    cout << "y=" << yv << endl;

    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) {
        cout << "hessian[" << i << "][" << j << "]=" << hessian[i][j] << endl;
      }
    }
  }
  return 0;
}

