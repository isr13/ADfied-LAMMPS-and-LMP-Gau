/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

#include <iostream>
#include <vector>
using namespace std;

#include "driver.hpp"

int main() {
  const int n=4, m=5; cout.precision(15);
  vector<double> xv(n), yv(m);
  vector<vector<double> > xa(n, vector<double>(m)), ya(m, vector<double>(m));
  for (int i=0;i<n;i++) { 
    xv[i]=1; 
    for (int j=0;j<m;j++) xa[i][j]=0; 
  }
  for (int i=0;i<m;i++) { 
    for (int j=0;j<m;j++) ya[i][j]=0; 
    ya[i][i]=1; 
  }
  driver(xv,xa,yv,ya);
  for (int i=0;i<m;i++)
    cout << "y[" << i << "]=" << yv[i] << endl;
  for (int i=0;i<n;i++) 
    for (int j=0;j<m;j++) 
      cout << "x_{(1)}[" << i << "][" << j << "]=" << xa[i][j] << endl;
  return 0;
}



