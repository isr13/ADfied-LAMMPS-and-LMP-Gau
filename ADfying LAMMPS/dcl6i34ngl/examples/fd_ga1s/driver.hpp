/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

void fd_driver(
    const vector<double>& xv, 
    vector<double>& xa,     
    vector<double>& xa1t2x,     
    vector<double>& yv,       
    vector<double>& yt2x,       
    vector<double>& ya      
);

