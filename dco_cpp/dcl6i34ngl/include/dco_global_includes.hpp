/*
dco/c++ v3.4.1
  - Algorithmic Differentiation by Operator Overloading in C++

COPYRIGHT
The Numerical Algorithms Group Limited and
Software and Tools for Computational Engineering @ RWTH Aachen University, 2017

This file is part of dco/c++.


*/

// This file includes all std headers required by dco/c++.
#pragma once

#include <sstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <complex>
#include <string>
#include <stack>
#include <exception>
#include <stdexcept>
#include <bitset>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <list>
#include <algorithm>
#include <sys/stat.h>
#include <cassert>
#include <limits>
#include <iomanip>
#include <cctype>
#include <fcntl.h>
#include <typeinfo>


#ifdef DCO_CPP11
#  include <type_traits>
#  include <memory>
#endif

#if !defined(_WIN32)
#  include <cxxabi.h>
#  include <unistd.h>
#  include <sys/mman.h>
#  include <sys/time.h>
#endif

#ifdef _WIN32
#  include <time.h>
#  ifndef DCO_SKIP_WINDOWS_H_INCLUDE
#    ifndef NOMINMAX
#      define NOMINMAX
#    endif
#    include <windows.h>
#  endif
#endif

#include "dco_version.hpp"


