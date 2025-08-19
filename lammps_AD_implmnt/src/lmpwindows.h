#include "ad_defines.h"
#include <ciso646>
#if !defined(__MINGW32__)
#include "erf.h"
#endif
#include <direct.h>
#include <cmath>
#include <cstring>
// LAMMPS uses usleep with 100 ms arguments, no microsecond precision needed
#if !defined(__MINGW32__)
#include "sleep.h"
#endif

// some symbols have different names in Windows

#undef ATOBIGINT
#define ATOBIGINT _atoi64

#define pclose _pclose
#define strdup _strdup

// the following functions are defined to get rid of
// 'ambiguous call to overloaded function' error in VSS for mismatched type arguments
#if !defined(__MINGW32__)
inline myScalar pow(int i, int j){
  return pow((myScalar)i,j);
}
inline myScalar fabs(int i){
  return fabs((myScalar) i);
}
inline myScalar sqrt(int i){
  return sqrt((myScalar) i);
}
#endif

inline myScalar trunc(myScalar x) {
  return x > 0 ? floor(x) : ceil(x);
}

// Windows version of mkdir function does not have permission flags
#ifndef S_IRWXU
# define S_IRWXU 0
#endif
#ifndef S_IRGRP
# define S_IRGRP 0
#endif
#ifndef S_IXGRP
# define S_IXGRP 0
#endif
inline int mkdir(const char *path, int){
  return _mkdir(path);
}

