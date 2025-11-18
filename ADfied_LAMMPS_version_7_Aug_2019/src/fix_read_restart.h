#include "ad_defines.h"
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(READ_RESTART,FixReadRestart)

#else

#ifndef LMP_FIX_READ_RESTART_H
#define LMP_FIX_READ_RESTART_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReadRestart : public Fix {
 public:
  int *count;
  myScalar **extra;

  FixReadRestart(class LAMMPS *, int, char **);
  ~FixReadRestart();
  int setmask();

  myScalar memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, myScalar *);
  int unpack_exchange(int, myScalar *);

 private:
  int nextra;          // max number of extra values for any atom
};

}

#endif
#endif
/* ERROR/WARNING messages:

*/
