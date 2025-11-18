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

#ifndef LMP_MIN_LSRCH_H
#define LMP_MIN_LSRCH_H

#include "min.h"

namespace LAMMPS_NS {

class MinLineSearch : public Min {
 public:
  MinLineSearch(class LAMMPS *);
  ~MinLineSearch();
  void init();
  void setup_style();
  void reset_vectors();

 protected:
  // vectors needed by linesearch minimizers
  // allocated and stored by fix_minimize
  // x,f are stored by parent or Atom class or Pair class

  myScalar *x0;                 // coords at start of linesearch
  myScalar *g;                  // old gradient vector
  myScalar *h;                  // search direction vector

  myScalar *gextra;             // g,h for extra global dof, x0 is stored by fix
  myScalar *hextra;

  myScalar **x0extra_atom;      // x0,g,h for extra per-atom dof
  myScalar **gextra_atom;
  myScalar **hextra_atom;

  typedef int (MinLineSearch::*FnPtr)(myScalar, myScalar &);
  FnPtr linemin;
  int linemin_backtrack(myScalar, myScalar &);
  int linemin_quadratic(myScalar, myScalar &);
  int linemin_forcezero(myScalar, myScalar &);

  myScalar alpha_step(myScalar, int);
  myScalar compute_dir_deriv(myScalar &);
};

}

#endif
