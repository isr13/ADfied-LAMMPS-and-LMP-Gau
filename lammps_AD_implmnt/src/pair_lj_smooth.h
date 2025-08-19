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

#ifdef PAIR_CLASS

PairStyle(lj/smooth,PairLJSmooth)

#else

#ifndef LMP_PAIR_LJ_SMOOTH_H
#define LMP_PAIR_LJ_SMOOTH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSmooth : public Pair {
 public:
  PairLJSmooth(class LAMMPS *);
  virtual ~PairLJSmooth();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  myScalar init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  myScalar single(int, int, int, int, myScalar, myScalar, myScalar, myScalar &);

 protected:
  myScalar cut_inner_global,cut_global;
  myScalar **cut,**cut_inner,**cut_inner_sq;
  myScalar **epsilon,**sigma;
  myScalar **lj1,**lj2,**lj3,**lj4;
  myScalar **ljsw0,**ljsw1,**ljsw2,**ljsw3,**ljsw4;
  myScalar **offset;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
