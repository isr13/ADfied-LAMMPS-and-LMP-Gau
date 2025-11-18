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

PairStyle(zbl,PairZBL)

#else

#ifndef LMP_PAIR_ZBL_H
#define LMP_PAIR_ZBL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairZBL : public Pair {
 public:
  PairZBL(class LAMMPS *);
  virtual ~PairZBL();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual myScalar init_one(int, int);
  myScalar single(int, int, int, int, myScalar, myScalar, myScalar, myScalar &);

 protected:
  myScalar cut_global, cut_inner;
  myScalar cut_globalsq, cut_innersq;
  myScalar *z;
  myScalar **d1a,**d2a,**d3a,**d4a,**zze;
  myScalar **sw1,**sw2,**sw3,**sw4,**sw5;

  virtual void allocate();
  myScalar e_zbl(myScalar, int, int);
  myScalar dzbldr(myScalar, int, int);
  myScalar d2zbldr2(myScalar, int, int);
  void set_coeff(int, int, myScalar, myScalar);
};

namespace PairZBLConstants {

  // ZBL constants

  static const myScalar pzbl = 0.23;
  static const myScalar a0 = 0.46850;
  static const myScalar c1 = 0.02817;
  static const myScalar c2 = 0.28022;
  static const myScalar c3 = 0.50986;
  static const myScalar c4 = 0.18175;
  static const myScalar d1 = 0.20162;
  static const myScalar d2 = 0.40290;
  static const myScalar d3 = 0.94229;
  static const myScalar d4 = 3.19980;
}

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
