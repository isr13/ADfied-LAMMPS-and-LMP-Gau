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

PairStyle(lj/charmmfsw/coul/charmmfsh,PairLJCharmmfswCoulCharmmfsh)

#else

#ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_CHARMMFSH_H
#define LMP_PAIR_LJ_CHARMMFSW_COUL_CHARMMFSH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCharmmfswCoulCharmmfsh : public Pair {
 public:
  PairLJCharmmfswCoulCharmmfsh(class LAMMPS *);
  virtual ~PairLJCharmmfswCoulCharmmfsh();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual myScalar init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual myScalar single(int, int, int, int, myScalar, myScalar, myScalar, myScalar &);
  virtual void *extract(const char *, int &);

 protected:
  int implicit;
  int dihedflag;

  myScalar cut_lj_inner,cut_lj,cut_coul,cut_coulinv,cut_ljinv,cut_lj_innerinv;
  myScalar cut_lj_innersq,cut_ljsq,cut_coulsq,cut_bothsq;
  myScalar cut_lj3inv,cut_lj_inner3inv,cut_lj3,cut_lj_inner3;
  myScalar cut_lj6inv,cut_lj_inner6inv,cut_lj6,cut_lj_inner6;
  myScalar denom_lj,denom_lj12,denom_lj6;

  myScalar **epsilon,**sigma,**eps14,**sigma14;
  myScalar **lj1,**lj2,**lj3,**lj4;
  myScalar **lj14_1,**lj14_2,**lj14_3,**lj14_4;

  virtual void allocate();
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

E: Pair style lj/charmmfsw/coul/charmmfsh requires atom attribute q

The atom style defined does not have these attributes.

E: Pair inner lj cutoff >= Pair outer lj cutoff

UNDOCUMENTED

U: Pair inner cutoff >= Pair outer cutoff

The specified cutoffs for the pair style are inconsistent.

*/
