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

#ifdef COMPUTE_CLASS

ComputeStyle(orientorder/atom,ComputeOrientOrderAtom)

#else

#ifndef LMP_COMPUTE_ORIENTORDER_ATOM_H
#define LMP_COMPUTE_ORIENTORDER_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOrientOrderAtom : public Compute {
 public:
  ComputeOrientOrderAtom(class LAMMPS *, int, char **);
  ~ComputeOrientOrderAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  myScalar memory_usage();
  myScalar cutsq;
  int iqlcomp, qlcomp, qlcompflag;
  int *qlist;
  int nqlist;

 private:
  int nmax,maxneigh,ncol,nnn;
  class NeighList *list;
  myScalar *distsq;
  int *nearest;
  myScalar **rlist;
  int qmax;
  myScalar **qnarray;
  myScalar **qnm_r;
  myScalar **qnm_i;

  void select3(int, int, myScalar *, int *, myScalar **);
  void calc_boop(myScalar **rlist, int numNeighbors,
                 myScalar qn[], int nlist[], int nnlist);
  myScalar dist(const myScalar r[]);

  myScalar polar_prefactor(int, int, myScalar);
  myScalar associated_legendre(int, int, myScalar);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute orientorder/atom requires a pair style be defined

Self-explanatory.

E: Compute orientorder/atom cutoff is longer than pairwise cutoff

Cannot compute order parameter beyond cutoff.

W: More than one compute orientorder/atom

It is not efficient to use compute orientorder/atom more than once.

*/
