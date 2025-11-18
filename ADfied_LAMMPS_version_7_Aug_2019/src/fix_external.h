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

FixStyle(external,FixExternal)

#else

#ifndef LMP_FIX_EXTERNAL_H
#define LMP_FIX_EXTERNAL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixExternal : public Fix {
 public:
  myScalar **fexternal;

  FixExternal(class LAMMPS *, int, char **);
  ~FixExternal();
  int setmask();
  void init();
  void setup(int);
  void setup_pre_reverse(int, int);
  void min_setup(int);
  void pre_reverse(int, int);
  void post_force(int);
  void min_post_force(int);
  myScalar compute_scalar();
  myScalar compute_vector(int);

  void set_energy_global(myScalar);
  void set_virial_global(myScalar *);
  void set_energy_peratom(myScalar *);
  void set_virial_peratom(myScalar **);
  void set_vector_length(int);
  void set_vector(int,myScalar);

  myScalar memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, myScalar *);
  int unpack_exchange(int, myScalar *);

  typedef void (*FnPtr)(void *, bigint, int, tagint *, myScalar **, myScalar **);
  void set_callback(FnPtr, void *);

 private:
  int mode,ncall,napply,eflag_caller;
  FnPtr callback;
  void *ptr_caller;
  myScalar user_energy;
  myScalar *caller_vector;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix external callback function not set

This must be done by an external program in order to use this fix.

E: Invalid set_vector index in fix external

UNDOCUMENTED

*/
