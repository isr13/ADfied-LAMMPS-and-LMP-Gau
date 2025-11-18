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

#ifdef ATOM_CLASS

AtomStyle(atomic,AtomVecAtomic)

#else

#ifndef LMP_ATOM_VEC_ATOMIC_H
#define LMP_ATOM_VEC_ATOMIC_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecAtomic : public AtomVec {
 public:
  AtomVecAtomic(class LAMMPS *);
  virtual ~AtomVecAtomic() {}
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  virtual int pack_comm(int, int *, myScalar *, int, int *);
  virtual int pack_comm_vel(int, int *, myScalar *, int, int *);
  virtual void unpack_comm(int, int, myScalar *);
  virtual void unpack_comm_vel(int, int, myScalar *);
  int pack_reverse(int, int, myScalar *);
  void unpack_reverse(int, int *, myScalar *);
  virtual int pack_border(int, int *, myScalar *, int, int *);
  virtual int pack_border_vel(int, int *, myScalar *, int, int *);
  virtual void unpack_border(int, int, myScalar *);
  virtual void unpack_border_vel(int, int, myScalar *);
  virtual int pack_exchange(int, myScalar *);
  virtual int unpack_exchange(myScalar *);
  int size_restart();
  int pack_restart(int, myScalar *);
  int unpack_restart(myScalar *);
  void create_atom(int, myScalar *);
  void data_atom(myScalar *, imageint, char **);
  void pack_data(myScalar **);
  void write_data(FILE *, int, myScalar **);
  bigint memory_usage();

 protected:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  myScalar **x,**v,**f;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

*/
