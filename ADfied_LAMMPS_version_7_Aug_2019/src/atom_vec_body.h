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

AtomStyle(body,AtomVecBody)

#else

#ifndef LMP_ATOM_VEC_BODY_H
#define LMP_ATOM_VEC_BODY_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecBody : public AtomVec {
 public:
  class Body *bptr;

  struct Bonus {
    myScalar quat[4];
    myScalar inertia[3];
    int ninteger,nmyScalar;
    int iindex,dindex;
    int *ivalue;
    myScalar *dvalue;
    int ilocal;
  };
  struct Bonus *bonus;

  AtomVecBody(class LAMMPS *);
  ~AtomVecBody();
  void process_args(int, char **);
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, myScalar *, int, int *);
  int pack_comm_vel(int, int *, myScalar *, int, int *);
  int pack_comm_hybrid(int, int *, myScalar *);
  void unpack_comm(int, int, myScalar *);
  void unpack_comm_vel(int, int, myScalar *);
  int unpack_comm_hybrid(int, int, myScalar *);
  int pack_reverse(int, int, myScalar *);
  int pack_reverse_hybrid(int, int, myScalar *);
  void unpack_reverse(int, int *, myScalar *);
  int unpack_reverse_hybrid(int, int *, myScalar *);
  int pack_border(int, int *, myScalar *, int, int *);
  int pack_border_vel(int, int *, myScalar *, int, int *);
  int pack_border_hybrid(int, int *, myScalar *);
  void unpack_border(int, int, myScalar *);
  void unpack_border_vel(int, int, myScalar *);
  int unpack_border_hybrid(int, int, myScalar *);
  int pack_exchange(int, myScalar *);
  int unpack_exchange(myScalar *);
  int size_restart();
  int pack_restart(int, myScalar *);
  int unpack_restart(myScalar *);
  void create_atom(int, myScalar *);
  void data_atom(myScalar *, imageint, char **);
  int data_atom_hybrid(int, char **);
  void data_vel(int, char **);
  int data_vel_hybrid(int, char **);
  void pack_data(myScalar **);
  int pack_data_hybrid(int, myScalar *);
  void write_data(FILE *, int, myScalar **);
  int write_data_hybrid(FILE *, myScalar *);
  void pack_vel(myScalar **);
  int pack_vel_hybrid(int, myScalar *);
  void write_vel(FILE *, int, myScalar **);
  int write_vel_hybrid(FILE *, myScalar *);
  bigint memory_usage();

  // manipulate Bonus data structure for extra atom info

  void clear_bonus();
  void data_body(int, int, int, int *, myScalar *);

  // methods used by other classes to query/set body info

  myScalar radius_body(int, int, int *, myScalar *);
  void set_quat(int, myScalar *);

  int nlocal_bonus;

 private:
  tagint *tag;
  int *type,*mask;
  imageint *image;
  myScalar **x,**v,**f;
  myScalar *radius;
  myScalar *rmass;
  myScalar **angmom,**torque;
  int *body;

  int nghost_bonus,nmax_bonus;
  int intmyScalarratio;       // sizeof(myScalar) / sizeof(int)

  MyPoolChunk<int> *icp;
  MyPoolChunk<myScalar> *dcp;

  void grow_bonus();
  void copy_bonus(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Internal error in atom_style body

This error should not occur.  Contact the developers.

E: Invalid atom_style body command

No body style argument was provided.

E: Unrecognized body style

The choice of body style is unknown.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

E: Invalid density in Atoms section of data file

Density value cannot be <= 0.0.

E: Assigning body parameters to non-body atom

Self-explanatory.

E: Assigning quat to non-body atom

Self-explanatory.

*/
