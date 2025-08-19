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

#ifndef LMP_ATOM_VEC_H
#define LMP_ATOM_VEC_H

#include "pointers.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class AtomVec : protected Pointers {
 public:
  int molecular;                       // 0 = atomic, 1 = molecular system
  int bonds_allow,angles_allow;        // 1 if bonds, angles are used
  int dihedrals_allow,impropers_allow; // 1 if dihedrals, impropers used
  int mass_type;                       // 1 if per-type masses
  int dipole_type;                     // 1 if per-type dipole moments
  int forceclearflag;                  // 1 if has forceclear() method

  int comm_x_only;                     // 1 if only exchange x in forward comm
  int comm_f_only;                     // 1 if only exchange f in reverse comm

  int size_forward;                    // # of values per atom in comm
  int size_reverse;                    // # in reverse comm
  int size_border;                     // # in border comm
  int size_velocity;                   // # of velocity based quantities
  int size_data_atom;                  // number of values in Atom line
  int size_data_vel;                   // number of values in Velocity line
  int size_data_bonus;                 // number of values in Bonus line
  int xcol_data;                       // column (1-N) where x is in Atom line
  int maxexchange;                     // max size of exchanged atom
                                       // only needs to be set if size > BUFEXTRA

  class Molecule **onemols;            // list of molecules for style template
  int nset;                            // # of molecules in list

  int kokkosable;                      // 1 if atom style is KOKKOS-enabled

  int nargcopy;          // copy of command-line args for atom_style command
  char **argcopy;        // used when AtomVec is realloced (restart,replicate)

  AtomVec(class LAMMPS *);
  virtual ~AtomVec();
  void store_args(int, char **);
  virtual void process_args(int, char **);
  virtual void init();

  virtual void grow(int) = 0;
  virtual void grow_reset() = 0;
  virtual void copy(int, int, int) = 0;
  virtual void clear_bonus() {}
  virtual void force_clear(int, size_t) {}

  virtual int pack_comm(int, int *, myScalar *, int, int *) = 0;
  virtual int pack_comm_vel(int, int *, myScalar *, int, int *) = 0;
  virtual int pack_comm_hybrid(int, int *, myScalar *) {return 0;}
  virtual void unpack_comm(int, int, myScalar *) = 0;
  virtual void unpack_comm_vel(int, int, myScalar *) = 0;
  virtual int unpack_comm_hybrid(int, int, myScalar *) {return 0;}

  virtual int pack_reverse(int, int, myScalar *) = 0;
  virtual int pack_reverse_hybrid(int, int, myScalar *) {return 0;}
  virtual void unpack_reverse(int, int *, myScalar *) = 0;
  virtual int unpack_reverse_hybrid(int, int *, myScalar *) {return 0;}

  virtual int pack_border(int, int *, myScalar *, int, int *) = 0;
  virtual int pack_border_vel(int, int *, myScalar *, int, int *) = 0;
  virtual int pack_border_hybrid(int, int *, myScalar *) {return 0;}
  virtual void unpack_border(int, int, myScalar *) = 0;
  virtual void unpack_border_vel(int, int, myScalar *) = 0;
  virtual int unpack_border_hybrid(int, int, myScalar *) {return 0;}

  virtual int pack_exchange(int, myScalar *) = 0;
  virtual int unpack_exchange(myScalar *) = 0;

  virtual int size_restart() = 0;
  virtual int pack_restart(int, myScalar *) = 0;
  virtual int unpack_restart(myScalar *) = 0;

  virtual void create_atom(int, myScalar *) = 0;

  virtual void data_atom(myScalar *, imageint, char **) = 0;
  virtual void data_atom_bonus(int, char **) {}
  virtual int data_atom_hybrid(int, char **) {return 0;}
  virtual void data_vel(int, char **);
  virtual int data_vel_hybrid(int, char **) {return 0;}

  virtual void pack_data(myScalar **) = 0;
  virtual int pack_data_hybrid(int, myScalar *) {return 0;}
  virtual void write_data(FILE *, int, myScalar **) = 0;
  virtual int write_data_hybrid(FILE *, myScalar *) {return 0;}
  virtual void pack_vel(myScalar **);
  virtual int pack_vel_hybrid(int, myScalar *) {return 0;}
  virtual void write_vel(FILE *, int, myScalar **);
  virtual int write_vel_hybrid(FILE *, myScalar *) {return 0;}

  int pack_bond(tagint **);
  void write_bond(FILE *, int, tagint **, int);
  int pack_angle(tagint **);
  void write_angle(FILE *, int, tagint **, int);
  int pack_dihedral(tagint **);
  void write_dihedral(FILE *, int, tagint **, int);
  int pack_improper(tagint **);
  void write_improper(FILE *, int, tagint **, int);

  virtual int property_atom(char *) {return -1;}
  virtual void pack_property_atom(int, myScalar *, int, int) {}

  virtual bigint memory_usage() = 0;

 protected:
  int nmax;                             // local copy of atom->nmax
  int deform_vremap;                    // local copy of domain properties
  int deform_groupbit;
  myScalar *h_rate;

  // union data struct for packing 32-bit and 64-bit ints into myScalar bufs
  // this avoids aliasing issues by having 2 pointers (myScalar,int)
  //   to same buf memory
  // constructor for 32-bit int prevents compiler
  //   from possibly calling the myScalar constructor when passed an int
  // copy to a myScalar *buf:
  //   buf[m++] = ubuf(foo).d, where foo is a 32-bit or 64-bit int
  // copy from a myScalar *buf:
  //   foo = (int) ubuf(buf[m++]).i;, where (int) or (tagint) match foo
  //   the cast prevents compiler warnings about possible truncation

  union ubuf {
    myScalar d;
    int64_t i;
    ubuf(myScalar arg) : d(arg) {}
    ubuf(int64_t arg) : i(arg) {}
    ubuf(int arg) : i(arg) {}
  };

  void grow_nmax();
  int grow_nmax_bonus(int);
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid atom_style command

Self-explanatory.

E: KOKKOS package requires a kokkos enabled atom_style

Self-explanatory.

*/
