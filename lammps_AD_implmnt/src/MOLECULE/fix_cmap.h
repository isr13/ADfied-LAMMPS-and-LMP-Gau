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

FixStyle(cmap,FixCMAP)

#else

#ifndef LMP_FIX_CMAP_H
#define LMP_FIX_CMAP_H

#include "fix.h"
namespace LAMMPS_NS {

class FixCMAP : public Fix {
 public:
  FixCMAP(class LAMMPS *, int, char **);
  ~FixCMAP();
  int setmask();
  void init();
  void setup(int);
  void setup_pre_neighbor();
  void setup_pre_reverse(int, int);
  void min_setup(int);
  void pre_neighbor();
  void pre_reverse(int, int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  myScalar compute_scalar();

  void read_data_header(char *);
  void read_data_section(char *, int, char *, tagint);
  bigint read_data_skip_lines(char *);
  void write_data_header(FILE *, int);
  void write_data_section_size(int, int &, int &);
  void write_data_section_pack(int, myScalar **);
  void write_data_section_keyword(int, FILE *);
  void write_data_section(int, FILE *, int, myScalar **, int);

  void write_restart(FILE *);
  void restart(char *);
  int pack_restart(int, myScalar *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, myScalar *);
  int unpack_exchange(int, myScalar *);

  myScalar memory_usage();

 private:
  int nprocs,me;
  int newton_bond,eflag_caller;
  int ctype,nlevels_respa;
  int ncrosstermtypes,crossterm_per_atom,maxcrossterm;
  int ncrosstermlist;
  bigint ncmap;

  int **crosstermlist;

  int nmax_previous;
  int *num_crossterm;
  int **crossterm_type;
  tagint **crossterm_atom1,**crossterm_atom2,**crossterm_atom3;
  tagint **crossterm_atom4,**crossterm_atom5;

  myScalar E,dEdPhi,dEdPsi;
  myScalar ecmap;
  myScalar fcmap[4],cij[4][4];
  myScalar *g_axis;

  // CMAP grid points obtained from external file

  myScalar ***cmapgrid;

  // partial derivatives and cross-derivatives of the grid data

  myScalar ***d1cmapgrid,***d2cmapgrid,***d12cmapgrid;

  // read map grid data

  void read_grid_map(char *);

  // read in CMAP cross terms from LAMMPS data file

  void read_cmap_data(int, char *);

  // pre-compute the partial and cross-derivatives of map grid points

  void set_map_derivatives(myScalar **, myScalar **, myScalar **, myScalar **);

  // cubic spline interpolation functions for derivatives of map grid points

  void spline(myScalar *, myScalar *, int);
  void spl_interpolate(myScalar, myScalar *, myScalar *, myScalar &, myScalar &);

  // calculate dihedral angles

  myScalar dihedral_angle_atan2(myScalar, myScalar, myScalar, myScalar, myScalar, myScalar,
                              myScalar, myScalar, myScalar, myScalar);

  // calculate bicubic interpolation coefficient matrix c_ij

  void bc_coeff(myScalar *, myScalar *, myScalar *, myScalar *);

  // perform bicubic interpolation at point of interest

  void bc_interpol(myScalar, myScalar, int, int, myScalar *, myScalar *, myScalar *,
                   myScalar *);
};
}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: CMAP atoms %d %d %d %d %d missing on proc %d at step %ld

UNDOCUMENTED

E: Invalid CMAP crossterm_type

UNDOCUMENTED

E: Cannot open fix cmap file %s

UNDOCUMENTED

E: CMAP: atan2 function cannot take 2 zero arguments

UNDOCUMENTED

E: Invalid read data header line for fix cmap

UNDOCUMENTED

E: Incorrect %s format in data file

UNDOCUMENTED

E: Too many CMAP crossterms for one atom

UNDOCUMENTED

*/
