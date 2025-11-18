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

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(qeq/reax,FixQEqReax)

#else

#ifndef LMP_FIX_QEQ_REAX_H
#define LMP_FIX_QEQ_REAX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixQEqReax : public Fix {
 public:
  FixQEqReax(class LAMMPS *, int, char **);
  ~FixQEqReax();
  int setmask();
  virtual void post_constructor();
  virtual void init();
  void init_list(int,class NeighList *);
  virtual void init_storage();
  void setup_pre_force(int);
  virtual void pre_force(int);

  void setup_pre_force_respa(int, int);
  void pre_force_respa(int, int, int);

  void min_setup_pre_force(int);
  void min_pre_force(int);

  int matvecs;
  myScalar qeq_time;

 protected:
  int nevery,reaxflag;
  int n, N, m_fill;
  int n_cap, nmax, m_cap;
  int pack_flag;
  int nlevels_respa;
  class NeighList *list;
  class PairReaxC *reaxc;

  myScalar swa, swb;      // lower/upper Taper cutoff radius
  myScalar Tap[8];        // Taper function
  myScalar tolerance;     // tolerance for the norm of the rel residual in CG

  myScalar *chi,*eta,*gamma;  // qeq parameters
  myScalar **shld;

  bigint ngroup;

  // fictitious charges

  myScalar *s, *t;
  myScalar **s_hist, **t_hist;
  int nprev;

  typedef struct{
    int n, m;
    int *firstnbr;
    int *numnbrs;
    int *jlist;
    myScalar *val;
  } sparse_matrix;

  sparse_matrix H;
  myScalar *Hdia_inv;
  myScalar *b_s, *b_t;
  myScalar *b_prc, *b_prm;

  //CG storage
  myScalar *p, *q, *r, *d;

  //GMRES storage
  //myScalar *g,*y;
  //myScalar **v;
  //myScalar **h;
  //myScalar *hc, *hs;

  char *pertype_option;  // argument to determine how per-type info is obtained
  virtual void pertype_parameters(char*);
  void init_shielding();
  void init_taper();
  virtual void allocate_storage();
  virtual void deallocate_storage();
  void reallocate_storage();
  virtual void allocate_matrix();
  void deallocate_matrix();
  void reallocate_matrix();

  virtual void init_matvec();
  void init_H();
  virtual void compute_H();
  myScalar calculate_H(myScalar,myScalar);
  virtual void calculate_Q();

  virtual int CG(myScalar*,myScalar*);
  //int GMRES(myScalar*,myScalar*);
  virtual void sparse_matvec(sparse_matrix*,myScalar*,myScalar*);

  virtual int pack_forward_comm(int, int *, myScalar *, int, int *);
  virtual void unpack_forward_comm(int, int, myScalar *);
  virtual int pack_reverse_comm(int, int, myScalar *);
  virtual void unpack_reverse_comm(int, int *, myScalar *);
  virtual myScalar memory_usage();
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  virtual int pack_exchange(int, myScalar *);
  virtual int unpack_exchange(int, myScalar *);

  virtual myScalar parallel_norm( myScalar*, int );
  virtual myScalar parallel_dot( myScalar*, myScalar*, int );
  virtual myScalar parallel_vector_acc( myScalar*, int );

  virtual void vector_sum(myScalar*,myScalar,myScalar*,myScalar,myScalar*,int);
  virtual void vector_add(myScalar*, myScalar, myScalar*,int);

  // dual CG support
  int dual_enabled;  // 0: Original, separate s & t optimization; 1: dual optimization
  int matvecs_s, matvecs_t; // Iteration count for each system
};

}

#endif
#endif
