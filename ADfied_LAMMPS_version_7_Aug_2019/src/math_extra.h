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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_H
#define LMP_MATH_EXTRA_H

#include <cmath>

namespace MathExtra {

  // 3 vector operations

  inline void copy3(const myScalar *v, myScalar *ans);
  inline void zero3(myScalar *v);
  inline void norm3(myScalar *v);
  inline void normalize3(const myScalar *v, myScalar *ans);
  inline void snormalize3(const myScalar, const myScalar *v, myScalar *ans);
  inline void negate3(myScalar *v);
  inline void scale3(myScalar s, myScalar *v);
  inline void add3(const myScalar *v1, const myScalar *v2, myScalar *ans);
  inline void scaleadd3(myScalar s, const myScalar *v1, const myScalar *v2,
                        myScalar *ans);
  inline void sub3(const myScalar *v1, const myScalar *v2, myScalar *ans);
  inline myScalar len3(const myScalar *v);
  inline myScalar lensq3(const myScalar *v);
  inline myScalar distsq3(const myScalar *v1, const myScalar *v2);
  inline myScalar dot3(const myScalar *v1, const myScalar *v2);
  inline void cross3(const myScalar *v1, const myScalar *v2, myScalar *ans);

  // 3x3 matrix operations

  inline void col2mat(const myScalar *ex, const myScalar *ey, const myScalar *ez,
                      myScalar m[3][3]);
  inline myScalar det3(const myScalar mat[3][3]);
  inline void diag_times3(const myScalar *d, const myScalar m[3][3],
                          myScalar ans[3][3]);
  inline void times3_diag(const myScalar m[3][3], const myScalar *d,
                          myScalar ans[3][3]);
  inline void plus3(const myScalar m[3][3], const myScalar m2[3][3],
                    myScalar ans[3][3]);
  inline void times3(const myScalar m[3][3], const myScalar m2[3][3],
                     myScalar ans[3][3]);
  inline void transpose_times3(const myScalar m[3][3], const myScalar m2[3][3],
                               myScalar ans[3][3]);
  inline void times3_transpose(const myScalar m[3][3], const myScalar m2[3][3],
                               myScalar ans[3][3]);
  inline void invert3(const myScalar mat[3][3], myScalar ans[3][3]);
  inline void matvec(const myScalar mat[3][3], const myScalar *vec, myScalar *ans);
  inline void matvec(const myScalar *ex, const myScalar *ey, const myScalar *ez,
                     const myScalar *vec, myScalar *ans);
  inline void transpose_matvec(const myScalar mat[3][3], const myScalar *vec,
                               myScalar *ans);
  inline void transpose_matvec(const myScalar *ex, const myScalar *ey,
                               const myScalar *ez, const myScalar *v,
                               myScalar *ans);
  inline void transpose_diag3(const myScalar m[3][3], const myScalar *d,
                              myScalar ans[3][3]);
  inline void vecmat(const myScalar *v, const myScalar m[3][3], myScalar *ans);
  inline void scalar_times3(const myScalar f, myScalar m[3][3]);

  void write3(const myScalar mat[3][3]);
  int mldivide3(const myScalar mat[3][3], const myScalar *vec, myScalar *ans);
  int jacobi(myScalar matrix[3][3], myScalar *evalues, myScalar evectors[3][3]);
  void rotate(myScalar matrix[3][3], int i, int j, int k, int l,
              myScalar s, myScalar tau);
  void richardson(myScalar *q, myScalar *m, myScalar *w, myScalar *moments, myScalar dtq);
  void no_squish_rotate(int k, myScalar *p, myScalar *q, myScalar *inertia,
                        myScalar dt);

  // shape matrix operations
  // upper-triangular 3x3 matrix stored in Voigt notation as 6-vector

  inline void multiply_shape_shape(const myScalar *one, const myScalar *two,
                                   myScalar *ans);

  // quaternion operations

  inline void qnormalize(myScalar *q);
  inline void qconjugate(myScalar *q, myScalar *qc);
  inline void vecquat(myScalar *a, myScalar *b, myScalar *c);
  inline void quatvec(myScalar *a, myScalar *b, myScalar *c);
  inline void quatquat(myScalar *a, myScalar *b, myScalar *c);
  inline void invquatvec(myScalar *a, myScalar *b, myScalar *c);
  inline void axisangle_to_quat(const myScalar *v, const myScalar angle,
                                myScalar *quat);

  void angmom_to_omega(myScalar *m, myScalar *ex, myScalar *ey, myScalar *ez,
                       myScalar *idiag, myScalar *w);
  void omega_to_angmom(myScalar *w, myScalar *ex, myScalar *ey, myScalar *ez,
                       myScalar *idiag, myScalar *m);
  void mq_to_omega(myScalar *m, myScalar *q, myScalar *moments, myScalar *w);
  void exyz_to_q(myScalar *ex, myScalar *ey, myScalar *ez, myScalar *q);
  void q_to_exyz(myScalar *q, myScalar *ex, myScalar *ey, myScalar *ez);
  void quat_to_mat(const myScalar *quat, myScalar mat[3][3]);
  void quat_to_mat_trans(const myScalar *quat, myScalar mat[3][3]);

  // rotation operations

  inline void rotation_generator_x(const myScalar m[3][3], myScalar ans[3][3]);
  inline void rotation_generator_y(const myScalar m[3][3], myScalar ans[3][3]);
  inline void rotation_generator_z(const myScalar m[3][3], myScalar ans[3][3]);

  void BuildRxMatrix(myScalar R[3][3], const myScalar angle);
  void BuildRyMatrix(myScalar R[3][3], const myScalar angle);
  void BuildRzMatrix(myScalar R[3][3], const myScalar angle);

  // moment of inertia operations

  void inertia_ellipsoid(myScalar *shape, myScalar *quat, myScalar mass,
                         myScalar *inertia);
  void inertia_line(myScalar length, myScalar theta, myScalar mass,
                    myScalar *inertia);
  void inertia_triangle(myScalar *v0, myScalar *v1, myScalar *v2,
                        myScalar mass, myScalar *inertia);
  void inertia_triangle(myScalar *idiag, myScalar *quat, myScalar mass,
                        myScalar *inertia);
}

/* ----------------------------------------------------------------------
   copy a vector, return in ans
------------------------------------------------------------------------- */

inline void MathExtra::copy3(const myScalar *v, myScalar *ans)
{
  ans[0] = v[0];
  ans[1] = v[1];
  ans[2] = v[2];
}

/* ----------------------------------------------------------------------
   set vector equal to zero
------------------------------------------------------------------------- */

inline void MathExtra::zero3(myScalar *v)
{
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */

inline void MathExtra::norm3(myScalar *v)
{
  myScalar scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

inline void MathExtra::normalize3(const myScalar *v, myScalar *ans)
{
  myScalar scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

inline void MathExtra::snormalize3(const myScalar length, const myScalar *v,
                                   myScalar *ans)
{
  myScalar scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   negate vector v
------------------------------------------------------------------------- */

inline void MathExtra::negate3(myScalar *v)
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

/* ----------------------------------------------------------------------
   scale vector v by s
------------------------------------------------------------------------- */

inline void MathExtra::scale3(myScalar s, myScalar *v)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

/* ----------------------------------------------------------------------
   ans = v1 + v2
------------------------------------------------------------------------- */

inline void MathExtra::add3(const myScalar *v1, const myScalar *v2, myScalar *ans)
{
  ans[0] = v1[0] + v2[0];
  ans[1] = v1[1] + v2[1];
  ans[2] = v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = s*v1 + v2
------------------------------------------------------------------------- */

inline void MathExtra::scaleadd3(myScalar s, const myScalar *v1,
                                 const myScalar *v2, myScalar *ans)
{
  ans[0] = s*v1[0] + v2[0];
  ans[1] = s*v1[1] + v2[1];
  ans[2] = s*v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = v1 - v2
------------------------------------------------------------------------- */

inline void MathExtra::sub3(const myScalar *v1, const myScalar *v2, myScalar *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

inline myScalar MathExtra::len3(const myScalar *v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

inline myScalar MathExtra::lensq3(const myScalar *v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/* ----------------------------------------------------------------------
   ans = distance squared between pts v1 and v2
------------------------------------------------------------------------- */

inline myScalar MathExtra::distsq3(const myScalar *v1, const myScalar *v2)
{
  myScalar dx = v1[0] - v2[0];
  myScalar dy = v1[1] - v2[1];
  myScalar dz = v1[2] - v2[2];
  return dx*dx + dy*dy + dz*dz;
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

inline myScalar MathExtra::dot3(const myScalar *v1, const myScalar *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

inline void MathExtra::cross3(const myScalar *v1, const myScalar *v2, myScalar *ans)
{
  ans[0] = v1[1]*v2[2] - v1[2]*v2[1];
  ans[1] = v1[2]*v2[0] - v1[0]*v2[2];
  ans[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   construct matrix from 3 column vectors
------------------------------------------------------------------------- */

void MathExtra::col2mat(const myScalar *ex, const myScalar *ey, const myScalar *ez,
                        myScalar m[3][3])
{
  m[0][0] = ex[0];
  m[1][0] = ex[1];
  m[2][0] = ex[2];
  m[0][1] = ey[0];
  m[1][1] = ey[1];
  m[2][1] = ey[2];
  m[0][2] = ez[0];
  m[1][2] = ez[1];
  m[2][2] = ez[2];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

inline myScalar MathExtra::det3(const myScalar m[3][3])
{
  myScalar ans = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] -
    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] +
    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

inline void MathExtra::diag_times3(const myScalar *d, const myScalar m[3][3],
                                   myScalar ans[3][3])
{
  ans[0][0] = d[0]*m[0][0];
  ans[0][1] = d[0]*m[0][1];
  ans[0][2] = d[0]*m[0][2];
  ans[1][0] = d[1]*m[1][0];
  ans[1][1] = d[1]*m[1][1];
  ans[1][2] = d[1]*m[1][2];
  ans[2][0] = d[2]*m[2][0];
  ans[2][1] = d[2]*m[2][1];
  ans[2][2] = d[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   full matrix times a diagonal matrix
------------------------------------------------------------------------- */

void MathExtra::times3_diag(const myScalar m[3][3], const myScalar *d,
                            myScalar ans[3][3])
{
  ans[0][0] = m[0][0]*d[0];
  ans[0][1] = m[0][1]*d[1];
  ans[0][2] = m[0][2]*d[2];
  ans[1][0] = m[1][0]*d[0];
  ans[1][1] = m[1][1]*d[1];
  ans[1][2] = m[1][2]*d[2];
  ans[2][0] = m[2][0]*d[0];
  ans[2][1] = m[2][1]*d[1];
  ans[2][2] = m[2][2]*d[2];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

inline void MathExtra::plus3(const myScalar m[3][3], const myScalar m2[3][3],
                             myScalar ans[3][3])
{
  ans[0][0] = m[0][0]+m2[0][0];
  ans[0][1] = m[0][1]+m2[0][1];
  ans[0][2] = m[0][2]+m2[0][2];
  ans[1][0] = m[1][0]+m2[1][0];
  ans[1][1] = m[1][1]+m2[1][1];
  ans[1][2] = m[1][2]+m2[1][2];
  ans[2][0] = m[2][0]+m2[2][0];
  ans[2][1] = m[2][1]+m2[2][1];
  ans[2][2] = m[2][2]+m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

inline void MathExtra::times3(const myScalar m[3][3], const myScalar m2[3][3],
                              myScalar ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[1][0] + m[0][2]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[0][1]*m2[1][1] + m[0][2]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[0][1]*m2[1][2] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[1][0] + m[1][2]*m2[2][0];
  ans[1][1] = m[1][0]*m2[0][1] + m[1][1]*m2[1][1] + m[1][2]*m2[2][1];
  ans[1][2] = m[1][0]*m2[0][2] + m[1][1]*m2[1][2] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[2][0]*m2[0][1] + m[2][1]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[2][0]*m2[0][2] + m[2][1]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

inline void MathExtra::transpose_times3(const myScalar m[3][3],
                                        const myScalar m2[3][3],myScalar ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[1][0]*m2[1][0] + m[2][0]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[1][0]*m2[1][1] + m[2][0]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[1][0]*m2[1][2] + m[2][0]*m2[2][2];
  ans[1][0] = m[0][1]*m2[0][0] + m[1][1]*m2[1][0] + m[2][1]*m2[2][0];
  ans[1][1] = m[0][1]*m2[0][1] + m[1][1]*m2[1][1] + m[2][1]*m2[2][1];
  ans[1][2] = m[0][1]*m2[0][2] + m[1][1]*m2[1][2] + m[2][1]*m2[2][2];
  ans[2][0] = m[0][2]*m2[0][0] + m[1][2]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[0][2]*m2[0][1] + m[1][2]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[0][2]*m2[0][2] + m[1][2]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times transpose of mat2
------------------------------------------------------------------------- */

inline void MathExtra::times3_transpose(const myScalar m[3][3],
                                        const myScalar m2[3][3],myScalar ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[0][1] + m[0][2]*m2[0][2];
  ans[0][1] = m[0][0]*m2[1][0] + m[0][1]*m2[1][1] + m[0][2]*m2[1][2];
  ans[0][2] = m[0][0]*m2[2][0] + m[0][1]*m2[2][1] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[0][1] + m[1][2]*m2[0][2];
  ans[1][1] = m[1][0]*m2[1][0] + m[1][1]*m2[1][1] + m[1][2]*m2[1][2];
  ans[1][2] = m[1][0]*m2[2][0] + m[1][1]*m2[2][1] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[0][1] + m[2][2]*m2[0][2];
  ans[2][1] = m[2][0]*m2[1][0] + m[2][1]*m2[1][1] + m[2][2]*m2[1][2];
  ans[2][2] = m[2][0]*m2[2][0] + m[2][1]*m2[2][1] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT check for singular or badly scaled matrix
------------------------------------------------------------------------- */

inline void MathExtra::invert3(const myScalar m[3][3], myScalar ans[3][3])
{
  myScalar den = m[0][0]*m[1][1]*m[2][2]-m[0][0]*m[1][2]*m[2][1];
  den += -m[1][0]*m[0][1]*m[2][2]+m[1][0]*m[0][2]*m[2][1];
  den += m[2][0]*m[0][1]*m[1][2]-m[2][0]*m[0][2]*m[1][1];

  ans[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1]) / den;
  ans[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1]) / den;
  ans[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1]) / den;
  ans[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0]) / den;
  ans[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0]) / den;
  ans[1][2] = -(m[0][0]*m[1][2]-m[0][2]*m[1][0]) / den;
  ans[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0]) / den;
  ans[2][1] = -(m[0][0]*m[2][1]-m[0][1]*m[2][0]) / den;
  ans[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) / den;
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

inline void MathExtra::matvec(const myScalar m[3][3], const myScalar *v,
                              myScalar *ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

inline void MathExtra::matvec(const myScalar *ex, const myScalar *ey,
                              const myScalar *ez, const myScalar *v, myScalar *ans)
{
  ans[0] = ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2];
  ans[1] = ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2];
  ans[2] = ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

inline void MathExtra::transpose_matvec(const myScalar m[3][3], const myScalar *v,
                                 myScalar *ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

inline void MathExtra::transpose_matvec(const myScalar *ex, const myScalar *ey,
                                 const myScalar *ez, const myScalar *v,
                                 myScalar *ans)
{
  ans[0] = ex[0]*v[0] + ex[1]*v[1] + ex[2]*v[2];
  ans[1] = ey[0]*v[0] + ey[1]*v[1] + ey[2]*v[2];
  ans[2] = ez[0]*v[0] + ez[1]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

inline void MathExtra::transpose_diag3(const myScalar m[3][3], const myScalar *d,
                                myScalar ans[3][3])
{
  ans[0][0] = m[0][0]*d[0];
  ans[0][1] = m[1][0]*d[1];
  ans[0][2] = m[2][0]*d[2];
  ans[1][0] = m[0][1]*d[0];
  ans[1][1] = m[1][1]*d[1];
  ans[1][2] = m[2][1]*d[2];
  ans[2][0] = m[0][2]*d[0];
  ans[2][1] = m[1][2]*d[1];
  ans[2][2] = m[2][2]*d[2];
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

inline void MathExtra::vecmat(const myScalar *v, const myScalar m[3][3],
                              myScalar *ans)
{
  ans[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0];
  ans[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1];
  ans[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   matrix times scalar, in place
------------------------------------------------------------------------- */

inline void MathExtra::scalar_times3(const myScalar f, myScalar m[3][3])
{
  m[0][0] *= f; m[0][1] *= f; m[0][2] *= f;
  m[1][0] *= f; m[1][1] *= f; m[1][2] *= f;
  m[2][0] *= f; m[2][1] *= f; m[2][2] *= f;
}

/* ----------------------------------------------------------------------
   multiply 2 shape matrices
   upper-triangular 3x3, stored as 6-vector in Voigt notation
------------------------------------------------------------------------- */

inline void MathExtra::multiply_shape_shape(const myScalar *one,
                                            const myScalar *two, myScalar *ans)
{
  ans[0] = one[0]*two[0];
  ans[1] = one[1]*two[1];
  ans[2] = one[2]*two[2];
  ans[3] = one[1]*two[3] + one[3]*two[2];
  ans[4] = one[0]*two[4] + one[5]*two[3] + one[4]*two[2];
  ans[5] = one[0]*two[5] + one[5]*two[1];
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

inline void MathExtra::qnormalize(myScalar *q)
{
  myScalar norm = 1.0 / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

inline void MathExtra::qconjugate(myScalar *q, myScalar *qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   vector-quaternion multiply: c = a*b, where a = (0,a)
------------------------------------------------------------------------- */

inline void MathExtra::vecquat(myScalar *a, myScalar *b, myScalar *c)
{
  c[0] = -a[0]*b[1] - a[1]*b[2] - a[2]*b[3];
  c[1] = b[0]*a[0] + a[1]*b[3] - a[2]*b[2];
  c[2] = b[0]*a[1] + a[2]*b[1] - a[0]*b[3];
  c[3] = b[0]*a[2] + a[0]*b[2] - a[1]*b[1];
}

/* ----------------------------------------------------------------------
   quaternion-vector multiply: c = a*b, where b = (0,b)
------------------------------------------------------------------------- */

inline void MathExtra::quatvec(myScalar *a, myScalar *b, myScalar *c)
{
  c[0] = -a[1]*b[0] - a[2]*b[1] - a[3]*b[2];
  c[1] = a[0]*b[0] + a[2]*b[2] - a[3]*b[1];
  c[2] = a[0]*b[1] + a[3]*b[0] - a[1]*b[2];
  c[3] = a[0]*b[2] + a[1]*b[1] - a[2]*b[0];
}

/* ----------------------------------------------------------------------
   quaternion-quaternion multiply: c = a*b
------------------------------------------------------------------------- */

inline void MathExtra::quatquat(myScalar *a, myScalar *b, myScalar *c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  c[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
}

/* ----------------------------------------------------------------------
   quaternion multiply: c = inv(a)*b
   a is a quaternion
   b is a four component vector
   c is a three component vector
------------------------------------------------------------------------- */

inline void MathExtra::invquatvec(myScalar *a, myScalar *b, myScalar *c)
{
  c[0] = -a[1]*b[0] + a[0]*b[1] + a[3]*b[2] - a[2]*b[3];
  c[1] = -a[2]*b[0] - a[3]*b[1] + a[0]*b[2] + a[1]*b[3];
  c[2] = -a[3]*b[0] + a[2]*b[1] - a[1]*b[2] + a[0]*b[3];
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */

inline void MathExtra::axisangle_to_quat(const myScalar *v, const myScalar angle,
                                         myScalar *quat)
{
  myScalar halfa = 0.5*angle;
  myScalar sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about x to rotation matrix m
------------------------------------------------------------------------- */

inline void MathExtra::rotation_generator_x(const myScalar m[3][3],
                                            myScalar ans[3][3])
{
  ans[0][0] = 0;
  ans[0][1] = -m[0][2];
  ans[0][2] = m[0][1];
  ans[1][0] = 0;
  ans[1][1] = -m[1][2];
  ans[1][2] = m[1][1];
  ans[2][0] = 0;
  ans[2][1] = -m[2][2];
  ans[2][2] = m[2][1];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about y to rotation matrix m
------------------------------------------------------------------------- */

inline void MathExtra::rotation_generator_y(const myScalar m[3][3],
                                            myScalar ans[3][3])
{
  ans[0][0] = m[0][2];
  ans[0][1] = 0;
  ans[0][2] = -m[0][0];
  ans[1][0] = m[1][2];
  ans[1][1] = 0;
  ans[1][2] = -m[1][0];
  ans[2][0] = m[2][2];
  ans[2][1] = 0;
  ans[2][2] = -m[2][0];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about z to rotation matrix m
------------------------------------------------------------------------- */

inline void MathExtra::rotation_generator_z(const myScalar m[3][3],
                                            myScalar ans[3][3])
{
  ans[0][0] = -m[0][1];
  ans[0][1] = m[0][0];
  ans[0][2] = 0;
  ans[1][0] = -m[1][1];
  ans[1][1] = m[1][0];
  ans[1][2] = 0;
  ans[2][0] = -m[2][1];
  ans[2][1] = m[2][0];
  ans[2][2] = 0;
}

#endif
