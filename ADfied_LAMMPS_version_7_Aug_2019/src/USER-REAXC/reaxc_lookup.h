#include "ad_defines.h"
/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#ifndef __LOOKUP_H_
#define __LOOKUP_H_

#include "reaxc_types.h"
namespace LAMMPS_NS { class Error; }

void Tridiagonal_Solve( const myScalar *a, const myScalar *b,
                        myScalar *c, myScalar *d, myScalar *x, unsigned int n);

void Natural_Cubic_Spline( LAMMPS_NS::Error*, const myScalar *h, const myScalar *f,
                           cubic_spline_coef *coef, unsigned int n );

void Complete_Cubic_Spline( LAMMPS_NS::Error*, const myScalar *h, const myScalar *f,
                            myScalar v0, myScalar vlast, cubic_spline_coef *coef,
                            unsigned int n );

int Init_Lookup_Tables( reax_system*, control_params*, storage*,
                        mpi_datatypes*, char* );

void Deallocate_Lookup_Tables( reax_system* );

#endif
