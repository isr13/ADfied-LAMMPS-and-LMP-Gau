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

#ifndef LMP_REGION_H
#define LMP_REGION_H

#include "pointers.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class Region : protected Pointers {
 public:
  char *id,*style;
  int interior;                     // 1 for interior, 0 for exterior
  int scaleflag;                    // 1 for lattice, 0 for box
  myScalar xscale,yscale,zscale;      // scale factors for box/lattice units
  myScalar extent_xlo,extent_xhi;     // bounding box on region
  myScalar extent_ylo,extent_yhi;
  myScalar extent_zlo,extent_zhi;
  int bboxflag;                     // 1 if bounding box is computable
  int varshape;                     // 1 if region shape changes over time
  int dynamic;                      // 1 if position/orient changes over time
  int moveflag,rotateflag;          // 1 if position/orientation changes
  int openflag;             // 1 if any face is open
  int open_faces[6];            // flags for which faces are open

  int copymode;                     // 1 if copy of original class

  // contact = particle near region surface (for soft interactions)
  // touch = particle touching region surface (for granular interactions)

  struct Contact {
    myScalar r;                 // distance between particle & surf, r > 0.0
    myScalar delx,dely,delz;    // vector from surface pt to particle
    myScalar radius;            // curvature of region at contact point
    int iwall;            // unique id of wall for storing shear history
    int varflag;              // 1 if wall can be variable-controlled
  };
  Contact *contact;           // list of contacts
  int cmax;                   // max # of contacts possible with region
  int tmax;           // max # of touching contacts possible

  // motion attributes of region
  // public so can be accessed by other classes

  myScalar dx,dy,dz,theta;      // current displacement and orientation
  myScalar v[3];            // translational velocity
  myScalar rpoint[3];       // current origin of rotation axis
  myScalar omega[3];        // angular velocity
  myScalar rprev;               // speed of time-dependent radius, if applicable
  myScalar xcenter[3];          // translated/rotated center of cylinder/sphere (only used if varshape)
  myScalar prev[5];             // stores displacement (X3), angle and if
                              //  necessary, region variable size (e.g. radius)
                              //  at previous time step
  int vel_timestep;           // store timestep at which set_velocity was called
                              //   prevents multiple fix/wall/gran/region calls
  int nregion;                // For union and intersect
  int size_restart;
  int *list;

  Region(class LAMMPS *, int, char **);
  virtual ~Region();
  virtual void init();
  int dynamic_check();

  // called by other classes to check point versus region

  void prematch();
  int match(myScalar, myScalar, myScalar);
  int surface(myScalar, myScalar, myScalar, myScalar);

  virtual void set_velocity();
  void velocity_contact(myScalar *, myScalar *, int);
  virtual void write_restart(FILE *);
  virtual int restart(char *, int&);
  virtual void length_restart_string(int&);
  virtual void reset_vel();

  // implemented by each region, not called by other classes

  virtual int inside(myScalar, myScalar, myScalar) = 0;
  virtual int surface_interior(myScalar *, myScalar) = 0;
  virtual int surface_exterior(myScalar *, myScalar) = 0;
  virtual void shape_update() {}
  virtual void pretransform();
  virtual void set_velocity_shape() {}
  virtual void velocity_contact_shape(myScalar*, myScalar*) {}

 protected:
  void add_contact(int, myScalar *, myScalar, myScalar, myScalar);
  void options(int, char **);
  void point_on_line_segment(myScalar *, myScalar *, myScalar *, myScalar *);
  void forward_transform(myScalar &, myScalar &, myScalar &);
  myScalar point[3],runit[3];

 private:
  char *xstr,*ystr,*zstr,*tstr;
  int xvar,yvar,zvar,tvar;
  myScalar axis[3];

  void inverse_transform(myScalar &, myScalar &, myScalar &);
  void rotate(myScalar &, myScalar &, myScalar &, myScalar);
};

}

#endif

/* ERROR/WARNING messages:

E: Variable name for region does not exist

Self-explanatory.

E: Variable for region is invalid style

Only equal-style variables can be used.

E: Variable for region is not equal style

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region union or intersect cannot be dynamic

The sub-regions can be dynamic, but not the combined region.

E: Region cannot have 0 length rotation vector

Self-explanatory.

*/
