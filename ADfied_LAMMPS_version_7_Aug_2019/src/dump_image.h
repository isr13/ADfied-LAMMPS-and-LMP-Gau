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

#ifdef DUMP_CLASS

DumpStyle(image,DumpImage)

#else

#ifndef LMP_DUMP_IMAGE_H
#define LMP_DUMP_IMAGE_H

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpImage : public DumpCustom {
 public:
  int multifile_override;          // used by write_dump command

  DumpImage(class LAMMPS *, int, char**);
  virtual ~DumpImage();
  int pack_forward_comm(int, int *, myScalar *, int, int *);
  void unpack_forward_comm(int, int, myScalar *);

 protected:
  int filetype;
  enum{PPM,JPG,PNG};

  int atomflag;                    // 0/1 for draw atoms
  int acolor,adiam;                // what determines color/diam of atoms
  myScalar adiamvalue;               // atom diameter value

  int lineflag;                    // 0/1 for draw atoms as lines
  int lcolor,ldiam;                // what determines color/diam of lines
  myScalar ldiamvalue;               // line diameter value
  int triflag;                     // 0/1 for draw atoms as triangles
  int tcolor,tstyle;               // what determines color/style of tris
  myScalar tdiamvalue;               // tri edge diameter value
  int bodyflag;                    // 0/1 for draw atoms as bodies
  int bodycolor;                   // what determines color of bodies
  myScalar bodyflag1,bodyflag2;      // user-specified params for drawing bodies
  int fixflag;                     // 0/1 to draw what fix provides
  int fixcolor;                    // what determines color of fix objects
  myScalar fixflag1,fixflag2;        // user-specified params for fix objects

  int bondflag;                    // 0/1 for draw bonds
  int bcolor,bdiam;                // what determines color/diam of bonds
  myScalar bdiamvalue;               // bond diameter value

  int extraflag;                   // 0/1 for any of line/tri/body flag set
  char *thetastr,*phistr;          // variables for view theta,phi
  int thetavar,phivar;             // index to theta,phi vars
  int cflag;                       // static/dynamic box center
  myScalar cx,cy,cz;                 // fractional box center
  char *cxstr,*cystr,*czstr;       // variables for box center
  int cxvar,cyvar,czvar;           // index to box center vars
  char *upxstr,*upystr,*upzstr;    // view up vector variables
  int upxvar,upyvar,upzvar;        // index to up vector vars
  char *zoomstr,*perspstr;         // view zoom and perspective variables
  int zoomvar,perspvar;            // index to zoom,persp vars
  int boxflag,axesflag;            // 0/1 for draw box and axes
  myScalar boxdiam,axeslen,axesdiam; // params for drawing box and axes
  int subboxflag;
  myScalar subboxdiam;

  int viewflag;                    // overall view is static or dynamic

  myScalar *diamtype,*diamelement,*bdiamtype;         // per-type diameters
  myScalar **colortype,**colorelement,**bcolortype;   // per-type colors

  class AtomVecLine *avec_line;    // ptrs to atom style (sub)classes
  class AtomVecTri *avec_tri;
  class AtomVecBody *avec_body;

  class Fix *fixptr;               // ptr to Fix that provides image data

  class Image *image;              // class that renders each image

  int *chooseghost;                // extended choose array for comm
  myScalar **bufcopy;                // buffer for communicating bond/atom info
  int maxbufcopy;

  virtual void init_style();
  int modify_param(int, char **);
  void write();

  void box_center();
  void view_params();
  void box_bounds();

  void create_image();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump image filename

The file produced by dump image cannot be binary and must
be for a single processor.

E: Support for writing images in JPEG format not included

LAMMPS was not built with the -DLAMMPS_JPEG switch in the Makefile.

E: Support for writing images in PNG format not included

LAMMPS was not built with the -DLAMMPS_PNG switch in the Makefile.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump image bond not allowed with no bond types

Self-explanatory.

E: Invalid dump image theta value

Theta must be between 0.0 and 180.0 inclusive.

E: Dump image persp option is not yet supported

Self-explanatory.

E: Dump image line requires atom style line

Self-explanatory.

E: Dump image tri requires atom style tri

Self-explanatory.

E: Dump image body yes requires atom style body

Self-explanatory.

E: Fix ID for dump image does not exist

UNDOCUMENTED

E: Dump image requires one snapshot per file

Use a "*" in the filename.

E: Dump image cannot perform sorting

Self-explanatory.

E: Variable name for dump image theta does not exist

Self-explanatory.

E: Variable for dump image theta is invalid style

Must be an equal-style variable.

E: Variable name for dump image phi does not exist

Self-explanatory.

E: Variable for dump image phi is invalid style

Must be an equal-style variable.

E: Variable name for dump image center does not exist

Self-explanatory.

E: Variable for dump image center is invalid style

Must be an equal-style variable.

E: Variable name for dump image zoom does not exist

Self-explanatory.

E: Variable for dump image zoom is invalid style

Must be an equal-style variable.

E: Variable name for dump image persp does not exist

Self-explanatory.

E: Variable for dump image persp is invalid style

Must be an equal-style variable.

E: Invalid dump image element name

The specified element name was not in the standard list of elements.
See the dump_modify doc page.

E: Invalid color map min/max values

The min/max values are not consistent with either each other or
with values in the color map.

E: Invalid dump image zoom value

Zoom value must be > 0.0.

E: Invalid dump image persp value

Persp value must be >= 0.0.

E: Invalid color in dump_modify command

The specified color name was not in the list of recognized colors.
See the dump_modify doc page.

E: Dump modify bcolor not allowed with no bond types

Self-explanatory.

E: Dump modify bdiam not allowed with no bond types

Self-explanatory.

*/
