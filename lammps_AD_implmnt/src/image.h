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

#ifndef LMP_IMAGE_H
#define LMP_IMAGE_H

#include <cmath>
#include "pointers.h"

namespace LAMMPS_NS {

class Image : protected Pointers {
 public:
  int width,height;             // size of image
  myScalar theta,phi;             // view image from theta,phi
  myScalar xctr,yctr,zctr;        // center of image in user coords
  myScalar up[3];                 // up direction in image
  myScalar zoom;                  // zoom factor
  myScalar persp;                 // perspective factor
  myScalar shiny;                 // shininess of objects
  int ssao;                     // SSAO on or off
  int seed;                     // RN seed for SSAO
  myScalar ssaoint;               // strength of shading from 0 to 1
  myScalar *boxcolor;             // color to draw box outline with
  int background[3];            // RGB values of background

  Image(class LAMMPS *, int);
  ~Image();
  void buffers();
  void clear();
  void merge();
  void write_JPG(FILE *);
  void write_PNG(FILE *);
  void write_PPM(FILE *);
  void view_params(myScalar, myScalar, myScalar, myScalar, myScalar, myScalar);

  void draw_sphere(myScalar *, myScalar *, myScalar);
  void draw_cube(myScalar *, myScalar *, myScalar);
  void draw_cylinder(myScalar *, myScalar *, myScalar *, myScalar, int);
  void draw_triangle(myScalar *, myScalar *, myScalar *, myScalar *);
  void draw_box(myScalar (*)[3], myScalar);
  void draw_axes(myScalar (*)[3], myScalar);

  int map_dynamic(int);
  int map_reset(int, int, char **);
  int map_minmax(int, myScalar, myScalar);
  myScalar *map_value2color(int, myScalar);

  int addcolor(char *, myScalar, myScalar, myScalar);
  myScalar *element2color(char *);
  myScalar element2diam(char *);
  myScalar *color2rgb(const char *, int index=0);
  int default_colors();

 private:
  int me,nprocs;
  int npixels;

  class ColorMap **maps;
  int nmap;

  myScalar *depthBuffer,*surfaceBuffer;
  myScalar *depthcopy,*surfacecopy;
  unsigned char *imageBuffer,*rgbcopy,*writeBuffer;

  // constant view params

  myScalar FOV;
  myScalar ambientColor[3];

  myScalar keyLightTheta;
  myScalar keyLightPhi;
  myScalar keyLightColor[3];

  myScalar fillLightTheta;
  myScalar fillLightPhi;
  myScalar fillLightColor[3];

  myScalar backLightTheta;
  myScalar backLightPhi;
  myScalar backLightColor[3];

  myScalar specularHardness;
  myScalar specularIntensity;

  myScalar SSAORadius;
  int SSAOSamples;
  myScalar SSAOJitter;

  // dynamic view params

  myScalar zdist;
  myScalar tanPerPixel;
  myScalar camDir[3],camUp[3],camRight[4],camPos[3];
  myScalar keyLightDir[3],fillLightDir[3],backLightDir[3];
  myScalar keyHalfDir[3];

  // color values

  int ncolors;
  char **username;
  myScalar **userrgb;

  // SSAO RNG

  class RanMars *random;

  // internal methods

  void draw_pixel(int, int, myScalar, myScalar *, myScalar*);
  void compute_SSAO();

  // inline functions

  inline myScalar saturate(myScalar v) {
    if (v < 0.0) return 0.0;
    else if (v > 1.0) return 1.0;
    else return v;
  }

  inline myScalar distance(myScalar* a, myScalar* b) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) +
                (a[1] - b[1]) * (a[1] - b[1]) +
                (a[2] - b[2]) * (a[2] - b[2]));
  }
};

// ColorMap class

class ColorMap : protected Pointers {
 public:
  int dynamic;                     // 0/1 if lo/hi bounds are static/dynamic

  ColorMap(class LAMMPS *, class Image*);
  ~ColorMap();
  int reset(int, char **);
  int minmax(myScalar, myScalar);
  myScalar *value2color(myScalar);

 private:
  class Image *image;              // caller with color2rgb() method
  int mstyle,mrange;               // 2-letter style/range of color map
  int mlo,mhi;                     // bounds = NUMERIC or MINVALUE or MAXVALUE
  myScalar mlovalue,mhivalue;        // user bounds if NUMERIC
  myScalar locurrent,hicurrent;      // current bounds for this snapshot
  myScalar mbinsize,mbinsizeinv;     // bin size for sequential color map
  myScalar interpolate[3];           // local storage for returned RGB color

  struct MapEntry {
    int single,lo,hi;              // NUMERIC or MINVALUE or MAXVALUE
    myScalar svalue,lvalue,hvalue;   // actual value
    myScalar *color;                 // RGB values
  };

  MapEntry *mentry;
  int nentry;
};

}

#endif

/* ERROR/WARNING messages:

E: Invalid image up vector

Up vector cannot be (0,0,0).

*/
