/* ----------------------------------------------------------------------
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
   Contributing author:  Aidan Thompson (SNL)
                         Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "compute_orientorder_atom.h"
#include <cstring>
#include <cstdlib>
#include <cmath>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtom::ComputeOrientOrderAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  qlist(NULL), distsq(NULL), nearest(NULL), rlist(NULL),
  qnarray(NULL), qnm_r(NULL), qnm_i(NULL)
{
  if (narg < 3 ) error->all(FLERR,"Illegal compute orientorder/atom command");

  // set default values for optional args

  nnn = 12;
  cutsq = 0.0;
  qlcompflag = 0;

  // specify which orders to request

  nqlist = 5;
  memory->create(qlist,nqlist,"orientorder/atom:qlist");
  qlist[0] = 4;
  qlist[1] = 6;
  qlist[2] = 8;
  qlist[3] = 10;
  qlist[4] = 12;
  qmax = 12;

  // process optional args

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"nnn") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      if (strcmp(arg[iarg+1],"NULL") == 0) {
        nnn = 0;
      } else {
        nnn = ad::passive_value(force->numeric(FLERR,arg[iarg+1]));
        if (nnn <= 0)
          error->all(FLERR,"Illegal compute orientorder/atom command");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"degrees") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      nqlist = ad::passive_value(force->numeric(FLERR,arg[iarg+1]));
      if (nqlist <= 0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      memory->destroy(qlist);
      memory->create(qlist,nqlist,"orientorder/atom:qlist");
      iarg += 2;
      if (iarg+nqlist > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      qmax = 0;
      for (int iw = 0; iw < nqlist; iw++) {
        qlist[iw] = ad::passive_value(force->numeric(FLERR,arg[iarg+iw]));
        if (qlist[iw] < 0)
          error->all(FLERR,"Illegal compute orientorder/atom command");
        if (qlist[iw] > qmax) qmax = qlist[iw];
      }
      iarg += nqlist;
    } else if (strcmp(arg[iarg],"components") == 0) {
      qlcompflag = 1;
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      qlcomp = ad::passive_value(force->numeric(FLERR,arg[iarg+1]));
      if (qlcomp <= 0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      iqlcomp = -1;
      for (int iw = 0; iw < nqlist; iw++)
        if (qlcomp == qlist[iw]) {
          iqlcomp = iw;
          break;
        }
      if (iqlcomp < 0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      myScalar cutoff = force->numeric(FLERR,arg[iarg+1]);
      if (cutoff <= 0.0)
        error->all(FLERR,"Illegal compute orientorder/atom command");
      cutsq = cutoff*cutoff;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute orientorder/atom command");
  }

  if (qlcompflag) ncol = nqlist + 2*(2*qlcomp+1);
  else ncol = nqlist;

  peratom_flag = 1;
  size_peratom_cols = ncol;

  nmax = 0;
  maxneigh = 0;
}

/* ---------------------------------------------------------------------- */

ComputeOrientOrderAtom::~ComputeOrientOrderAtom()
{
  memory->destroy(qnarray);
  memory->destroy(distsq);
  memory->destroy(rlist);
  memory->destroy(nearest);
  memory->destroy(qlist);
  memory->destroy(qnm_r);
  memory->destroy(qnm_i);

}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute orientorder/atom requires a "
               "pair style be defined");
  if (cutsq == 0.0) cutsq = force->pair->cutforce * force->pair->cutforce;
  else if (sqrt(cutsq) > force->pair->cutforce)
    error->all(FLERR,"Compute orientorder/atom cutoff is "
               "longer than pairwise cutoff");

  memory->create(qnm_r,qmax,2*qmax+1,"orientorder/atom:qnm_r");
  memory->create(qnm_i,qmax,2*qmax+1,"orientorder/atom:qnm_i");

  // need an occasional full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"orientorder/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute orientorder/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::compute_peratom()
{
  int i,j,ii,jj,inum,jnum;
  myScalar xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow order parameter array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(qnarray);
    nmax = atom->nmax;
    memory->create(qnarray,nmax,ncol,"orientorder/atom:qnarray");
    array_atom = qnarray;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute order parameter for each atom in group
  // use full neighbor list to count atoms less than cutoff

  myScalar **x = atom->x;
  int *mask = atom->mask;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    myScalar* qn = qnarray[i];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // insure distsq and nearest arrays are long enough

      if (jnum > maxneigh) {
        memory->destroy(distsq);
        memory->destroy(rlist);
        memory->destroy(nearest);
        maxneigh = jnum;
        memory->create(distsq,maxneigh,"orientorder/atom:distsq");
        memory->create(rlist,maxneigh,3,"orientorder/atom:rlist");
        memory->create(nearest,maxneigh,"orientorder/atom:nearest");
      }

      // loop over list of all neighbors within force cutoff
      // distsq[] = distance sq to each
      // rlist[] = distance vector to each
      // nearest[] = atom indices of neighbors

      int ncount = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          distsq[ncount] = rsq;
          rlist[ncount][0] = delx;
          rlist[ncount][1] = dely;
          rlist[ncount][2] = delz;
          nearest[ncount++] = j;
        }
      }

      // if not nnn neighbors, order parameter = 0;

      if ((ncount == 0) || (ncount < nnn)) {
        for (int iw = 0; iw < nqlist; iw++)
          qn[iw] = 0.0;
        continue;
      }

      // if nnn > 0, use only nearest nnn neighbors

      if (nnn > 0) {
        select3(nnn,ncount,distsq,nearest,rlist);
        ncount = nnn;
      }

      calc_boop(rlist, ncount, qn, qlist, nqlist);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

myScalar ComputeOrientOrderAtom::memory_usage()
{
  myScalar bytes = ncol*nmax * sizeof(myScalar);
  bytes += (qmax*(2*qmax+1)+maxneigh*4) * sizeof(myScalar);
  bytes += (nqlist+maxneigh) * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   select3 routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   sort auxiliary arrays at same time
------------------------------------------------------------------------- */

// Use no-op do while to create single statement

#define SWAP(a,b) do {       \
    tmp = a; a = b; b = tmp; \
  } while(0)

#define ISWAP(a,b) do {        \
    itmp = a; a = b; b = itmp; \
  } while(0)

#define SWAP3(a,b) do {                  \
    tmp = a[0]; a[0] = b[0]; b[0] = tmp; \
    tmp = a[1]; a[1] = b[1]; b[1] = tmp; \
    tmp = a[2]; a[2] = b[2]; b[2] = tmp; \
  } while(0)

/* ---------------------------------------------------------------------- */

void ComputeOrientOrderAtom::select3(int k, int n, myScalar *arr, int *iarr, myScalar **arr3)
{
  int i,ir,j,l,mid,ia,itmp;
  myScalar a,tmp,a3[3];

  arr--;
  iarr--;
  arr3--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir]);
        ISWAP(iarr[l],iarr[ir]);
        SWAP3(arr3[l],arr3[ir]);
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1]);
      ISWAP(iarr[mid],iarr[l+1]);
      SWAP3(arr3[mid],arr3[l+1]);
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
        ISWAP(iarr[l],iarr[ir]);
        SWAP3(arr3[l],arr3[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
        ISWAP(iarr[l+1],iarr[ir]);
        SWAP3(arr3[l+1],arr3[ir]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
        ISWAP(iarr[l],iarr[l+1]);
        SWAP3(arr3[l],arr3[l+1]);
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      a3[0] = arr3[l+1][0];
      a3[1] = arr3[l+1][1];
      a3[2] = arr3[l+1][2];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j]);
        ISWAP(iarr[i],iarr[j]);
        SWAP3(arr3[i],arr3[j]);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      arr3[l+1][0] = arr3[j][0];
      arr3[l+1][1] = arr3[j][1];
      arr3[l+1][2] = arr3[j][2];
      arr3[j][0] = a3[0];
      arr3[j][1] = a3[1];
      arr3[j][2] = a3[2];
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate the bond orientational order parameters
------------------------------------------------------------------------- */

void ComputeOrientOrderAtom::calc_boop(myScalar **rlist,
                                       int ncount, myScalar qn[],
                                       int qlist[], int nqlist) {
  for (int iw = 0; iw < nqlist; iw++) {
    int n = qlist[iw];

    qn[iw] = 0.0;
    for(int m = 0; m < 2*n+1; m++) {
      qnm_r[iw][m] = 0.0;
      qnm_i[iw][m] = 0.0;
    }
  }

  for(int ineigh = 0; ineigh < ncount; ineigh++) {
    const myScalar * const r = rlist[ineigh];
    myScalar rmag = dist(r);
    if(rmag <= MY_EPSILON) {
      return;
    }

    myScalar costheta = r[2] / rmag;
    myScalar expphi_r = r[0];
    myScalar expphi_i = r[1];
    myScalar rxymag = sqrt(expphi_r*expphi_r+expphi_i*expphi_i);
    if(rxymag <= MY_EPSILON) {
      expphi_r = 1.0;
      expphi_i = 0.0;
    } else {
      myScalar rxymaginv = 1.0/rxymag;
      expphi_r *= rxymaginv;
      expphi_i *= rxymaginv;
    }

    for (int iw = 0; iw < nqlist; iw++) {
      int n = qlist[iw];

      qnm_r[iw][n] += polar_prefactor(n, 0, costheta);
      myScalar expphim_r = expphi_r;
      myScalar expphim_i = expphi_i;
      for(int m = 1; m <= +n; m++) {
        myScalar prefactor = polar_prefactor(n, m, costheta);
        myScalar c_r = prefactor * expphim_r;
        myScalar c_i = prefactor * expphim_i;
        qnm_r[iw][m+n] += c_r;
        qnm_i[iw][m+n] += c_i;
        if(m & 1) {
          qnm_r[iw][-m+n] -= c_r;
          qnm_i[iw][-m+n] += c_i;
        } else {
          qnm_r[iw][-m+n] += c_r;
          qnm_i[iw][-m+n] -= c_i;
        }
        myScalar tmp_r = expphim_r*expphi_r - expphim_i*expphi_i;
        myScalar tmp_i = expphim_r*expphi_i + expphim_i*expphi_r;
        expphim_r = tmp_r;
        expphim_i = tmp_i;
      }

    }
  }

  myScalar fac = sqrt(MY_4PI) / ncount;
  myScalar normfac = 0.0;
  for (int iw = 0; iw < nqlist; iw++) {
    int n = qlist[iw];
    myScalar qm_sum = 0.0;
    for(int m = 0; m < 2*n+1; m++) {
      qm_sum += qnm_r[iw][m]*qnm_r[iw][m] + qnm_i[iw][m]*qnm_i[iw][m];
      //      printf("Ylm^2 = %d %d %g\n",n,m,
      //     qnm_r[iw][m]*qnm_r[iw][m] + qnm_i[iw][m]*qnm_i[iw][m]);
    }
    qn[iw] = fac * sqrt(qm_sum / (2*n+1));
    if (qlcompflag && iqlcomp == iw) normfac = 1.0/sqrt(qm_sum);

  }

  // output of the complex vector

  if (qlcompflag) {
    int j = nqlist;
    for(int m = 0; m < 2*qlcomp+1; m++) {
      qn[j++] = qnm_r[iqlcomp][m] * normfac;
      qn[j++] = qnm_i[iqlcomp][m] * normfac;
    }
  }
}

/* ----------------------------------------------------------------------
   calculate scalar distance
------------------------------------------------------------------------- */

myScalar ComputeOrientOrderAtom::dist(const myScalar r[])
{
  return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

/* ----------------------------------------------------------------------
   polar prefactor for spherical harmonic Y_l^m, where
   Y_l^m (theta, phi) = prefactor(l, m, cos(theta)) * exp(i*m*phi)
------------------------------------------------------------------------- */

myScalar ComputeOrientOrderAtom::polar_prefactor(int l, int m, myScalar costheta)
{
  const int mabs = abs(m);

  myScalar prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= static_cast<myScalar>(i);

  prefactor = sqrt(static_cast<myScalar>(2*l+1)/(MY_4PI*prefactor))
    * associated_legendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
------------------------------------------------------------------------- */

myScalar ComputeOrientOrderAtom::associated_legendre(int l, int m, myScalar x)
{
  if (l < m) return 0.0;

  myScalar p(1.0), pm1(0.0), pm2(0.0);

  if (m != 0) {
    const myScalar sqx = sqrt(1.0-x*x);
    for (int i=1; i < m+1; ++i)
      p *= static_cast<myScalar>(2*i-1) * sqx;
  }

  for (int i=m+1; i < l+1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = (static_cast<myScalar>(2*i-1)*x*pm1
         - static_cast<myScalar>(i+m-1)*pm2) / static_cast<myScalar>(i-m);
  }

  return p;
}
