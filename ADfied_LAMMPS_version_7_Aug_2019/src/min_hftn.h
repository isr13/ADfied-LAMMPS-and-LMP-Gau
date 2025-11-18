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

#ifdef MINIMIZE_CLASS

MinimizeStyle(hftn,MinHFTN)

#else

#ifndef LMP_MIN_HFTN_H
#define LMP_MIN_HFTN_H

#include "min.h"

namespace LAMMPS_NS
{

class MinHFTN : public Min
{
  public:

  MinHFTN (LAMMPS *);
  ~MinHFTN (void);
  void init();
  void setup_style();
  void reset_vectors();
  int iterate (int);

  private:

  //---- THE ALGORITHM NEEDS TO STORE THIS MANY ATOM-BASED VECTORS,
  //---- IN ADDITION TO ATOM POSITIONS AND THE FORCE VECTOR.
  enum {
    VEC_XK=0,   //-- ATOM POSITIONS AT SUBITER START
    VEC_CG_P,   //-- STEP p IN CG SUBITER
    VEC_CG_D,   //-- DIRECTION d IN CG SUBITER
    VEC_CG_HD,  //-- HESSIAN-VECTOR PRODUCT Hd
    VEC_CG_R,   //-- RESIDUAL r IN CG SUBITER
    VEC_DIF1,   //-- FOR FINITE DIFFERENCING
    VEC_DIF2,   //-- FOR FINITE DIFFERENCING
    NUM_HFTN_ATOM_BASED_VECTORS
  };

  //---- ATOM-BASED STORAGE VECTORS.
  myScalar *   _daAVectors[NUM_HFTN_ATOM_BASED_VECTORS];
  myScalar **  _daExtraAtom[NUM_HFTN_ATOM_BASED_VECTORS];

  //---- GLOBAL DOF STORAGE.  ELEMENT [0] IS X0 (XK), NOT USED IN THIS ARRAY.
  myScalar *   _daExtraGlobal[NUM_HFTN_ATOM_BASED_VECTORS];

  int     _nNumUnknowns;
  FILE *  _fpPrint;

  int   execute_hftn_ (const bool      bPrintProgress,
                       const myScalar    dInitialEnergy,
                       const myScalar    dInitialForce2,
                       myScalar &  dFinalEnergy,
                       myScalar &  dFinalForce2);
  bool    compute_inner_cg_step_ (const myScalar    dTrustRadius,
                                  const myScalar    dForceTol,
                                  const int       nMaxEvals,
                                  const bool      bHaveEvalAtXin,
                                  const myScalar    dEnergyAtXin,
                                  const myScalar    dForce2AtXin,
                                  myScalar &  dEnergyAtXout,
                                  myScalar &  dForce2AtXout,
                                  int    &  nStepType,
                                  myScalar &  dStepLength2,
                                  myScalar &  dStepLengthInf);
  myScalar  calc_xinf_using_mpi_ (void) const;
  myScalar  calc_dot_prod_using_mpi_ (const int  nIx1,
                                    const int  nIx2) const;
  myScalar  calc_grad_dot_v_using_mpi_ (const int  nIx) const;
  void    calc_dhd_dd_using_mpi_ (myScalar &  dDHD,
                                  myScalar &  dDD) const;
  void    calc_ppnew_pdold_using_mpi_ (myScalar &  dPnewDotPnew,
                                       myScalar &  dPoldDotD) const;
  void    calc_plengths_using_mpi_ (myScalar &  dStepLength2,
                                    myScalar &  dStepLengthInf) const;
  bool    step_exceeds_TR_ (const myScalar    dTrustRadius,
                            const myScalar    dPP,
                            const myScalar    dPD,
                            const myScalar    dDD,
                            myScalar &  dTau) const;
  bool    step_exceeds_DMAX_ (void) const;
  void    adjust_step_to_tau_ (const myScalar  tau);
  myScalar  compute_to_tr_ (const myScalar  dPP,
                          const myScalar  dPD,
                          const myScalar  dDD,
                          const myScalar  dTrustRadius,
                          const bool    bConsiderBothRoots,
                          const myScalar  dDHD,
                          const myScalar  dPdotHD,
                          const myScalar  dGradDotD) const;
  void    evaluate_dir_der_ (const bool      bUseForwardDiffs,
                             const int       nIxDir,
                             const int       nIxResult,
                             const bool      bEvaluateAtX,
                             myScalar &  dNewEnergy);
  void  open_hftn_print_file_ (void);
  void  hftn_print_line_ (const bool    bIsStepAccepted,
                          const int     nIteration,
                          const int     nTotalEvals,
                          const myScalar  dEnergy,
                          const myScalar  dForce2,
                          const int     nStepType,
                          const myScalar  dTrustRadius,
                          const myScalar  dStepLength2,
                          const myScalar  dActualRed,
                          const myScalar  dPredictedRed) const;
  void  close_hftn_print_file_ (void);
};

}

#endif
#endif
