#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SETOPERATORCOEFFICIENTS_H_
#define _SETOPERATORCOEFFICIENTS_H_

#include <iostream>
using std::cerr;
#include "ParmParse.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "BRMeshRefine.H"
#include "FABView.H"
#include "DebugDump.H"
#include "VClocalFuncs.H"
#include "PoissonParameters.H"
#include "MultilevelLinearOp.H"
#include "BiCGStabSolver.H"
#include "AMRIO.H"
#include "UsingNamespace.H"

/* see BoxTools/BaseFabMacros.H  for ForAllXBNN loop macro */
void
setACoef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
    RealVect pos;
    int num;

    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // cell centered
        RealVect ccOffset = 0.5*a_dx*RealVect::Unit; 

        FArrayBox& aCoef = a_aCoef[dit];
        FArrayBox& chi = a_chi[dit];
        Box thisBox = aCoef.box();
     
        BoxIterator bit(thisBox);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv);
            loc *= a_dx;
            loc += ccOffset;

            Real M = M_value(loc, a_params);

            Real chi_0 = chi(iv,0);

            aCoef(iv,0) = (-0.625) * M * (chi_0 * chi_0 * chi_0 * chi_0);
        }
    }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
void
setBCoef(LevelData<FArrayBox>& a_bCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
  RealVect pos;
  int num;

  DataIterator dit = a_bCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& bCoef = a_bCoef[dit];
      ForAllXBNN(Real, bCoef, bCoef.box(), 0, bCoef.nComp());
      {
        bCoefR = 1.0; // constant 
      }
    }EndFor;
}

#endif
