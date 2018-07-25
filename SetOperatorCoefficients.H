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
#include "SetOperatorCoefficients.H"
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

  if (a_params.ACoeftype == 0) // A(x) = x
    {
    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& aCoef = a_aCoef[dit];
        ForAllXBNN(Real, aCoef, aCoef.box(), 0, aCoef.nComp());  // Macros in BoxTools/BaseFabMacros.H
        {
          num = nR;
          D_TERM(pos[0]=a_dx[0]*(iR+0.5);,
                pos[1]=a_dx[1]*(jR+0.5);,
                pos[2]=a_dx[2]*(kR+0.5));
          aCoefR = pos[0];// this just fixed the aCoef to be x (which is what was also set in ChiFunctionsF.ChF
        }EndFor;
      } // end loop over grids
    }
  else if (a_params.ACoeftype == 1) // ACoef is constant 1
    {
      DataIterator dit = a_aCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& aCoef = a_aCoef[dit];
          ForAllXBNN(Real, aCoef, aCoef.box(), 0, aCoef.nComp());
          {
            aCoefR = 1.0; // constant 
          }
        }EndFor;
    } 
  else if (a_params.ACoeftype == 2) // ACoef = -5/8 * chi_0^4 M(chi_0,K) 
    {
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
    } // end ACoef
}

void
setBCoef(LevelData<FluxBox>& a_bCoef,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
  if (a_params.BCoeftype == 0 ) // original B(x)=(x,y,z)
    {
    DataIterator dit = a_bCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = a_bCoef[dit];
        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& dirFlux = thisBCoef[dir];
            const Box& dirBox = dirFlux.box();
            // this sets up a vector which is 0 in the dir
            // direct and 0.5 in the other (cell-centered) directions
            RealVect offsets = BASISREALV(dir);
            RealVect pos;
            offsets -= RealVect::Unit;
            offsets *= -0.5;
            int n;
            ForAllXBNN(Real, dirFlux, dirBox, 0, dirFlux.nComp())
              {
                n = nR;
                D_TERM(pos[0] = a_dx[0]*(iR+offsets[0]);,
                      pos[1] = a_dx[1]*(jR+offsets[1]);,
                      pos[2] = a_dx[2]*(kR+offsets[2]));
                dirFluxR = D_TERM(pos[0], +pos[1], +pos[2]);
                // constant-coefficient
                //dirFluxR = 1.0;
              }EndFor
          } // end loop over directions
      }
    }
  else if (a_params.BCoeftype == 1) // B(x) = (1,1,1)
    {
    DataIterator dit = a_bCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = a_bCoef[dit];
        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& dirFlux = thisBCoef[dir];
            const Box& dirBox = dirFlux.box();
            // this sets up a vector which is 0 in the dir
            // direct and 0.5 in the other (cell-centered) directions
            ForAllXBNN(Real, dirFlux, dirBox, 0, dirFlux.nComp())
              {
                dirFluxR = 1.0; // constant direction
              }EndFor
          } // end loop over directions
      }
    } // end BCoeftype
}

void
setCCoef(LevelData<FArrayBox>& a_cCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
  RealVect pos;
  int num;

  if (a_params.CCoeftype == 0) // C(x) = x
    {
    DataIterator dit = a_cCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& cCoef = a_cCoef[dit];
        ForAllXBNN(Real, cCoef, cCoef.box(), 0, cCoef.nComp());  // Macros in BoxTools/BaseFabMacros.H
        {
          num = nR;
          D_TERM(pos[0]=a_dx[0]*(iR+0.5);,
                pos[1]=a_dx[1]*(jR+0.5);,
                pos[2]=a_dx[2]*(kR+0.5));
          cCoefR = pos[0];// this just fixed the aCoef to be x (which is what was also set in ChiFunctionsF.Chf:w
          // Eugene change 
  //        aCoefR = D_TERM(pos[0]*,pos[1]*,pos[2]);
          // constant-coefficient
          //aCoefR = 1.0;
        }EndFor;
      } // end loop over grids
    }
  else if (a_params.CCoeftype == 1) // CCoef is constant 1
    {
      DataIterator dit = a_cCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& cCoef = a_cCoef[dit];
          ForAllXBNN(Real, cCoef, cCoef.box(), 0, cCoef.nComp());
          {
            cCoefR = 1.0; // constant 
          }
        }EndFor;
    } // end ACoeftype
}


#endif
