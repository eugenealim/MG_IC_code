#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GRIDFUNCS_H_
#define _GRIDFUNCS_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "HamiltonianPoissonOperator.H"
#include "HamiltonianPoissonOperatorFactory.H"
#include "ParmParse.H"
#include "PoissonParameters.H"
#include "FArrayBox.H"

#include "UsingNamespace.H"

// functions operating on grids
extern int setGrids(Vector<DisjointBoxLayout>& vectGrids,
                    PoissonParameters&         a_params);

extern void setRHS(LevelData<FArrayBox>&    a_rhs,
                   LevelData<FArrayBox>&    a_chi,
                   LevelData<FArrayBox>&    a_phi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_initial_chi(LevelData<FArrayBox>&    a_chi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_initial_phi(LevelData<FArrayBox>&    a_phi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

/* see BoxTools/BaseFabMacros.H  for ForAllXBNN loop macro */
extern void
setACoef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx);

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
extern void
setBCoef(LevelData<FArrayBox>& a_bCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx);

extern Real calculate_M_value(RealVect&          loc,
                    const PoissonParameters& a_params);

extern void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         PoissonParameters&         a_params);

extern void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels) ;

#endif
