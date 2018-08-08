#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SETLEVELDATA_H_
#define _SETLEVELDATA_H_

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

// KC TODO - make this a namespace SetLevelData?
// various functions operating on grids
extern void set_rhs(LevelData<FArrayBox>&    a_rhs,
                   LevelData<FArrayBox>&    a_psi,
                   LevelData<FArrayBox>&    a_phi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_initial_psi(LevelData<FArrayBox>&    a_psi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_initial_phi(LevelData<FArrayBox>&    a_phi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_m_value(Real& m, RealVect& loc,
                    const PoissonParameters& a_params);

/* see BoxTools/BaseFabMacros.H  for ForAllXBNN loop macro */
extern void
set_a_coef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_psi,
         const PoissonParameters& a_params,
         const RealVect& a_dx);

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
extern void
set_b_coef(LevelData<FArrayBox>& a_bCoef,
         const PoissonParameters& a_params,
         const RealVect& a_dx);

extern void set_update_psi0(LevelData<FArrayBox>& a_psi, LevelData<FArrayBox>& a_dpsi);

#endif
