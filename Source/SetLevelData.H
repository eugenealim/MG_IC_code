/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _SETLEVELDATA_H_
#define _SETLEVELDATA_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "VariableCoeffPoissonOperator.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "ParmParse.H"
#include "PoissonParameters.H"
#include "FArrayBox.H"
#include "UserVariables.hpp"

#include "UsingNamespace.H"

// various functions operating on grids to set up data for initial
// conditions, spatially dependent coefficients and rhs etc

extern void set_rhs(LevelData<FArrayBox>&    a_rhs,
                   LevelData<FArrayBox>&    a_psi,
                   LevelData<FArrayBox>&    a_phi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_initial_psi(LevelData<FArrayBox>&    a_psi,
                   LevelData<FArrayBox>&    a_dpsi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_initial_phi(LevelData<FArrayBox>&    a_phi,
                   const RealVect&          a_dx,
                   const PoissonParameters& a_params);

extern void set_m_value(Real& m, const Real& phi_here,
                    const PoissonParameters& a_params);

extern void set_A2_value(Real &A2, const PoissonParameters &a_params);

extern void
set_a_coef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_psi,
         LevelData<FArrayBox>& a_phi,
         const PoissonParameters& a_params,
         const RealVect& a_dx);

extern void
set_b_coef(LevelData<FArrayBox>& a_bCoef,
         const PoissonParameters& a_params,
         const RealVect& a_dx);

extern void set_update_psi0(LevelData<FArrayBox>& a_psi, LevelData<FArrayBox>& a_dpsi,
                            const Copier& a_exchange_copier);

extern void set_output_data(LevelData<FArrayBox>& a_vars, LevelData<FArrayBox>& a_psi, 
                     LevelData<FArrayBox>& a_phi, const PoissonParameters &a_params);

#endif
