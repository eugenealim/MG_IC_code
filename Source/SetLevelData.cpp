/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SetLevelData.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "PoissonParameters.H"
#include "SetLevelDataF_F.H"
#include "UsingNamespace.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij, scalar field phi
// and \bar Aij = psi^2 A_ij
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars, LevelData<FArrayBox> &a_dpsi, 
                     const RealVect &a_dx,
                     const PoissonParameters &a_params) {

  RealVect cc_offset = 0.5 * a_dx * RealVect::Unit;

  DataIterator dit = a_multigrid_vars.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
    FArrayBox &dpsi_box = a_dpsi[dit()];
    Box b = multigrid_vars_box.box();
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {

      // work out location on the grid
      IntVect iv = bit();
      RealVect loc(iv);
      loc *= a_dx;
      loc += cc_offset;
      RealVect center =
          RealVect(D_DECL(a_params.phi_center[0], a_params.phi_center[1],
                          a_params.phi_center[2]));
      loc -= center;

      // distance from centre squared
      Real r2 = D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]);

      // a gaussian
      multigrid_vars_box(iv, c_phi_0) =
          a_params.phi_strength * exp(-r2 / a_params.phi_scale);

      multigrid_vars_box(iv, c_psi) = a_params.initial_psi;
      dpsi_box(iv, 0) = 0.0;
    }
  }
} // end set_initial_conditions

void set_rhs(LevelData<FArrayBox> &a_rhs, LevelData<FArrayBox> &a_multigrid_vars,
             const RealVect &a_dx,
             const PoissonParameters &a_params) {

  DataIterator dit = a_multigrid_vars.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit];
    FArrayBox &rhs_box = a_rhs[dit()];
    rhs_box.setVal(0.0, 0);
    Box this_box = rhs_box.box(); // no ghost cells

    // calculate the laplacian of psi across the box
    FArrayBox laplacian_of_psi(this_box, 1);
    FORT_GETLAPLACIANPSIF(CHF_FRA1(laplacian_of_psi, 0),
                          CHF_CONST_FRA1(multigrid_vars_box, c_psi),
                          CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

    // calculate the rho contribution from gradients of phi
    FArrayBox rho_gradient(this_box, 1);
    FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, 0),
                        CHF_CONST_FRA1(multigrid_vars_box, c_phi_0),
                        CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();

      // rhs = m/8 psi_0^5 - 2 pi rho_grad psi_0  - laplacian(psi_0)
      Real m = 0;
      set_m_value(m, multigrid_vars_box(iv, c_phi), a_params);

      // Also \bar  A_ij \bar A^ij
      Real A2 = 0.0;
      A2 =      multigrid_vars_box(iv, c_A11_0) * multigrid_vars_box(iv, c_A11_0)
           +  2*multigrid_vars_box(iv, c_A12_0) * multigrid_vars_box(iv, c_A12_0)
           +  2*multigrid_vars_box(iv, c_A13_0) * multigrid_vars_box(iv, c_A13_0)
           +    multigrid_vars_box(iv, c_A22_0) * multigrid_vars_box(iv, c_A22_0)
           +  2*multigrid_vars_box(iv, c_A23_0) * multigrid_vars_box(iv, c_A23_0)
           +    multigrid_vars_box(iv, c_A33_0) * multigrid_vars_box(iv, c_A33_0);

      Real psi_0 = multigrid_vars_box(iv, c_psi);
      rhs_box(iv, 0) =
          0.125 * m * pow(psi_0, 5.0) - 0.125 * A2 * pow(psi_0, -7.0) -
          2.0 * M_PI * a_params.G_Newton * rho_gradient(iv) * psi_0 -
          laplacian_of_psi(iv);
    }
  }
} // end set_rhs

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox> &a_multigrid_vars, LevelData<FArrayBox> &a_dpsi,
                     const Copier &a_exchange_copier) {

  // first exchange ghost cells for dpsi so they are filled with the correct
  // values
  a_dpsi.exchange(a_dpsi.interval(), a_exchange_copier);

  DataIterator dit = a_multigrid_vars.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit];
    FArrayBox &dpsi_box = a_dpsi[dit];

    Box this_box = multigrid_vars_box.box();
    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      multigrid_vars_box(iv, c_psi) += dpsi_box(iv, 0);
    }
  }
}

// m(K, rho) = 2/3K^2 - 16piG rho
void set_m_value(Real &m, const Real &phi_here,
                 const PoissonParameters &a_params) {

  // KC TODO:
  // For now rho is just the gradient term which is kept separate
  // ... may want to add V(phi) and phidot/Pi here later though
  Real Pi_field = 0.0;
  Real V_of_phi = 0.0;
  Real rho = 0.5*Pi_field*Pi_field + V_of_phi; 

  m = (2.0 / 3.0) * (a_params.constant_K * a_params.constant_K) -
      16.0 * M_PI * a_params.G_Newton * rho;
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef, 
                LevelData<FArrayBox> &a_multigrid_vars,
                const PoissonParameters &a_params, const RealVect &a_dx) {

  DataIterator dit = a_aCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &aCoef_box = a_aCoef[dit];
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit];
    Box this_box = aCoef_box.box();

    // calculate the rho contribution from gradients of phi
    FArrayBox rho_gradient(this_box, 1);
    FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, 0),
                        CHF_CONST_FRA1(multigrid_vars_box, c_phi_0),
                        CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      // m(K, phi) = 2/3 K^2 - 16 pi G rho
      Real m;
      set_m_value(m, multigrid_vars_box(iv, c_phi_0), a_params);

      // Also \bar  A_ij \bar A^ij
      Real A2 = 0.0;
      A2 =      multigrid_vars_box(iv, c_A11_0) * multigrid_vars_box(iv, c_A11_0)
           +  2*multigrid_vars_box(iv, c_A12_0) * multigrid_vars_box(iv, c_A12_0)
           +  2*multigrid_vars_box(iv, c_A13_0) * multigrid_vars_box(iv, c_A13_0)
           +    multigrid_vars_box(iv, c_A22_0) * multigrid_vars_box(iv, c_A22_0)
           +  2*multigrid_vars_box(iv, c_A23_0) * multigrid_vars_box(iv, c_A23_0)
           +    multigrid_vars_box(iv, c_A33_0) * multigrid_vars_box(iv, c_A33_0);

      Real psi_0 = multigrid_vars_box(iv, c_psi);
      aCoef_box(iv, 0) = -0.625 * m * pow(psi_0, 4.0) - A2 * pow(psi_0, -8.0) +
                     2.0 * M_PI * a_params.G_Newton * rho_gradient(iv);
    }
  }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
  // the rhs source of the Poisson eqn
void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                const PoissonParameters &a_params, const RealVect &a_dx) {

  CH_assert(a_bCoef.nComp() == 1);
  int comp_number = 0;

  for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &bCoef_box = a_bCoef[dit()];
    bCoef_box.setVal(1.0, comp_number);
  }
}

// used to set output data for all ADM Vars for GRChombo restart
void set_output_data(LevelData<FArrayBox> &a_grchombo_vars,
                     LevelData<FArrayBox> &a_multigrid_vars,
                     const PoissonParameters &a_params) {

  CH_assert(a_grchombo_vars.nComp() == NUM_GRCHOMBO_VARS);
  CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

  DataIterator dit = a_grchombo_vars.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &grchombo_vars_box = a_grchombo_vars[dit];
    FArrayBox &multigrid_vars_box = a_multigrid_vars[dit];

    // first set everything to zero
    for (int comp = 0; comp < NUM_GRCHOMBO_VARS; comp++) {
      grchombo_vars_box.setVal(0.0, comp);
    }

    // now set non zero terms - const across whole box
    // Conformally flat, and lapse = 1
    grchombo_vars_box.setVal(1.0, c_h11);
    grchombo_vars_box.setVal(1.0, c_h22);
    grchombo_vars_box.setVal(1.0, c_h33);
    grchombo_vars_box.setVal(1.0, c_lapse);

    // constant K
    grchombo_vars_box.setVal(a_params.constant_K, c_K);

    // now non constant terms by location
    Box this_box = grchombo_vars_box.box();
    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      // Copy phi across
      grchombo_vars_box(iv, c_phi) = multigrid_vars_box(iv, c_phi_0);

      // GRChombo conformal factor chi = psi^-4
      grchombo_vars_box(iv, c_chi) = pow(multigrid_vars_box(iv, c_psi), -4.0);
    }
  }
}
