#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SetLevelData.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "HamiltonianPoissonOperatorFactory.H"
#include "LoadBalance.H"
#include "PoissonParameters.H"
#include "SetLevelDataF_F.H"
#include "UsingNamespace.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij
// Usually just psi = 1 for flat space but may want to change this for BHs
// e.g. to put in schwarzschild
void set_initial_psi(LevelData<FArrayBox> &a_psi, const RealVect &a_dx,
                     const PoissonParameters &a_params) {

  CH_assert(a_psi.nComp() == 1);

  for (DataIterator dit = a_psi.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_psi_box = a_psi[dit()];
    Box b = this_psi_box.box();
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      for (int comp = 0; comp < a_psi.nComp(); ++comp) {
        this_psi_box(iv, comp) = a_params.initial_psi;
      }
    }
  }
} // end set_initial_psi

// set initial value for the scalar field phi
void set_initial_phi(LevelData<FArrayBox> &a_phi, const RealVect &a_dx,
                     const PoissonParameters &a_params) {

  CH_assert(a_phi.nComp() == 1);
  int comp_num = 1;
  RealVect cc_offset = 0.5 * a_dx * RealVect::Unit;

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_phi_box = a_phi[dit()];
    Box b = this_phi_box.box();
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {

      // work out location on the grid
      IntVect iv = bit();
      RealVect loc(iv);
      loc *= a_dx;
      loc += cc_offset;
      RealVect center =
          RealVect(D_DECL(a_params.rho_center_1[0], a_params.rho_center_1[1],
                          a_params.rho_center_1[2]));
      loc -= center;

      // distance from centre squared
      Real r2 = D_TERM(loc[0] * loc[0], +loc[1] * loc[1], +loc[2] * loc[2]);

      // a gaussian
      this_phi_box(iv, comp_num) =
          a_params.rho_strength * exp(-r2 / a_params.rho_scale);
    }
  }
} // end set_initial_phi

void set_rhs(LevelData<FArrayBox> &a_rhs, LevelData<FArrayBox> &a_psi,
             LevelData<FArrayBox> &a_phi, const RealVect &a_dx,
             const PoissonParameters &a_params) {

  CH_assert(a_rhs.nComp() == 1);
  int comp_number = 0;

  // rhs is cell-centered
  RealVect cc_offset = 0.5 * a_dx * RealVect::Unit;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_rhs = a_rhs[dit()];
    this_rhs.setVal(0, comp_number);
    Box this_box = this_rhs.box();

    // first get the scalar field and psi config
    FArrayBox &this_phi = a_phi[dit()];
    FArrayBox &this_psi = a_psi[dit()];
    BoxIterator bit(this_box);

    // calculate the laplacian of psi across the box
    FArrayBox laplacian_of_psi(this_box, 1);;
    FORT_GETLAPLACIANPSIF(CHF_FRA1(laplacian_of_psi, comp_number),
                          CHF_CONST_FRA1(this_psi, comp_number),
                          CHF_CONST_REAL(a_dx[0]),
                          CHF_BOX(this_box));

    // calculate the rho contribution from gradients of phi
    FArrayBox rho_gradient(this_box, 1);;
    FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, comp_number),
                          CHF_CONST_FRA1(this_phi, comp_number),
                          CHF_CONST_REAL(a_dx[0]),
                          CHF_BOX(this_box));

    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();

      // rhs = m/8 psi_0^5 - laplacian(psi_0)
      Real m = 0;
      set_m_value(m, this_phi(iv, comp_number), a_params);
      Real psi_0 = this_psi(iv, comp_number);
      this_rhs(iv, comp_number) =   0.125 * m * pow(psi_0, 5.0)
                                  + 0.25 * a_params.kappa_sq * rho_gradient(iv) * psi_0
                                  - laplacian_of_psi(iv);
    }
  }
} // end set_rhs

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox> &a_psi,
                     LevelData<FArrayBox> &a_dpsi) {
  DataIterator dit = a_psi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    FArrayBox &psi_box = a_psi[dit];
    FArrayBox &dpsi_box = a_dpsi[dit];
    psi_box += dpsi_box;
  }
}

// m(K, rho) = 2/3K^2 - 16piG rho
void set_m_value(Real &m, Real &phi, const PoissonParameters &a_params) {

  // KC TODO:
  // For now rho is just the gradient term which is kept separate
  // ... may want to add V(phi) and phidot here though
  Real rho = 0.0;

  m = (2.0 / 3.0) * (a_params.constant_K * a_params.constant_K) -
      2.0 * (a_params.kappa_sq) * rho;
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef, LevelData<FArrayBox> &a_psi,
                LevelData<FArrayBox> &a_phi, const PoissonParameters &a_params,
                const RealVect &a_dx)
{
  CH_assert(a_phi.nComp() == 1);
  int comp_number = 0;

  DataIterator dit = a_aCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    // cell centered
    RealVect cc_offset = 0.5 * a_dx * RealVect::Unit;

    FArrayBox &aCoef = a_aCoef[dit];
    FArrayBox &this_psi = a_psi[dit];
    FArrayBox &this_phi = a_phi[dit];
    Box this_box = aCoef.box();

    // calculate the rho contribution from gradients of phi
    FArrayBox rho_gradient(this_box, 1);;
    FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, comp_number),
                          CHF_CONST_FRA1(this_phi, comp_number),
                          CHF_CONST_REAL(a_dx[0]),
                          CHF_BOX(this_box));

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv);
      loc *= a_dx;
      loc += cc_offset;

      // m(K, phi) = 2/3 K^2 - 16 pi G rho
      Real m;
      set_m_value(m, this_phi(iv, comp_number), a_params);
      Real psi_0 = this_psi(iv, 0);
      aCoef(iv, 0) = - 0.625 * m * pow(psi_0, 4.0) 
                     - 0.25 * a_params.kappa_sq * rho_gradient(iv);
    }
  }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                const PoissonParameters &a_params, const RealVect &a_dx) {
  int comp_number = 0;
  DataIterator dit = a_bCoef.dataIterator();
  for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_bCoef = a_bCoef[dit()];
    this_bCoef.setVal(1.0, comp_number);
  }
}