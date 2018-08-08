#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PoissonParameters.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "SetLevelData.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "HamiltonianPoissonOperatorFactory.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

void set_initial_psi(LevelData<FArrayBox> &a_psi,
                     const RealVect &a_dx, const PoissonParameters &a_params) {

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

void set_initial_phi(LevelData<FArrayBox> &a_phi, const RealVect &a_dx,
                     const PoissonParameters &a_params) {

  CH_assert(a_phi.nComp() == 1);

  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_phi_box = a_phi[dit()];
    Box b = this_phi_box.box();
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      for (int comp = 0; comp < a_phi.nComp(); ++comp) {
        this_phi_box(iv, comp) = a_params.initial_phi;
      }
    }
  }
} // end set_initial_phi

void set_rhs(LevelData<FArrayBox> &a_rhs, LevelData<FArrayBox> &a_psi,
             LevelData<FArrayBox> &a_phi, const RealVect &a_dx,
             const PoissonParameters &a_params) {

  CH_assert(a_rhs.nComp() == 1);
  int comp_number = 0;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_rhs = a_rhs[dit()];
    this_rhs.setVal(0, comp_number);
    Box this_box = this_rhs.box();

    // first set up the scalar field config
    FArrayBox &this_phi = a_phi[dit()];
    this_phi.setVal(0, comp_number);

    // rhs is cell-centered...
    RealVect cc_offset = 0.5 * a_dx * RealVect::Unit;

    // put in a gaussian for now
    Real strength = 1.0;
    Real scale = 1.0e-2; // variance
    RealVect center = 0.5 * RealVect::Unit;

    BoxIterator bit(this_box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv);
      loc *= a_dx;
      loc += cc_offset;

      RealVect dist = loc - center;
      Real radSqr =
            D_TERM(dist[0] * dist[0], +dist[1] * dist[1], +dist[2] * dist[2]);

      Real val = strength * exp(-radSqr / scale);
      this_phi(iv, 0) += val;
      this_rhs(iv, 0) += val;
    }
  }
} // end set_rhs

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox>& a_psi, LevelData<FArrayBox>& a_dpsi)
{
    DataIterator dit = a_psi.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox& psi_box  = a_psi[dit];
        FArrayBox& dpsi_box = a_dpsi[dit];
        psi_box += dpsi_box;
    }
}

// m(K, rho) = 2/3K^2 - 16piG rho
// for now we put rho to be a gaussian with center (0.5,0.5,0.5)
void set_m_value(Real &m, RealVect &loc, const PoissonParameters &a_params) {

  Real rho;
  RealVect center;
  center = RealVect(
      D_DECL(a_params.rho_center_1[0], a_params.rho_center_1[1],
             a_params.rho_center_1[2])); // put it in the middle for now

  RealVect dist = loc - center;
  Real radSqr =
      D_TERM(dist[0] * dist[0], +dist[1] * dist[1], +dist[2] * dist[2]);

  rho = a_params.rho_strength * exp(-radSqr / a_params.rho_scale);

  m = (2.0 / 3.0) * (a_params.constant_K * a_params.constant_K) -
          2.0 * (a_params.kappa_sq) * rho;
}

void
set_a_coef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_psi,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // cell centered
        RealVect cc_offset = 0.5*a_dx*RealVect::Unit; 

        FArrayBox& aCoef = a_aCoef[dit];
        FArrayBox& psi = a_psi[dit];
        Box thisBox = aCoef.box();
     
        BoxIterator bit(thisBox);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv);
            loc *= a_dx;
            loc += cc_offset;

            // m(K, phi) = 2/3 K^2 - 16 pi G rho
            Real m;
            set_m_value(m, loc, a_params);

            Real psi_0 = psi(iv,0);

            aCoef(iv,0) = (-0.625) * m * (psi_0 * psi_0 * psi_0 * psi_0);
        }
    }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
void
set_b_coef(LevelData<FArrayBox>& a_bCoef,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
  int comp_number = 0;
  DataIterator dit = a_bCoef.dataIterator();
  for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &this_bCoef = a_bCoef[dit()];
    this_bCoef.setVal(1.0, comp_number);
  }
}
