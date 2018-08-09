#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "HamiltonianPoissonOperator.H"
#include "AMRMultiGrid.H"
#include "AMRPoissonOpF_F.H"
#include "AverageF_F.H"
#include "BoxIterator.H"
#include "CoarseAverage.H"
#include "CoarseAverageFace.H"
#include "DebugOut.H"
#include "FORT_PROTO.H"
#include "FineInterp.H"
#include "HamiltonianPoissonOperatorF_F.H"
#include "InterpF_F.H"
#include "LayoutIterator.H"
#include "Misc.H"

#include "NamespaceHeader.H"

void HamiltonianPoissonOperator::residualI(LevelData<FArrayBox> &a_lhs,
                                           const LevelData<FArrayBox> &a_dpsi,
                                           const LevelData<FArrayBox> &a_rhs,
                                           bool a_homogeneous) {
  CH_TIME("HamiltonianPoissonOperator::residualI");

  LevelData<FArrayBox> &dpsi = (LevelData<FArrayBox> &)a_dpsi;
  Real dx = m_dx;
  const DisjointBoxLayout &dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = dpsi.dataIterator();
  {
    CH_TIME("HamiltonianPoissonOperator::residualIBC");

    for (dit.begin(); dit.ok(); ++dit) {
      m_bc(dpsi[dit], dbl[dit()], m_domain, dx, a_homogeneous);
    }
  }

  dpsi.exchange(dpsi.interval(), m_exchangeCopier);

  for (dit.begin(); dit.ok(); ++dit) {
    const Box &region = dbl[dit()];

#if CH_SPACEDIM == 1
    FORT_VCCOMPUTERES1D
#elif CH_SPACEDIM == 2
    FORT_VCCOMPUTERES2D
#elif CH_SPACEDIM == 3
    FORT_VCCOMPUTERES3D
#else
    This_will_not_compile !
#endif
        (CHF_FRA(a_lhs[dit]), CHF_CONST_FRA(dpsi[dit]),
         CHF_CONST_FRA(a_rhs[dit]), CHF_CONST_REAL(m_alpha),
         CHF_CONST_FRA((*m_aCoef)[dit]), CHF_CONST_REAL(m_beta),
         CHF_CONST_FRA((*m_bCoef)[dit]), CHF_BOX(region), CHF_CONST_REAL(m_dx));
  } // end loop over boxes
}

/**************************/
// this preconditioner first initializes dpsihat to (IA)dpsihat = rhshat
// (diagonization of L -- A is the matrix version of L)
// then smooths with a couple of passes of levelGSRB
void HamiltonianPoissonOperator::preCond(LevelData<FArrayBox> &a_dpsi,
                                         const LevelData<FArrayBox> &a_rhs) {
  CH_TIME("HamiltonianPoissonOperator::preCond");

  // diagonal term of this operator in:
  //
  //       alpha * a(i)
  //     + beta  * sum_over_dir (b(i-1/2*e_dir) + b(i+1/2*e_dir)) / (dx*dx)
  //
  // The inverse of this is our initial multiplier.

  int ncomp = a_dpsi.nComp();

  CH_assert(m_lambda.isDefined());
  CH_assert(a_rhs.nComp() == ncomp);
  CH_assert(m_bCoef->nComp() == ncomp);

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  // don't need to use a Copier -- plain copy will do
  DataIterator dit = a_dpsi.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    // also need to average and sum face-centered bCoefs to cell-centers
    Box gridBox = a_rhs[dit].box();

    // approximate inverse
    a_dpsi[dit].copy(a_rhs[dit]);
    a_dpsi[dit].mult(m_lambda[dit], gridBox, 0, 0, ncomp);
  }

  relax(a_dpsi, a_rhs, 2);
}

void HamiltonianPoissonOperator::applyOpI(LevelData<FArrayBox> &a_lhs,
                                          const LevelData<FArrayBox> &a_dpsi,
                                          bool a_homogeneous) {
  CH_TIME("HamiltonianPoissonOperator::applyOpI");
  LevelData<FArrayBox> &dpsi = (LevelData<FArrayBox> &)a_dpsi;
  Real dx = m_dx;
  const DisjointBoxLayout &dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = dpsi.dataIterator();

  for (dit.begin(); dit.ok(); ++dit) {
    m_bc(dpsi[dit], dbl[dit()], m_domain, dx, a_homogeneous);
  }

  applyOpNoBoundary(a_lhs, a_dpsi);
}

void HamiltonianPoissonOperator::applyOpNoBoundary(
    LevelData<FArrayBox> &a_lhs, const LevelData<FArrayBox> &a_dpsi) {
  CH_TIME("HamiltonianPoissonOperator::applyOpNoBoundary");

  LevelData<FArrayBox> &dpsi = (LevelData<FArrayBox> &)a_dpsi;

  const DisjointBoxLayout &dbl = a_lhs.disjointBoxLayout();
  DataIterator dit = dpsi.dataIterator();

  dpsi.exchange(dpsi.interval(), m_exchangeCopier);

  for (dit.begin(); dit.ok(); ++dit) {
    const Box &region = dbl[dit()];

#if CH_SPACEDIM == 1
    FORT_VCCOMPUTEOP1D
#elif CH_SPACEDIM == 2
    FORT_VCCOMPUTEOP2D
#elif CH_SPACEDIM == 3
    FORT_VCCOMPUTEOP3D
#else
    This_will_not_compile !
#endif
        (CHF_FRA(a_lhs[dit]), CHF_CONST_FRA(dpsi[dit]), CHF_CONST_REAL(m_alpha),
         CHF_CONST_FRA((*m_aCoef)[dit]), CHF_CONST_REAL(m_beta),
         CHF_CONST_FRA((*m_bCoef)[dit]), CHF_BOX(region), CHF_CONST_REAL(m_dx));
  } // end loop over boxes
}

void HamiltonianPoissonOperator::restrictResidual(
    LevelData<FArrayBox> &a_resCoarse, LevelData<FArrayBox> &a_dpsiFine,
    const LevelData<FArrayBox> &a_rhsFine) {
  CH_TIME("HamiltonianPoissonOperator::restrictResidual");

  homogeneousCFInterp(a_dpsiFine);
  const DisjointBoxLayout &dblFine = a_dpsiFine.disjointBoxLayout();
  for (DataIterator dit = a_dpsiFine.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &dpsi = a_dpsiFine[dit];
    m_bc(dpsi, dblFine[dit()], m_domain, m_dx, true);
  }

  a_dpsiFine.exchange(a_dpsiFine.interval(), m_exchangeCopier);

  for (DataIterator dit = a_dpsiFine.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &dpsi = a_dpsiFine[dit];
    const FArrayBox &rhs = a_rhsFine[dit];
    FArrayBox &res = a_resCoarse[dit];

    const FArrayBox &thisACoef = (*m_aCoef)[dit];
    const FArrayBox &thisBCoef = (*m_bCoef)[dit];

    Box region = dblFine.get(dit());
    const IntVect &iv = region.smallEnd();
    IntVect civ = coarsen(iv, 2);

    res.setVal(0.0);

#if CH_SPACEDIM == 1
    FORT_RESTRICTRESVC1D
#elif CH_SPACEDIM == 2
    FORT_RESTRICTRESVC2D
#elif CH_SPACEDIM == 3
    FORT_RESTRICTRESVC3D
#else
    This_will_not_compile !
#endif
        (CHF_FRA_SHIFT(res, civ), CHF_CONST_FRA_SHIFT(dpsi, iv),
         CHF_CONST_FRA_SHIFT(rhs, iv), CHF_CONST_REAL(m_alpha),
         CHF_CONST_FRA_SHIFT(thisACoef, iv), CHF_CONST_REAL(m_beta),
         CHF_CONST_FRA_SHIFT(thisBCoef, iv), CHF_BOX_SHIFT(region, iv),
         CHF_CONST_REAL(m_dx));
  }
}

void HamiltonianPoissonOperator::setAlphaAndBeta(const Real &a_alpha,
                                                 const Real &a_beta) {
  m_alpha = a_alpha;
  m_beta = a_beta;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void HamiltonianPoissonOperator::setCoefs(
    const RefCountedPtr<LevelData<FArrayBox>> &a_aCoef,
    const RefCountedPtr<LevelData<FArrayBox>> &a_bCoef, const Real &a_alpha,
    const Real &a_beta) {

  m_alpha = a_alpha;
  m_beta = a_beta;

  m_aCoef = a_aCoef;
  m_bCoef = a_bCoef;

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;
}

void HamiltonianPoissonOperator::resetLambda() {

  if (m_lambdaNeedsResetting) {

    Real scale = 1.0 / (m_dx * m_dx);

    // Compute it box by box, point by point
    for (DataIterator dit = m_lambda.dataIterator(); dit.ok(); ++dit) {
      FArrayBox &lambdaFab = m_lambda[dit];
      const FArrayBox &aCoefFab = (*m_aCoef)[dit];
      const FArrayBox &bCoefFab = (*m_bCoef)[dit];
      const Box &curBox = lambdaFab.box();

      // Compute the diagonal term
      lambdaFab.copy(aCoefFab);
      lambdaFab.mult(m_alpha);

      // EUGENE: Add in the Laplacian term 6.0*m_beta/(m_dx*m_dx)
      lambdaFab.plus(2.0 * SpaceDim * m_beta / (m_dx * m_dx));

      // Take its reciprocal
      lambdaFab.invert(1.0);
    }

    // Lambda is reset.
    m_lambdaNeedsResetting = false;
  }
}

// Compute the reciprocal of the diagonal entry of the operator matrix
void HamiltonianPoissonOperator::computeLambda() {
  CH_TIME("HamiltonianPoissonOperator::computeLambda");

  CH_assert(!m_lambda.isDefined());

  // Define lambda
  m_lambda.define(m_aCoef->disjointBoxLayout(), m_aCoef->nComp());
  resetLambda();
}

// KC TODO: Don't need this any more?
// reflux operator
void HamiltonianPoissonOperator::reflux(
    const LevelData<FArrayBox> &a_dpsiFine, const LevelData<FArrayBox> &a_dpsi,
    LevelData<FArrayBox> &a_residual,
    AMRLevelOp<LevelData<FArrayBox>> *a_finerOp) {}
/*
  CH_TIMERS("HamiltonianPoissonOperator::reflux");

  m_levfluxreg.setToZero();
  Interval interv(0, a_dpsi.nComp() - 1);

  CH_TIMER("HamiltonianPoissonOperator::reflux::incrementCoarse", t2);
  CH_START(t2);

  DataIterator dit = a_dpsi.dataIterator();
  for (dit.reset(); dit.ok(); ++dit) {
    const FArrayBox &coarfab = a_dpsi[dit];
    const FluxBox &coarBCoef = (*m_bCoef)[dit];
    const Box &gridBox = a_dpsi.getBoxes()[dit];

    if (m_levfluxreg.hasCF(dit())) {
      for (int idir = 0; idir < SpaceDim; idir++) {
        FArrayBox coarflux;
        Box faceBox = surroundingNodes(gridBox, idir);

        getFlux(coarflux, coarfab, coarBCoef, faceBox, idir);

        Real scale = 1.0;
        m_levfluxreg.incrementCoarse(coarflux, scale, dit(), interv, interv,
                                     idir);
      }
    }
  }

  CH_STOP(t2);

  // const cast:  OK because we're changing ghost cells only
  LevelData<FArrayBox> &dpsiFineRef = (LevelData<FArrayBox> &)a_dpsiFine;

  HamiltonianPoissonOperator *finerAMRPOp =
      (HamiltonianPoissonOperator *)a_finerOp;
  QuadCFInterp &quadCFI = finerAMRPOp->m_interpWithCoarser;

  quadCFI.coarseFineInterp(dpsiFineRef, a_dpsi);
  // I'm pretty sure this is not necessary. bvs -- flux calculations use
  // outer ghost cells, but not inner ones
  // dpsiFineRef.exchange(a_dpsiFine.interval());
  IntVect dpsiGhost = dpsiFineRef.ghostVect();
  int ncomps = a_dpsiFine.nComp();

  CH_TIMER("HamiltonianPoissonOperator::reflux::incrementFine", t3);
  CH_START(t3);

  DataIterator ditf = a_dpsiFine.dataIterator();
  const DisjointBoxLayout &dblFine = a_dpsiFine.disjointBoxLayout();
  for (ditf.reset(); ditf.ok(); ++ditf) {
    const FArrayBox &dpsifFab = a_dpsiFine[ditf];
    const FluxBox &fineBCoef = (*(finerAMRPOp->m_bCoef))[ditf];
    const Box &gridbox = dblFine.get(ditf());

    for (int idir = 0; idir < SpaceDim; idir++) {
      // int normalGhost = dpsiGhost[idir];
      SideIterator sit;
      for (sit.begin(); sit.ok(); sit.next()) {
        if (m_levfluxreg.hasCF(ditf(), sit())) {
          Side::LoHiSide hiorlo = sit();
          Box fluxBox = bdryBox(gridbox, idir, hiorlo, 1);

          FArrayBox fineflux(fluxBox, ncomps);
          getFlux(fineflux, dpsifFab, fineBCoef, fluxBox, idir, m_refToFiner);

          Real scale = 1.0;
          m_levfluxreg.incrementFine(fineflux, scale, ditf(), interv, interv,
                                     idir, hiorlo);
        }
      }
    }
  }

  CH_STOP(t3);

  Real scale = 1.0 / m_dx;
  m_levfluxreg.reflux(a_residual, scale);
}
*/

void HamiltonianPoissonOperator::levelGSRB(LevelData<FArrayBox> &a_dpsi,
                                           const LevelData<FArrayBox> &a_rhs) {
  CH_TIME("HamiltonianPoissonOperator::levelGSRB");

  CH_assert(a_dpsi.isDefined());
  CH_assert(a_rhs.isDefined());
  CH_assert(a_dpsi.ghostVect() >= IntVect::Unit);
  CH_assert(a_dpsi.nComp() == a_rhs.nComp());

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  const DisjointBoxLayout &dbl = a_dpsi.disjointBoxLayout();

  DataIterator dit = a_dpsi.dataIterator();

  // do first red, then black passes
  for (int whichPass = 0; whichPass <= 1; whichPass++) {
    CH_TIMERS("HamiltonianPoissonOperator::levelGSRB::Compute");

    // fill in intersection of ghostcells and a_dpsi's boxes
    {
      CH_TIME("HamiltonianPoissonOperator::levelGSRB::homogeneousCFInterp");
      homogeneousCFInterp(a_dpsi);
    }

    {
      CH_TIME("HamiltonianPoissonOperator::levelGSRB::exchange");
      a_dpsi.exchange(a_dpsi.interval(), m_exchangeCopier);
    }

    {
      CH_TIME("HamiltonianPoissonOperator::levelGSRB::BCs");
      // now step through grids...
      for (dit.begin(); dit.ok(); ++dit) {
        // invoke physical BC's where necessary
        m_bc(a_dpsi[dit], dbl[dit()], m_domain, m_dx, true);
      }
    }

    for (dit.begin(); dit.ok(); ++dit) {
      const Box &region = dbl.get(dit());

#if CH_SPACEDIM == 1
      FORT_GSRBHELMHOLTZVC1D
#elif CH_SPACEDIM == 2
      FORT_GSRBHELMHOLTZVC2D
#elif CH_SPACEDIM == 3
      FORT_GSRBHELMHOLTZVC3D
#else
      This_will_not_compile !
#endif
          (CHF_FRA(a_dpsi[dit]), CHF_CONST_FRA(a_rhs[dit]), CHF_BOX(region),
           CHF_CONST_REAL(m_dx), CHF_CONST_REAL(m_alpha),
           CHF_CONST_FRA((*m_aCoef)[dit]), CHF_CONST_REAL(m_beta),
           CHF_CONST_FRA((*m_bCoef)[dit]), CHF_CONST_FRA(m_lambda[dit]),
           CHF_CONST_INT(whichPass));
    } // end loop through grids
  }   // end loop through red-black
}

void HamiltonianPoissonOperator::levelMultiColor(
    LevelData<FArrayBox> &a_dpsi, const LevelData<FArrayBox> &a_rhs) {
  CH_TIME("HamiltonianPoissonOperator::levelMultiColor");
  MayDay::Abort(
      "HamiltonianPoissonOperator::levelMultiColor - Not implemented");
}

void HamiltonianPoissonOperator::looseGSRB(LevelData<FArrayBox> &a_dpsi,
                                           const LevelData<FArrayBox> &a_rhs) {
  CH_TIME("HamiltonianPoissonOperator::looseGSRB");
  MayDay::Abort("HamiltonianPoissonOperator::looseGSRB - Not implemented");
}

void HamiltonianPoissonOperator::overlapGSRB(
    LevelData<FArrayBox> &a_dpsi, const LevelData<FArrayBox> &a_rhs) {
  CH_TIME("HamiltonianPoissonOperator::overlapGSRB");
  MayDay::Abort("HamiltonianPoissonOperator::overlapGSRB - Not implemented");
}

void HamiltonianPoissonOperator::levelGSRBLazy(
    LevelData<FArrayBox> &a_dpsi, const LevelData<FArrayBox> &a_rhs) {
  CH_TIME("HamiltonianPoissonOperator::levelGSRBLazy");
  MayDay::Abort("HamiltonianPoissonOperator::levelGSRBLazy - Not implemented");
}

void HamiltonianPoissonOperator::levelJacobi(
    LevelData<FArrayBox> &a_dpsi, const LevelData<FArrayBox> &a_rhs) {

  CH_TIME("HamiltonianPoissonOperator::levelJacobi");

  // Recompute the relaxation coefficient if needed.
  resetLambda();

  LevelData<FArrayBox> resid;
  create(resid, a_rhs);

  // Get the residual
  residual(resid, a_dpsi, a_rhs, true);

  // Multiply by the weights
  DataIterator dit = m_lambda.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    resid[dit].mult(m_lambda[dit]);
  }

  // Do the Jacobi relaxation
  incr(a_dpsi, resid, 0.5);

  // exchange ghost cells
  a_dpsi.exchange(a_dpsi.interval(), m_exchangeCopier);
}

// KC TODO: Think we don't need this any more
/*
void HamiltonianPoissonOperator::getFlux(FArrayBox &a_flux,
                                         const FArrayBox &a_data,
                                         const FArrayBox &b_data,
                                         const Box &a_facebox, int a_dir,
                                         int a_ref) const {

  CH_TIME("HamiltonianPoissonOperator::getFlux");

  CH_assert(a_dir >= 0);
  CH_assert(a_dir < SpaceDim);
  CH_assert(!a_data.box().isEmpty());
  CH_assert(!a_facebox.isEmpty());

  // probably the simplest way to test centering
  // a_box needs to be face-centered in the a_dir
  Box faceTestBox(IntVect::Zero, IntVect::Unit);
  faceTestBox.surroundingNodes(a_dir);
  CH_assert(a_facebox.type() == faceTestBox.type());

  const FArrayBox &bCoefDir = a_bCoef[a_dir];

  // reality check for bCoef
  CH_assert(bCoefDir.box().contains(a_facebox));

  a_flux.resize(a_facebox, a_data.nComp());
  BoxIterator bit(a_facebox);

  Real scale = m_beta * a_ref / m_dx;

  for (bit.begin(); bit.ok(); bit.next()) {
    IntVect iv = bit();
    IntVect shiftiv = BASISV(a_dir);
    IntVect ivlo = iv - shiftiv;
    IntVect ivhi = iv;

    CH_assert(a_data.box().contains(ivlo));
    CH_assert(a_data.box().contains(ivhi));

    for (int ivar = 0; ivar < a_data.nComp(); ivar++) {
      Real dpsihi = a_data(ivhi, ivar);
      Real dpsilo = a_data(ivlo, ivar);
      Real graddpsi = (dpsihi - dpsilo) * scale;

      a_flux(iv, ivar) = -bCoefDir(iv, ivar) * graddpsi;
    }
  }
}
*/

//-----------------------------------------------------------------------
void HamiltonianPoissonOperator::setTime(Real a_time) {
  // Jot down the time.
  m_time = a_time;

  // Interpolate the b coefficient data if necessary / possible. If
  // the B coefficient depends upon the solution, the operator is nonlinear
  // and the integrator must decide how to treat it.
  if (!m_bCoefInterpolator.isNull() &&
      !m_bCoefInterpolator->dependsUponSolution())
    m_bCoefInterpolator->interpolate(*m_bCoef, a_time);

  // Our relaxation parameter is officially out of date!
  m_lambdaNeedsResetting = true;

  // Set the time on the boundary holder.
  m_bc.setTime(a_time);

  // Notify our observers that the time has been set.
  // FIXME: Must implement response of multigrid operators!
  notifyObserversOfChange();
}

#include "NamespaceFooter.H"
