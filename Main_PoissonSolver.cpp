#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AMRIO.H"
#include "BClocalFuncs.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "DebugDump.H"
#include "FABView.H"
#include "FArrayBox.H"
#include "HamiltonianPoissonOperatorFactory.H"
#include "LevelData.H"
#include "LoadBalance.H"
#include "MultilevelLinearOp.H"
#include "ParmParse.H"
#include "PoissonParameters.H"
#include "SetGrids.H"
#include "SetLevelData.H"
#include "UsingNamespace.H"
#include "WriteOutput.H"
#include "computeNorm.H"
#include <iostream>

#ifdef CH_Linux
// Should be undefined by default
//#define TRAP_FPE
#undef TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif

using std::cerr;

// Sets up and runs the solver
// The equation solved is: [aCoef*I + bCoef*Laplacian](dpsi) = rhs
// We assume conformal flatness, K=const and Momentum constraint satisfied
// trivially  lapse = 1 shift = 0, phi is the scalar field and is used to
// calculate the rhs
int poissonSolve(Vector<LevelData<FArrayBox> *> &a_dpsi,
                 Vector<LevelData<FArrayBox> *> &a_psi,
                 Vector<LevelData<FArrayBox> *> &a_phi,
                 Vector<LevelData<FArrayBox> *> &a_rhs,
                 const Vector<DisjointBoxLayout> &a_grids,
                 const PoissonParameters &a_params) {
  ParmParse pp;

  int nlevels = a_params.numLevels;
  Vector<RefCountedPtr<LevelData<FArrayBox>>> aCoef(nlevels);
  Vector<RefCountedPtr<LevelData<FArrayBox>>> bCoef(nlevels);
  Vector<ProblemDomain> vectDomains(nlevels);
  Vector<RealVect> vectDx(nlevels);

  RealVect dxLev = RealVect::Unit;
  dxLev *= a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);

  // Declare variables here, with num comps = 1 and ghosts for sources
  // want output data to have 3 ghost cells to match GRChombo, although
  // not currently needed for 2nd order stencils used here
  for (int ilev = 0; ilev < nlevels; ilev++) {
    IntVect ghosts = 3 * IntVect::Unit;
    a_rhs[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero);
    a_dpsi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1, ghosts);
    a_psi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1, ghosts);
    a_phi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1, ghosts);
    aCoef[ilev] = RefCountedPtr<LevelData<FArrayBox>>(
        new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
    bCoef[ilev] = RefCountedPtr<LevelData<FArrayBox>>(
        new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
    vectDomains[ilev] = domLev;
    vectDx[ilev] = dxLev;

    // set initial guess for psi and zero the difference
    set_initial_psi(*a_psi[ilev], *a_dpsi[ilev], vectDx[ilev], a_params);

    // set source - scalar field
    set_initial_phi(*a_phi[ilev], vectDx[ilev], a_params);

    // prepare dx, domain for next level
    dxLev /= a_params.refRatio[ilev];
    domLev.refine(a_params.refRatio[ilev]);
  }

  // set up linear operator
  int lBase = 0;
  MultilevelLinearOp<FArrayBox> mlOp;

  int numMGIter = 1;
  pp.query("numMGIterations", numMGIter);
  mlOp.m_num_mg_iterations = numMGIter;

  int numMGSmooth = 4;
  pp.query("numMGsmooth", numMGSmooth);
  mlOp.m_num_mg_smooth = numMGSmooth;

  int preCondSolverDepth = -1;
  pp.query("preCondSolverDepth", preCondSolverDepth);
  mlOp.m_preCondSolverDepth = preCondSolverDepth;

  Real tolerance = 1.0e-7;
  pp.query("tolerance", tolerance);

  int max_iter = 10;
  pp.query("max_iterations", max_iter);

  BiCGStabSolver<Vector<LevelData<FArrayBox> *>> solver; // define solver object

  // Main loop
  int NL_iter = 0;
  int max_NL_iter = 1;
  pp.query("max_NL_iterations", max_NL_iter);
  Real dpsi_norm = 0.0;

  for (NL_iter; NL_iter < max_NL_iter; NL_iter++) {

    pout() << "Main Loop Iteration " << (NL_iter + 1) << " out of "
           << max_NL_iter << endl;

    // Assign values for coefficients here
    for (int ilev = 0; ilev < nlevels; ilev++) {
      set_a_coef(*aCoef[ilev], *a_psi[ilev], *a_phi[ilev], a_params,
                 vectDx[ilev]);
      set_b_coef(*bCoef[ilev], a_params, vectDx[ilev]);
      set_rhs(*a_rhs[ilev], *a_psi[ilev], *a_phi[ilev], vectDx[ilev], a_params);
    }

    // set up solver
    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox>>> opFactory =
        RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox>>>(
            defineOperatorFactory(a_grids, vectDomains, aCoef, bCoef,
                                  a_params));

    mlOp.define(a_grids, a_params.refRatio, vectDomains, vectDx, opFactory,
                lBase);

    bool homogeneousBC = false;
    solver.define(&mlOp, homogeneousBC);
    solver.m_verbosity = a_params.verbosity;
    solver.m_normType = 0;
    solver.m_eps = tolerance;
    solver.m_imax = max_iter;

    output_solver_data(a_dpsi, a_psi, a_phi, a_rhs, a_grids, a_params, NL_iter);

    solver.solve(a_dpsi, a_rhs);

    // Add the solution to the linearised eqn to the previous iteration
    // ie psi -> psi + dpsi
    // need to fill interlevel and intralevel ghosts first in dpsi
    for (int ilev = 0; ilev < nlevels; ilev++) {

      // For interlevel ghosts
      if (ilev > 0) {
        int num_comps = 1;
        QuadCFInterp quadCFI(a_grids[ilev], &a_grids[ilev - 1], vectDx[ilev][0],
                             a_params.refRatio[ilev], num_comps,
                             vectDomains[ilev]);
        quadCFI.coarseFineInterp(*a_dpsi[ilev], *a_dpsi[ilev - 1]);
      }

      // For intralevel ghosts
      Copier exchange_copier;
      exchange_copier.exchangeDefine(a_grids[ilev], IntVect::Unit);

      // now the update
      set_update_psi0(*a_psi[ilev], *a_dpsi[ilev], exchange_copier);
    }

    // check if converged and if so exit NL iteration for loop
    dpsi_norm = computeNorm(a_dpsi, a_params.refRatio, a_params.coarsestDx,
                            Interval(0, 0));
    pout() << "The norm of dpsi at step " << NL_iter << " is " << dpsi_norm
           << endl;

    if (dpsi_norm < tolerance) {
      break;
    }

  } // end NL iteration loop

  pout() << "The norm of dpsi at the final step was " << dpsi_norm << endl;

  // Mayday if result not converged at all - using a fairly generous threshold
  // for this as usually non convergence means everything goes nuts
  if (dpsi_norm > 1e-4) {
    MayDay::Error(
        "NL iterations did not converge - may need a better initial guess");
  }

  // now output final data in a form which can be read as a checkpoint file
  // for the AMR time dependent runs
  output_final_data(a_psi, a_phi, a_grids, vectDx, vectDomains, a_params);

  int exitStatus = solver.m_exitStatus;
  // note that for AMRMultiGrid, success = 1.
  exitStatus -= 1;
  return exitStatus;
}

int main(int argc, char *argv[]) {
  int status = 0;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    if (argc < 2) {
      cerr << " usage " << argv[0] << " <input_file_name> " << endl;
      exit(0);
    }

    char *inFile = argv[1];
    ParmParse pp(argc - 2, argv + 2, NULL, inFile);

    PoissonParameters params;
    Vector<DisjointBoxLayout> grids;

    // read params from file
    getPoissonParameters(params);
    int nlevels = params.numLevels;
    Vector<LevelData<FArrayBox> *> dpsi(
        nlevels, NULL); // the correction to the conformal factor
    Vector<LevelData<FArrayBox> *> rhs(nlevels, NULL); // rhs
    Vector<LevelData<FArrayBox> *> psi(nlevels, NULL); // the conformal factor
    Vector<LevelData<FArrayBox> *> phi(nlevels, NULL); // scalar field

    set_grids(grids, params);

    status = poissonSolve(dpsi, psi, phi, rhs, grids, params);

    // clear memory
    for (int level = 0; level < dpsi.size(); level++) {
      if (dpsi[level] != NULL) {
        delete dpsi[level];
        dpsi[level] = NULL;
      }
      if (psi[level] != NULL) {
        delete psi[level];
        psi[level] = NULL;
      }
      if (phi[level] != NULL) {
        delete phi[level];
        phi[level] = NULL;
      }
      if (rhs[level] != NULL) {
        delete rhs[level];
        rhs[level] = NULL;
      }
    }

  } // End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
