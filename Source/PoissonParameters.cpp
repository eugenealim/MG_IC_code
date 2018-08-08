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
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

void getPoissonParameters(PoissonParameters &a_params) {
  ParmParse pp;

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells", nCellsArray, 0, SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++) {
    a_params.nCells[idir] = nCellsArray[idir];
  }

  Vector<int> is_periodic(SpaceDim, false);
  pp.queryarr("periodic", is_periodic, 0, SpaceDim);

  pp.get("refine_threshold", a_params.refineThresh);
  pp.get("block_factor", a_params.blockFactor);
  pp.get("fill_ratio", a_params.fillRatio);
  pp.get("buffer_size", a_params.bufferSize);
  pp.get("alpha", a_params.alpha);
  pp.get("beta", a_params.beta);
  pp.get("kappa_sq", a_params.kappa_sq); // this is 8piG
  pp.get("initial_psi", a_params.initial_psi);
  pp.get("initial_phi", a_params.initial_phi);
  pp.get("constant_K", a_params.constant_K);

  // Get rho params
  pp.get("rho_type", a_params.rho_type); // 0 = gaussian, 1 = waves
  // Gaussian rho = strength * Exp[-r^2/scale]
  pp.get("rho_scale", a_params.rho_scale);
  pp.get("rho_strength", a_params.rho_strength);
  pp.getarr("rho_center_1", a_params.rho_center_1, 0, SpaceDim);
  pp.getarr("rho_center_2", a_params.rho_center_2, 0, SpaceDim);
  pp.getarr("rho_center_3", a_params.rho_center_3, 0, SpaceDim);

  // waves rho = SUM_i rho_amplitude_i * Sin[2pi rho_k_i x_i]
  pp.getarr("rho_k", a_params.rho_k, 0, SpaceDim);
  pp.getarr("rho_amplitude", a_params.rho_amplitude, 0, SpaceDim);
  pp.get("rho_baseline", a_params.rho_baseline);

  pout() << "alpha, beta = " << a_params.alpha << ", " << a_params.beta << endl;

  // set to a bogus default value, so we only break from solver
  // default if it's set to something real
  a_params.coefficient_average_type = -1;
  if (pp.contains("coefficient_average_type")) {
    std::string tempString;
    pp.get("coefficient_average_type", tempString);
    if (tempString == "arithmetic") {
      a_params.coefficient_average_type = CoarseAverage::arithmetic;
    } else if (tempString == "harmonic") {
      a_params.coefficient_average_type = CoarseAverage::harmonic;
    } else {
      MayDay::Error("bad coefficient_average_type in input");
    }
  } // end if an average_type is present in inputs

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio", a_params.refRatio, 0, a_params.numLevels);
  a_params.verbosity = 3;
  pp.query("verbosity", a_params.verbosity);

  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  Box crseDomBox(lo, hi);

  // Eugene : let's read in the periodic info first, implement later
  //  bool is_periodic[SpaceDim];
  Vector<int> is_periodic_int;
  pp.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir = 0; dir < SpaceDim; dir++) {
    is_periodic[dir] = (is_periodic_int[dir] == 1);
  }
  a_params.periodic = is_periodic_int;

  ProblemDomain crseDom(crseDomBox);
  for (int dir = 0; dir < SpaceDim; dir++) {
    crseDom.setPeriodic(dir, is_periodic[dir]);
  }
  a_params.coarsestDomain = crseDom;

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length", dLArray, 0, SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++) {
    a_params.domainLength[idir] = dLArray[idir];
  }

  pp.get("max_grid_size", a_params.maxGridSize);

  // derived stuff
  a_params.coarsestDx = a_params.domainLength[0] / a_params.nCells[0];

  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;
}
