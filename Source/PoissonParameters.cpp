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
  pp.get("G_Newton", a_params.G_Newton);
  pp.get("initial_psi", a_params.initial_psi);
  pp.get("constant_K", a_params.constant_K);

  // Gaussian phi = strength * Exp[-r^2/scale]
  pp.get("phi_scale", a_params.phi_scale);
  pp.get("phi_strength", a_params.phi_strength);
  pp.getarr("phi_center", a_params.phi_center, 0, SpaceDim);

  // print out the overall coeffs just to be sure we have selected them
  // correctly
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
