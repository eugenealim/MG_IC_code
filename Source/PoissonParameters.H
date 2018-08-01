#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSONPARAMETERS_H_
#define _POISSONPARAMETERS_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "HamiltonianPoissonOperator.H"
#include "ParmParse.H"
#include "FArrayBox.H"

#include "UsingNamespace.H"

///
/**
 */

class PoissonParameters
{
public:
  IntVect      nCells;
  int          maxGridSize;
  int          blockFactor;
  int          bufferSize;
  Real         fillRatio;
  Real         refineThresh;
  int          domBcType;
  int          maxLevel;
  int          numLevels;
  Vector<int>  refRatio;
  ProblemDomain coarsestDomain;
  Real         coarsestDx;
  RealVect     domainLength;
  RealVect     probLo;
  RealVect     probHi;
  Real         alpha;
  Real         beta;
  Real         kappa_sq;
  Real         initial_chi; // initial constant chi
  Real         initial_phi; // initial constant phi
  Real         constant_K; // initial constant K
  int          rho_type; // 0 = gaussian, 1 = wavemodes
  Real         rho_scale; // scale of rho gaussian
  Real         rho_strength; // strength of rho gaussian, rho = strength*Exp(-r^2/scale)
  Vector<Real> rho_center_1;
  Vector<Real> rho_center_2;
  Vector<Real> rho_center_3;
  Vector<Real> rho_k; // wave lengths for wave modes
  Vector<Real> rho_amplitude; // three amplitude for wave modes
  Real         rho_baseline; // baseline
  int          coefficient_average_type;
  int          verbosity;
  Vector<int>  periodic; // 1 1 1 for periodic, 0 0 0 for non-periodic
};

///get stuff from input file
extern void getPoissonParameters (PoissonParameters&  a_params);

#endif
