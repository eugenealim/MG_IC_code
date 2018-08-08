#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SETGRIDS_H_
#define _SETGRIDS_H_

#include <iostream>
using std::cerr;

#include "DisjointBoxLayout.H"
#include "AMRMultiGrid.H"
#include "HamiltonianPoissonOperator.H"
#include "HamiltonianPoissonOperatorFactory.H"
#include "ParmParse.H"
#include "SetLevelData.H"
#include "PoissonParameters.H"
#include "FArrayBox.H"

#include "UsingNamespace.H"

// KC TODO - make this a namespace?
// various functions operating to set up grids
extern int set_grids(Vector<DisjointBoxLayout>& vectGrids,
                    PoissonParameters&         a_params);

extern void set_domains_and_dx(Vector<ProblemDomain>&     vectDomain,
                               Vector<Real>&              vectDx,
                               PoissonParameters&         a_params);

extern void
set_tag_cells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels) ;

#endif
