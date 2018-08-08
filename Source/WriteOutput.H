#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _WRITEOUTPUT_H_
#define _WRITEOUTPUT_H_

#include <iostream>
#include "ParmParse.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "BRMeshRefine.H"
#include "FABView.H"
#include "DebugDump.H"
#include "PoissonParameters.H"
#include "MultilevelLinearOp.H"
#include "BiCGStabSolver.H"
#include "AMRIO.H"
#include "UsingNamespace.H"

#ifdef CH_Linux
// Should be undefined by default
//#define TRAP_FPE
#undef  TRAP_FPE
#endif

#ifdef TRAP_FPE
static void enableFpExceptions();
#endif


using std::cerr;

/***************/
void outputData(const Vector<LevelData<FArrayBox>* >&   a_dpsi,
                const Vector<LevelData<FArrayBox>* >&   a_rhs,
                const Vector< DisjointBoxLayout >&      a_grids,
                const PoissonParameters&              a_params,
                const int iter)
{
#ifdef CH_USE_HDF5

#if CH_SPACEDIM==2
    string fileName("vcPoissonOut.2d");
#else
    string fileName("vcPoissonOut.3d");
#endif

    char suffix[30];
    sprintf(suffix, "_%d.hdf5",iter);
    fileName += suffix;

    int ndpsiComp = a_dpsi[0]->nComp();
    int nRhsComp = a_rhs[0]->nComp();
    int totalComp = ndpsiComp + nRhsComp;

    Vector<string> variable_names(totalComp);
    // hardwire to single-component
    CH_assert(totalComp == 2);
    variable_names[0] = "dpsi";
    variable_names[1] = "rhs";


    CH_assert(a_dpsi.size() == a_rhs.size());
    Vector<LevelData<FArrayBox>* > tempData(a_dpsi.size(), NULL);
    for (int level=0; level<a_dpsi.size(); level++)
      {
        tempData[level] = new LevelData<FArrayBox>(a_grids[level], totalComp);
        Interval dpsiComps(0, ndpsiComp-1);
        Interval rhsComps(ndpsiComp, totalComp-1);
        a_dpsi[level]->copyTo(a_dpsi[level]->interval(),
                             *tempData[level], dpsiComps);
        a_rhs[level]->copyTo(a_rhs[level]->interval(),
                             *tempData[level], rhsComps);
      }


    Real fakeTime = 0.0;
    Real fakeDt = 1.0;
    WriteAMRHierarchyHDF5(fileName, a_grids,
                          tempData, variable_names,
                          a_params.coarsestDomain.domainBox(),
                          a_params.coarsestDx,
                          fakeDt, fakeTime,
                          a_params.refRatio,
                          a_params.numLevels);

    // clean up temporary storage
    for (int level=0; level<a_dpsi.size(); level++)
      {
        delete tempData[level];
        tempData[level] = NULL;
      }
#endif

}

#endif