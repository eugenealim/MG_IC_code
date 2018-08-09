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
#include "CH_HDF5.H"
#include "SetLevelData.H"
#include "UsingNamespace.H"
#include "UserVariables.hpp"

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
void output_solver_data(const Vector<LevelData<FArrayBox>* >&   a_dpsi,
                      const Vector<LevelData<FArrayBox>* >&   a_psi,
                      const Vector<LevelData<FArrayBox>* >&   a_phi,
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
    int npsiComp = a_psi[0]->nComp();
    int nphiComp = a_phi[0]->nComp();
    int totalComp = ndpsiComp + nRhsComp + npsiComp + nphiComp;

    Vector<string> variable_names(totalComp);
    // hardwire to single-component
    CH_assert(totalComp == 4);
    variable_names[0] = "dpsi";
    variable_names[1] = "psi";
    variable_names[2] = "phi";
    variable_names[3] = "rhs";


    CH_assert(a_dpsi.size() == a_rhs.size());
    CH_assert(a_dpsi.size() == a_phi.size());
    CH_assert(a_dpsi.size() == a_psi.size());

    Vector<LevelData<FArrayBox>* > tempData(a_dpsi.size(), NULL);
    for (int level=0; level<a_dpsi.size(); level++)
      {
        tempData[level] = new LevelData<FArrayBox>(a_grids[level], totalComp);
        Interval dpsiComps(0, 0);
        Interval psiComps(1, 1);
        Interval phiComps(2, 2);
        Interval rhsComps(3, 3);
        a_dpsi[level]->copyTo(a_dpsi[level]->interval(),
                             *tempData[level], dpsiComps);
        a_psi[level]->copyTo(a_psi[level]->interval(),
                             *tempData[level], psiComps);
        a_phi[level]->copyTo(a_phi[level]->interval(),
                             *tempData[level], phiComps);
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

void output_final_data( const Vector<LevelData<FArrayBox>* >&   a_psi,
                      const Vector<LevelData<FArrayBox>* >&   a_phi,
                      const Vector< DisjointBoxLayout >&      a_grids,
                      const Vector<RealVect>&               a_vectDx,
                      const Vector<ProblemDomain>&          a_vectDomains,
                      const PoissonParameters&              a_params)
{
#ifdef CH_USE_HDF5

#if CH_SPACEDIM==2
    string filename("vcPoissonFinal.2d.hdf5");
#else
    string filename("vcPoissonFinal.3d.hdf5");
#endif

    // set up full data and write to hdf5 file ready for GRChombo
    HDF5Handle handle(filename, HDF5Handle::CREATE);
    Vector<LevelData<FArrayBox>* > vars_data(a_psi.size(), NULL);
    HDF5HeaderData header;

    header.m_int["max_level"] = a_params.maxLevel;
    header.m_int["num_levels"] = a_params.maxLevel + 1;
    header.m_int["iteration"] = 0;
    header.m_real["time"] = 0.0;

    for (int level = 0; level < a_psi.size(); ++level)
    {
      char headername1[100];
      sprintf(headername1, "regrid_interval_%d", level);
      header.m_int[headername1] = 1;

      char headername2[100];
      sprintf(headername2, "steps_since_regrid_%d", level);
      header.m_int[headername2] = 0;
    }

    header.m_int["num_components"] = NUM_VARS;

    char comp_str[30];
    for (int comp = 0; comp < NUM_VARS; ++comp)
    {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = UserVariables::variable_names[comp];
    }

    // write header data to file
    header.writeToFile(handle);

    // loop over levels
    for (int level=0; level<a_psi.size(); level++)
    {
        HDF5HeaderData level_header;
        // want output data to have 3 ghost cells to match GRChombo
        vars_data[level] = new LevelData<FArrayBox>(a_grids[level], NUM_VARS, 3*IntVect::Unit);

        //Add in the rest of the data in the correct form
        set_output_data(*vars_data[level], *a_psi[level], *a_phi[level], a_params);

        //Now continue the write out
        char level_str[20];
        sprintf(level_str, "%d", level);
        const std::string label = std::string("level_") + level_str;
        handle.setGroup(label);

        // assign arbitrary values to the params
        // most can be overwritten in GRChombo setup
        level_header.m_int["ref_ratio"] = a_params.refRatio[level];
        level_header.m_int["tag_buffer_size"] = 3;
        level_header.m_real["dx"] = a_vectDx[level][0];
        level_header.m_real["dt"] = 0.25*a_vectDx[level][0];
        level_header.m_real["time"] = 0.0;
        level_header.m_box["prob_domain"] = a_vectDomains[level].domainBox();

        // Setup the periodicity info - for GRChombo we always treat it as periodic
        for (int dir = 0; dir < SpaceDim; ++dir)
        {
            char dir_str[20];
            sprintf(dir_str, "%d", dir);
            const std::string periodic_label =
                std::string("is_periodic_") + dir_str;
            level_header.m_int[periodic_label] = true;
        }

        level_header.writeToFile(handle);
        write(handle, a_grids[level]);
        write(handle, *vars_data[level], "data");
    }

    //shut the file
    handle.close();

    // clean up temporary storage
    for (int level=0; level<a_psi.size(); level++)
    {
        delete vars_data[level];
        vars_data[level] = NULL;
    }

#endif

}


#endif