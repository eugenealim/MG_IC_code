#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include "ParmParse.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "BRMeshRefine.H"
#include "FABView.H"
#include "DebugDump.H"
#include "VClocalFuncs.H"
#include "PoissonParameters.H"
#include "MultilevelLinearOp.H"
#include "BiCGStabSolver.H"
#include "AMRIO.H"
#include "SetOperatorCoefficients.H"
#include "WriteOutput.H"
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

//Sets up and runs the solver
int poissonSolve(Vector<LevelData<FArrayBox>* >& a_chi,
                 Vector<LevelData<FArrayBox>* >& a_rhs,
                 Vector<LevelData<FArrayBox>* >& a_psi,
                 const Vector< DisjointBoxLayout >&   a_grids,
                 const PoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_chi.resize(nlevels);
  a_rhs.resize(nlevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > aCoef(nlevels);
  Vector<RefCountedPtr<LevelData<FArrayBox> > > cCoef(nlevels);
  Vector<RefCountedPtr<LevelData<FluxBox> > > bCoef(nlevels);
  Vector<ProblemDomain> vectDomains(nlevels);
  Vector<RealVect> vectDx(nlevels);

  RealVect dxLev = RealVect::Unit;
  dxLev *=a_params.coarsestDx;
  ProblemDomain domLev(a_params.coarsestDomain);

  // Declare variables here
  for (int ilev = 0; ilev < nlevels; ilev++)
    {
      a_rhs[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Zero);
      a_chi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      a_psi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      aCoef[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
      cCoef[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
      bCoef[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      vectDomains[ilev] = domLev;
      vectDx[ilev] = dxLev;

      set_initial_chi (*a_chi[ilev], vectDx[ilev], a_params); 
      set_initial_psi (*a_psi[ilev], vectDx[ilev], a_params); 

      //prepare dx, domain for next level
      dxLev /=      a_params.refRatio[ilev];
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

  BiCGStabSolver<Vector<LevelData<FArrayBox>* > > solver; //define solver object

  // Main loop

  int NL_iter = 0;
  int max_NL_iter;
  pp.query("max_NL_iterations", max_NL_iter);

  for (NL_iter ; NL_iter < max_NL_iter; NL_iter++)
  {

    pout() << "Main Loop Iteration " << (NL_iter+1) << " out of " << max_NL_iter << endl;

    // Assign values here
    for (int ilev = 0; ilev < nlevels; ilev++)
      {
        setACoef(*aCoef[ilev], *a_chi[ilev], a_params, vectDx[ilev]);
        setBCoef(*bCoef[ilev], a_params, vectDx[ilev]);
        setCCoef(*cCoef[ilev], *a_chi[ilev], a_params, vectDx[ilev]);

//        for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
//          {
//            (*a_chi[ilev])[dit()].setVal(0.);
//          }

        setRHS (*a_rhs[ilev], *a_chi[ilev], *a_psi[ilev], vectDx[ilev], a_params);

      }

    // set up solver
    RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > > opFactory
      = RefCountedPtr<AMRLevelOpFactory<LevelData<FArrayBox> > >
          (defineOperatorFactory(a_grids, vectDomains, aCoef, bCoef, cCoef, a_params));

    mlOp.define(a_grids, a_params.refRatio, vectDomains,
                vectDx, opFactory, lBase);


    bool homogeneousBC = false;
    solver.define(&mlOp, homogeneousBC);
    solver.m_verbosity = a_params.verbosity;
    solver.m_normType = 0;
    solver.m_eps = tolerance;
    solver.m_imax = max_iter;

    outputData(a_chi, a_rhs, a_grids, a_params, NL_iter);

    solver.solve(a_chi, a_rhs);

  } // end NL_iter loop


  int exitStatus = solver.m_exitStatus;
  // note that for AMRMultiGrid, success = 1.
  exitStatus -= 1;
  return exitStatus;

}


int main(int argc, char* argv[])
{
  int status = 0;
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  // Scoping trick
  {
    if (argc < 2)
      {
        cerr << " usage " << argv[0] << " <input_file_name> " << endl;
        exit(0);
      }

    char* inFile = argv[1];
    ParmParse pp(argc-2,argv+2,NULL,inFile);

    PoissonParameters param;
    Vector<DisjointBoxLayout> grids;

    //read params from file
    getPoissonParameters(param);
    int nlevels = param.numLevels;
    Vector<LevelData<FArrayBox>* > chi(nlevels, NULL); // really chi -- will fixed eventually EUGENE
    Vector<LevelData<FArrayBox>* > psi(nlevels, NULL); // scalar field 
    Vector<LevelData<FArrayBox>* > rhs(nlevels, NULL); // rhs
    
    setGrids(grids,  param);

    status = poissonSolve(chi, rhs, psi, grids,  param);

    int dofileout;
    pp.get("write_output", dofileout);
    if (dofileout == 1)
      {
        outputData(chi, rhs, grids, param, 999);
      }

    // clear memory
    for (int level = 0; level<chi.size(); level++)
      {
        if (chi[level] != NULL)
          {
            delete chi[level];
            chi[level] = NULL;
          }
        if (rhs[level] != NULL)
          {
            delete rhs[level];
            rhs[level] = NULL;
          }
      }

  }// End scoping trick

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return status;
}
