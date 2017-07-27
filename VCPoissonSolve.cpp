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


/* see BoxTools/BaseFabMacros.H  for ForAllXBNN loop macro */
void
setACoef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_phi,
         const VCPoissonParameters& a_params,
         const RealVect& a_dx)
{
  RealVect pos;
  int num;

  if (a_params.ACoeftype == 0) // A(x) = x
    {
    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& aCoef = a_aCoef[dit];
        ForAllXBNN(Real, aCoef, aCoef.box(), 0, aCoef.nComp());  // Macros in BoxTools/BaseFabMacros.H
        {
          num = nR;
          D_TERM(pos[0]=a_dx[0]*(iR+0.5);,
                pos[1]=a_dx[1]*(jR+0.5);,
                pos[2]=a_dx[2]*(kR+0.5));
          aCoefR = pos[0];// this just fixed the aCoef to be x (which is what was also set in functionsF.ChF
        }EndFor;
      } // end loop over grids
    }
  else if (a_params.ACoeftype == 1) // ACoef is constant 1
    {
      DataIterator dit = a_aCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& aCoef = a_aCoef[dit];
          ForAllXBNN(Real, aCoef, aCoef.box(), 0, aCoef.nComp());
          {
            aCoefR = 1.0; // constant 
          }
        }EndFor;
    } 
  else if (a_params.ACoeftype == 2) // ACoef = -5/8 * phi_0^4 M(phi_0,K) 
    {
      DataIterator dit = a_aCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          // cell centered
          RealVect ccOffset = 0.5*a_dx*RealVect::Unit; 

          FArrayBox& aCoef = a_aCoef[dit];
          FArrayBox& phi = a_phi[dit];
          Box thisBox = aCoef.box();
         
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit)
          {
            IntVect iv = bit();
            RealVect loc(iv);
            loc *= a_dx;
            loc += ccOffset;

            Real M = M_value(loc, a_params);

            Real phi_0 = phi(iv,0);

            aCoef(iv,0) = (-0.625) * M * (phi_0 * phi_0 * phi_0 * phi_0);
          }
        }
    } // end ACoef

}


void
setBCoef(LevelData<FluxBox>& a_bCoef,
         const VCPoissonParameters& a_params,
         const RealVect& a_dx)
{
  if (a_params.BCoeftype == 0 ) // original B(x)=(x,y,z)
    {
    DataIterator dit = a_bCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = a_bCoef[dit];
        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& dirFlux = thisBCoef[dir];
            const Box& dirBox = dirFlux.box();
            // this sets up a vector which is 0 in the dir
            // direct and 0.5 in the other (cell-centered) directions
            RealVect offsets = BASISREALV(dir);
            RealVect pos;
            offsets -= RealVect::Unit;
            offsets *= -0.5;
            int n;
            ForAllXBNN(Real, dirFlux, dirBox, 0, dirFlux.nComp())
              {
                n = nR;
                D_TERM(pos[0] = a_dx[0]*(iR+offsets[0]);,
                      pos[1] = a_dx[1]*(jR+offsets[1]);,
                      pos[2] = a_dx[2]*(kR+offsets[2]));
                dirFluxR = D_TERM(pos[0], +pos[1], +pos[2]);
                // constant-coefficient
                //dirFluxR = 1.0;
              }EndFor
          } // end loop over directions
      }
    }
  else if (a_params.BCoeftype == 1) // B(x) = (1,1,1)
    {
    DataIterator dit = a_bCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FluxBox& thisBCoef = a_bCoef[dit];
        for (int dir=0; dir<SpaceDim; dir++)
          {
            FArrayBox& dirFlux = thisBCoef[dir];
            const Box& dirBox = dirFlux.box();
            // this sets up a vector which is 0 in the dir
            // direct and 0.5 in the other (cell-centered) directions
            ForAllXBNN(Real, dirFlux, dirBox, 0, dirFlux.nComp())
              {
                dirFluxR = 1.0; // constant direction
              }EndFor
          } // end loop over directions
      }
    } // end BCoeftype
}

void
setCCoef(LevelData<FArrayBox>& a_cCoef,
         LevelData<FArrayBox>& a_phi,
         const VCPoissonParameters& a_params,
         const RealVect& a_dx)
{
  RealVect pos;
  int num;

  if (a_params.CCoeftype == 0) // C(x) = x
    {
    DataIterator dit = a_cCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& cCoef = a_cCoef[dit];
        ForAllXBNN(Real, cCoef, cCoef.box(), 0, cCoef.nComp());  // Macros in BoxTools/BaseFabMacros.H
        {
          num = nR;
          D_TERM(pos[0]=a_dx[0]*(iR+0.5);,
                pos[1]=a_dx[1]*(jR+0.5);,
                pos[2]=a_dx[2]*(kR+0.5));
          cCoefR = pos[0];// this just fixed the aCoef to be x (which is what was also set in functionsF.Chf:w
          // Eugene change 
  //        aCoefR = D_TERM(pos[0]*,pos[1]*,pos[2]);
          // constant-coefficient
          //aCoefR = 1.0;
        }EndFor;
      } // end loop over grids
    }
  else if (a_params.CCoeftype == 1) // CCoef is constant 1
    {
      DataIterator dit = a_cCoef.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& cCoef = a_cCoef[dit];
          ForAllXBNN(Real, cCoef, cCoef.box(), 0, cCoef.nComp());
          {
            cCoefR = 1.0; // constant 
          }
        }EndFor;
    } // end ACoeftype
}

/***************/
void outputData(const Vector<LevelData<FArrayBox>* >&   a_phi,
                const Vector<LevelData<FArrayBox>* >&   a_rhs,
                const Vector< DisjointBoxLayout >&       a_grids,
                const VCPoissonParameters&                 a_params,
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

    int nPhiComp = a_phi[0]->nComp();
    int nRhsComp = a_rhs[0]->nComp();
    int totalComp = nPhiComp + nRhsComp;

    Vector<string> phiNames(totalComp);
    // hardwire to single-component
    CH_assert(totalComp == 2);
    phiNames[0] = "phi";
    phiNames[1] = "rhs";


    CH_assert(a_phi.size() == a_rhs.size());
    Vector<LevelData<FArrayBox>* > tempData(a_phi.size(), NULL);
    for (int level=0; level<a_phi.size(); level++)
      {
        tempData[level] = new LevelData<FArrayBox>(a_grids[level], totalComp);
        Interval phiComps(0, nPhiComp-1);
        Interval rhsComps(nPhiComp, totalComp-1);
        a_phi[level]->copyTo(a_phi[level]->interval(),
                             *tempData[level], phiComps);
        a_rhs[level]->copyTo(a_rhs[level]->interval(),
                             *tempData[level], rhsComps);
      }


    Real fakeTime = 1.0;
    Real fakeDt = 1.0;
    WriteAMRHierarchyHDF5(fileName, a_grids,
                          tempData, phiNames,
                          a_params.coarsestDomain.domainBox(),
                          a_params.coarsestDx,
                          fakeDt, fakeTime,
                          a_params.refRatio,
                          a_params.numLevels);

    // clean up temporary storage
    for (int level=0; level<a_phi.size(); level++)
      {
        delete tempData[level];
        tempData[level] = NULL;
      }
#endif

}

/******/
int poissonSolve(Vector<LevelData<FArrayBox>* >& a_phi,
                 Vector<LevelData<FArrayBox>* >& a_rhs,
                 Vector<LevelData<FArrayBox>* >& a_psi,
                 const Vector< DisjointBoxLayout >&   a_grids,
                 const VCPoissonParameters&             a_params)
{
  ParmParse pp;

  int nlevels = a_params.numLevels;
  a_phi.resize(nlevels);
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
      a_phi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      a_psi[ilev] = new LevelData<FArrayBox>(a_grids[ilev], 1,IntVect::Unit);
      aCoef[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
      cCoef[ilev] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(a_grids[ilev], 1, IntVect::Zero));
      bCoef[ilev] = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(a_grids[ilev], 1, IntVect::Zero));
      vectDomains[ilev] = domLev;
      vectDx[ilev] = dxLev;

      set_initial_phi (*a_phi[ilev], vectDx[ilev], a_params); 
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
        setACoef(*aCoef[ilev], *a_phi[ilev], a_params, vectDx[ilev]);
        setBCoef(*bCoef[ilev], a_params, vectDx[ilev]);
        setCCoef(*cCoef[ilev], *a_phi[ilev], a_params, vectDx[ilev]);

//        for (DataIterator dit = a_grids[ilev].dataIterator(); dit.ok(); ++dit)
//          {
//            (*a_phi[ilev])[dit()].setVal(0.);
//          }

        setRHS (*a_rhs[ilev], *a_phi[ilev], *a_psi[ilev], vectDx[ilev], a_params);

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

    solver.solve(a_phi, a_rhs);

    outputData(a_phi, a_rhs, a_grids, a_params, NL_iter);

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

    VCPoissonParameters param;
    Vector<DisjointBoxLayout> grids;

    //read params from file
    getPoissonParameters(param);
    int nlevels = param.numLevels;
    Vector<LevelData<FArrayBox>* > phi(nlevels, NULL); // really chi -- will fixed eventually EUGENE
    Vector<LevelData<FArrayBox>* > psi(nlevels, NULL); // scalar field 
    Vector<LevelData<FArrayBox>* > rhs(nlevels, NULL); // rhs
    
    setGrids(grids,  param);

    status = poissonSolve(phi, rhs, psi, grids,  param);

    int dofileout;
    pp.get("write_output", dofileout);
    if (dofileout == 1)
      {
        outputData(phi, rhs, grids, param, 999);
      }

    // clear memory
    for (int level = 0; level<phi.size(); level++)
      {
        if (phi[level] != NULL)
          {
            delete phi[level];
            phi[level] = NULL;
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
