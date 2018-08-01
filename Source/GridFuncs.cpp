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
#include "GridFuncs.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "HamiltonianPoissonOperatorFactory.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

int setGrids(Vector<DisjointBoxLayout> &vectGrids,
             PoissonParameters &a_params) {
  Vector<ProblemDomain> vectDomain;
  Vector<Real> vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);

  int numlevels = a_params.numLevels;

  ParmParse pp;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
  if (pp.contains("read_in_grids")) {
    pp.get("read_in_grids", readInGrids);
  }

  if (readInGrids) {

    ProblemDomain levDomain = a_params.coarsestDomain;
    for (int ilev = 0; ilev < a_params.numLevels; ilev++) {
      Vector<Box> boxes;
      char boxCountVar[100];
      int boxCount;
      sprintf(boxCountVar, "level_%d_box_count", ilev);
      pp.get(boxCountVar, boxCount);
      boxes.resize(boxCount);
      for (int ibox = 0; ibox < boxCount; ibox++) {
        char boxLoVar[100];
        char boxHiVar[100];
        sprintf(boxLoVar, "level_%d_box_%d_lo", ilev, ibox);
        sprintf(boxHiVar, "level_%d_box_%d_hi", ilev, ibox);
        Vector<int> boxLo, boxHi;
        pp.getarr(boxLoVar, boxLo, 0, SpaceDim);
        pp.getarr(boxHiVar, boxHi, 0, SpaceDim);
        IntVect ivLo(D_DECL(boxLo[0], boxLo[1], boxLo[2]));
        IntVect ivHi(D_DECL(boxHi[0], boxHi[1], boxHi[2]));
        boxes[ibox] = Box(ivLo, ivHi);
        if (!levDomain.contains(boxes[ibox])) {
          MayDay::Error("box outside of domain");
        }
      }
      // check to see if level 0 domain is covered
      if (ilev == 0) {
        IntVectSet ivDom(levDomain.domainBox());
        for (int ibox = 0; ibox < boxes.size(); ibox++) {
          ivDom -= boxes[ibox];
        }
        if (!ivDom.isEmpty()) {
          MayDay::Error("level 0 boxes must cover the domain");
        }
      }
      Vector<int> proc(a_params.numLevels);
      LoadBalance(proc, boxes);
      vectGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
      levDomain.refine(a_params.refRatio[ilev]);
    }

  } else // no read_grids
  {
    pout() << "tagging on gradient of RHS" << endl;
    int maxLevel = numlevels - 1;
    Vector<Vector<Box>> newBoxes(numlevels);
    Vector<Vector<Box>> oldBoxes(numlevels);

    // determine grids dynamically, based on grad(RHS)
    // will need temp storage for RHS
    Vector<LevelData<FArrayBox> *> vectRHS(maxLevel + 1, NULL);
    int ncomps = 1;

    // define base level first
    Vector<Vector<int>> procAssign(maxLevel + 1);
    domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize,
                a_params.blockFactor);
    procAssign[0].resize(oldBoxes[0].size());
    LoadBalance(procAssign[0], oldBoxes[0]);

    vectGrids[0].define(oldBoxes[0], procAssign[0], vectDomain[0]);

    vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomps, IntVect::Zero);

    int topLevel = 0;

    bool moreLevels = (maxLevel > 0);

    int nesting_radius = 2;
    // create grid generation object
    BRMeshRefine meshrefine(vectDomain[0], a_params.refRatio,
                            a_params.fillRatio, a_params.blockFactor,
                            nesting_radius, a_params.maxGridSize);

    while (moreLevels) {
      // default is moreLevels = false
      // (only repeat loop in the case where a new level
      // is generated which is still less than maxLevel)
      moreLevels = false;

      int baseLevel = 0;
      int oldTopLevel = topLevel;

      // now initialize RHS for this existing hierarchy

      for (int level = 0; level <= topLevel; level++) {
        RealVect dxLevel = vectDx[level] * RealVect::Unit;

        LevelData<FArrayBox> *DummyChi; // REMEMBER TO FIX THIS
        LevelData<FArrayBox> *DummyPhi; // REMEMBER TO FIX THIS

        DummyChi =
            new LevelData<FArrayBox>(vectGrids[level], ncomps, IntVect::Unit);
        DummyPhi =
            new LevelData<FArrayBox>(vectGrids[level], ncomps, IntVect::Unit);

        set_initial_chi(*DummyChi, dxLevel, a_params);
        set_initial_phi(*DummyPhi, dxLevel, a_params);
        setRHS(*vectRHS[level], *DummyChi, *DummyPhi, dxLevel, a_params);
      }

      Vector<IntVectSet> tagVect(topLevel + 1);
      int tags_grow = 1;
      tagCells(vectRHS, tagVect, vectDx, vectDomain, a_params.refineThresh,
               tags_grow, baseLevel, topLevel + 1);

      int new_finest =
          meshrefine.regrid(newBoxes, tagVect, baseLevel, topLevel, oldBoxes);

      if (new_finest > topLevel) {
        topLevel++;
      }

      oldBoxes = newBoxes;

      //  no need to do this for the base level (already done)
      for (int lev = 1; lev <= topLevel; lev++) {
        // do load balancing
        procAssign[lev].resize(newBoxes[lev].size());
        LoadBalance(procAssign[lev], newBoxes[lev]);
        const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                       vectDomain[lev]);
        vectGrids[lev] = newDBL;
        delete vectRHS[lev];
        vectRHS[lev] =
            new LevelData<FArrayBox>(vectGrids[lev], ncomps, IntVect::Zero);
      } // end loop over levels for initialization

      // figure out whether we need another pass through grid generation
      if ((topLevel < maxLevel) && (topLevel > oldTopLevel))
        moreLevels = true;

    } // end while moreLevels loop
    // clean up temp storage
    for (int ilev = 0; ilev < vectRHS.size(); ilev++) {
      if (vectRHS[ilev] != NULL) {
        delete vectRHS[ilev];
        vectRHS[ilev] = NULL;
      }
    }
  }
  return 0;
}

// EUGENE TODO : make it more general
void set_initial_chi(LevelData<FArrayBox> &a_chi,

                     const RealVect &a_dx, const PoissonParameters &a_params) {

  CH_assert(a_chi.nComp() == 1);

  int comp = 0;
  for (DataIterator dit = a_chi.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &thisCHIbox = a_chi[dit()];
    Box b = thisCHIbox.box();
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      for (int comp = 0; comp < a_chi.nComp(); ++comp) {
        thisCHIbox(iv, comp) = a_params.initial_chi;
      }
    }
  }
} // end set_initial_chi

// EUGENE : temporary
void set_initial_phi(LevelData<FArrayBox> &a_phi, const RealVect &a_dx,
                     const PoissonParameters &a_params) {

  CH_assert(a_phi.nComp() == 1);
  int comp = 0;
  for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &thisPHIbox = a_phi[dit()];
    Box b = thisPHIbox.box();
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      for (int comp = 0; comp < a_phi.nComp(); ++comp) {
        thisPHIbox(iv, comp) = a_params.initial_phi;
      }
    }
  }
} // end set_initial_phi

/********/
void setRHS(LevelData<FArrayBox> &a_rhs, LevelData<FArrayBox> &a_chi,
            LevelData<FArrayBox> &a_phi, const RealVect &a_dx,
            const PoissonParameters &a_params) {
  CH_assert(a_rhs.nComp() == 1);

  int comp = 0;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit) {
    FArrayBox &thisRHS = a_rhs[dit()];
    thisRHS.setVal(0,0);
    Box thisBox = thisRHS.box();

    // first set up the scalar field config
    FArrayBox &thisPhi = a_phi[dit()];

    // rhs is cell-centered...
    RealVect ccOffset = 0.5 * a_dx * RealVect::Unit;

    // put in a gaussian for now
    int numGaussians = 1;
    Vector<RealVect> center(numGaussians, RealVect::Zero);
    Vector<Real> scale(numGaussians, 1.0);
    Vector<Real> strength(numGaussians, 1.0);
    BoxIterator bit(thisRHS.box());

    // one central gaussian
    strength[0] = 1.0;
    scale[0] = 1.0e-2; // variance
    center[0] = 0.5 * RealVect::Unit;

    thisPhi.setVal(0, 0);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect iv = bit();
      RealVect loc(iv);
      loc *= a_dx;
      loc += ccOffset;

      for (int n = 0; n < numGaussians; n++) {
        RealVect dist = loc - center[n];
        Real radSqr =
            D_TERM(dist[0] * dist[0], +dist[1] * dist[1], +dist[2] * dist[2]);

        Real val = strength[n] * exp(-radSqr / scale[n]);
        thisPhi(iv, 0) += val;
        thisRHS(iv, 0) += val;
      }
    }
  }
} // end set_RHS

/*
  Set grid hierarchy from input file
 */
void getDomainsAndDxes(Vector<ProblemDomain> &vectDomain, Vector<Real> &vectDx,
                       PoissonParameters &a_params) {

  vectDomain.resize(a_params.numLevels);
  vectDx.resize(a_params.numLevels);
  vectDx[0] = a_params.coarsestDx;
  for (int ilev = 1; ilev < a_params.numLevels; ilev++) {
    vectDx[ilev] = vectDx[ilev - 1] / a_params.refRatio[ilev - 1];
  }

  vectDomain[0] = a_params.coarsestDomain;
  for (int ilev = 1; ilev < a_params.numLevels; ilev++) {
    vectDomain[ilev] =
        refine(vectDomain[ilev - 1], a_params.refRatio[ilev - 1]);
  }
}

void
setACoef(LevelData<FArrayBox>& a_aCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
    RealVect pos;
    int num;

    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        // cell centered
        RealVect ccOffset = 0.5*a_dx*RealVect::Unit; 

        FArrayBox& aCoef = a_aCoef[dit];
        FArrayBox& chi = a_chi[dit];
        Box thisBox = aCoef.box();
     
        BoxIterator bit(thisBox);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc(iv);
            loc *= a_dx;
            loc += ccOffset;

            Real M = calculate_M_value(loc, a_params);

            Real chi_0 = chi(iv,0);

            aCoef(iv,0) = (-0.625) * M * (chi_0 * chi_0 * chi_0 * chi_0);
        }
    }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
void
setBCoef(LevelData<FArrayBox>& a_bCoef,
         LevelData<FArrayBox>& a_chi,
         const PoissonParameters& a_params,
         const RealVect& a_dx)
{
  RealVect pos;
  int num;

  DataIterator dit = a_bCoef.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& bCoef = a_bCoef[dit];
      ForAllXBNN(Real, bCoef, bCoef.box(), 0, bCoef.nComp());
      {
        bCoefR = 1.0; // constant 
      }
    }EndFor;
}

// M(K, rho) = 2/3K^2 - 16piG rho
// for now we put rho to be a gaussian with center (0.5,0.5,0.5)
Real calculate_M_value(RealVect &loc, const PoissonParameters &a_params) {
  Real rho;

  RealVect center;
  center = RealVect(
      D_DECL(a_params.rho_center_1[0], a_params.rho_center_1[1],
             a_params.rho_center_1[2])); // put it in the middle for now

  RealVect dist = loc - center;
  Real radSqr =
      D_TERM(dist[0] * dist[0], +dist[1] * dist[1], +dist[2] * dist[2]);

  rho = a_params.rho_strength * exp(-radSqr / a_params.rho_scale);

  return ((2.0 / 3.0) * (a_params.constant_K * a_params.constant_K) -
          2.0 * (a_params.kappa_sq) * rho);
}

/*
  tag cells for refinement based on magnitude(RHS)
*/
void tagCells(Vector<LevelData<FArrayBox> *> &vectRHS,
              Vector<IntVectSet> &tagVect, Vector<Real> &vectDx,
              Vector<ProblemDomain> &vectDomain, const Real refine_thresh,
              const int tags_grow, const int baseLevel, int numLevels) {
  for (int lev = baseLevel; lev != numLevels; lev++) {
    IntVectSet local_tags;
    LevelData<FArrayBox> &levelRhs = *vectRHS[lev];
    DisjointBoxLayout level_domain = levelRhs.getBoxes();
    DataIterator dit = levelRhs.dataIterator();

    Real maxRHS = 0;

    maxRHS = norm(levelRhs, levelRhs.interval(), 0);

    Real tagVal = maxRHS * refine_thresh;

    // now loop through grids and tag cells where RHS > tagVal
    for (dit.reset(); dit.ok(); ++dit) {
      const Box thisBox = level_domain.get(dit());
      const FArrayBox &thisRhs = levelRhs[dit()];
      BoxIterator bit(thisBox);
      for (bit.begin(); bit.ok(); ++bit) {
        const IntVect &iv = bit();
        if (abs(thisRhs(iv)) >= tagVal)
          local_tags |= iv;
      }
    } // end loop over grids on this level

    local_tags.grow(tags_grow);
    const Box &domainBox = vectDomain[lev].domainBox();
    local_tags &= domainBox;

    tagVect[lev] = local_tags;

  } // end loop over levels
}
