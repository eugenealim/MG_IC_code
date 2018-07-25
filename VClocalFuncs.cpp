#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include "parstream.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "VClocalFuncs.H"
#include "functionsF_F.H"
#include "AMRIO.H"
#include "BoxIterator.H"
#include "BiCGStabSolver.H"
#include "computeNorm.H"
#include "BCFunc.H"
#include "CoarseAverage.H"
#include "CONSTANTS.H"

#include "UsingNamespace.H"

std::vector<bool> GlobalBCRS::s_printedThatLo = std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi = std::vector<bool>(SpaceDim, false);
std::vector<int>  GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int>  GlobalBCRS::s_bcHi = std::vector<int>();
RealVect          GlobalBCRS::s_trigvec = RealVect::Zero;
bool              GlobalBCRS::s_areBCsParsed= false;
bool              GlobalBCRS::s_valueParsed= false;
bool              GlobalBCRS::s_trigParsed= false;

// BCValueHolder class, which is a pointer to a void-type function with the 4 arguements
// given pos [x,y,z] position on center of cell edge int dir direction, x being 0 int side -1 for low, +1 = high, fill in the a_values array

void ParseValue(Real* pos,
                int* dir,
                Side::LoHiSide* side,
                Real* a_values)
{
  ParmParse pp;
  Real bcVal;
  pp.get("bc_value",bcVal);
  a_values[0]=bcVal;
}


void ParseBC(FArrayBox& a_state,
             const Box& a_valid,
             const ProblemDomain& a_domain,
             Real a_dx,
             bool a_homogeneous)
{
  if (!a_domain.domainBox().contains(a_state.box()))
    {

      if (!GlobalBCRS::s_areBCsParsed)
        {
          ParmParse pp;
          pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
          pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
          GlobalBCRS::s_areBCsParsed = true;
        }

      Box valid = a_valid;

      for (int i=0; i<CH_SPACEDIM; ++i)
        {

      // periodic?
          if ( !a_domain.isPeriodic(i))
          {
      // periodic?

            Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
            Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
            if (!a_domain.domainBox().contains(ghostBoxLo))
              {
                if (GlobalBCRS::s_bcLo[i] == 1)
                  {
                    if (!GlobalBCRS::s_printedThatLo[i])
                      {
                        GlobalBCRS::s_printedThatLo[i] = true;
                        pout() << "const neum bcs lo for direction " << i << endl;
                      }
                    NeumBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          ParseValue, //BCValueHolder class
                          i,
                          Side::Lo);
                  }
                else if (GlobalBCRS::s_bcLo[i] == 2)
                  {
                    if (!GlobalBCRS::s_printedThatLo[i])
                      {
                        GlobalBCRS::s_printedThatLo[i] = true;
                        pout() << "trig neum bcs lo for direction " << i << endl;
                      }
                    NeumBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          TrigValueNeum,
                          i,
                          Side::Lo);
                  }
                else if (GlobalBCRS::s_bcLo[i] == 0)
                  {
                    if (!GlobalBCRS::s_printedThatLo[i])
                      {
                        GlobalBCRS::s_printedThatLo[i] = true;
                        pout() << "const diri bcs lo for direction " << i << endl;
                      }
                    DiriBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          ParseValue,
                          i,
                          Side::Lo);

                  }
                else if (GlobalBCRS::s_bcLo[i] == 3)
                  {
                    if (!GlobalBCRS::s_printedThatLo[i])
                      {
                        GlobalBCRS::s_printedThatLo[i] = true;
                        pout() << "trig diri bcs lo for direction " << i << endl;
                      }
                    DiriBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          TrigValueDiri,
                          i,
                          Side::Lo);

                  }
                else if (GlobalBCRS::s_bcLo[i] == 4)
                  {
                    if (!GlobalBCRS::s_printedThatLo[i])
                      {
                        GlobalBCRS::s_printedThatLo[i] = true;
                        pout() << "periodic bcs lo for direction " << i << endl;
                      }

                  }
                else
                  {
                    MayDay::Error("bogus bc flag lo");
                  }
              }

            if (!a_domain.domainBox().contains(ghostBoxHi))
              {
                if (GlobalBCRS::s_bcHi[i] == 1)
                  {
                    if (!GlobalBCRS::s_printedThatHi[i])
                      {
                        GlobalBCRS::s_printedThatHi[i] = true;
                        pout() << "const neum bcs hi for direction " << i << endl;
                      }
                    NeumBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          ParseValue,
                          i,
                          Side::Hi);
                  }
                else if (GlobalBCRS::s_bcHi[i] == 2)
                  {
                    if (!GlobalBCRS::s_printedThatHi[i])
                      {
                        GlobalBCRS::s_printedThatHi[i] = true;
                        pout() << "trig neum bcs hi for direction " << i << endl;
                      }
                    NeumBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          TrigValueNeum,
                          i,
                          Side::Hi);
                  }
                else if (GlobalBCRS::s_bcHi[i] == 0)
                  {
                    if (!GlobalBCRS::s_printedThatHi[i])
                      {
                        GlobalBCRS::s_printedThatHi[i] = true;
                        pout() << "const diri bcs hi for direction " << i << endl;
                      }
                    DiriBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          ParseValue,
                          i,
                          Side::Hi);
                  }
                else if (GlobalBCRS::s_bcHi[i] == 3)
                  {
                    if (!GlobalBCRS::s_printedThatHi[i])
                      {
                        GlobalBCRS::s_printedThatHi[i] = true;
                        pout() << "trig diri bcs hi for direction " << i << endl;
                      }
                    DiriBC(a_state,
                          valid,
                          a_dx,
                          a_homogeneous,
                          TrigValueDiri,
                          i,
                          Side::Hi);
                  }
                else if (GlobalBCRS::s_bcHi[i] == 4)
                  {
                    if (!GlobalBCRS::s_printedThatHi[i])
                      {
                        GlobalBCRS::s_printedThatHi[i] = true;
                        pout() << "periodic bcs hi for direction " << i << endl;
                      }
                  }
                else
                  {
                    MayDay::Error("bogus bc flag hi");
                  }
              }
          // is periodic
           }
          // is periodic
        }
    }
}


/***************/
RealVect& getTrigRV()
{
  Real pi = 4.*atan(1.0);
  if (!GlobalBCRS::s_trigParsed)
    {
      GlobalBCRS::s_trigParsed = true;
      ParmParse pp;
      std::vector<Real> trigvec(SpaceDim);
      pp.getarr("trig",trigvec,0,SpaceDim);

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          GlobalBCRS::s_trigvec[idir] = pi*trigvec[idir];
        }
    }
  return GlobalBCRS::s_trigvec;
}

/********/
void getPoissonParameters(VCPoissonParameters&  a_params)
{
  ParmParse pp;

  std::vector<int> nCellsArray(SpaceDim);
  pp.getarr("n_cells",nCellsArray,0,SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.nCells[idir] = nCellsArray[idir];
    }

  Vector<int> is_periodic(SpaceDim, false);
  pp.queryarr("periodic", is_periodic, 0, SpaceDim);

  pp.get("refine_threshold",a_params.refineThresh);
  pp.get("block_factor",a_params.blockFactor);
  pp.get("fill_ratio",a_params.fillRatio);
  pp.get("buffer_size",a_params.bufferSize);
  pp.get("alpha",a_params.alpha);
  pp.get("beta", a_params.beta);
  pp.get("gamma", a_params.gamma);
  pp.get("kappa_sq", a_params.kappa_sq); // this is 8piG
  pp.get("initial_chi", a_params.initial_chi);
  pp.get("initial_psi", a_params.initial_psi);
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

  pout() << "alpha beta gamma = " << a_params.alpha << a_params.beta << a_params.gamma << endl;

  pp.query("probtype", a_params.probtype);
  pp.query("ACoeftype", a_params.ACoeftype);
  pp.query("BCoeftype", a_params.BCoeftype);
  pp.query("CCoeftype", a_params.CCoeftype);

  // set to a bogus default value, so we only break from solver
  // default if it's set to something real
  a_params.coefficient_average_type = -1;
  if (pp.contains("coefficient_average_type"))
    {
      std::string tempString;
      pp.get("coefficient_average_type", tempString);
      if (tempString == "arithmetic")
        {
          a_params.coefficient_average_type = CoarseAverage::arithmetic;
        }
      else if (tempString == "harmonic")
        {
          a_params.coefficient_average_type = CoarseAverage::harmonic;
        }
      else
        {
          MayDay::Error("bad coefficient_average_type in input");
        }
    } // end if an average_type is present in inputs

  pp.get("max_level", a_params.maxLevel);
  a_params.numLevels = a_params.maxLevel + 1;
  pp.getarr("ref_ratio",a_params.refRatio,0,a_params.numLevels);
  a_params.verbosity = 3;
  pp.query("verbosity", a_params.verbosity);

  IntVect lo = IntVect::Zero;
  IntVect hi = a_params.nCells;
  hi -= IntVect::Unit;

  Box crseDomBox(lo,hi);

  // Eugene : let's read in the periodic info first, implement later
  //  bool is_periodic[SpaceDim];
  Vector<int> is_periodic_int;
  pp.getarr("is_periodic", is_periodic_int, 0, SpaceDim);
  for (int dir=0; dir<SpaceDim; dir++)
  {
    is_periodic[dir] = (is_periodic_int[dir] == 1);
  }
  a_params.periodic = is_periodic_int;

  ProblemDomain crseDom(crseDomBox);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      crseDom.setPeriodic(dir, is_periodic[dir]);
    }
  a_params.coarsestDomain = crseDom;

  std::vector<Real> dLArray(SpaceDim);
  pp.getarr("domain_length",dLArray,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      a_params.domainLength[idir] = dLArray[idir];
    }

  pp.get("max_grid_size",a_params.maxGridSize);

  //derived stuff
  a_params.coarsestDx = a_params.domainLength[0]/a_params.nCells[0];

  a_params.probLo = RealVect::Zero;
  a_params.probHi = RealVect::Zero;
  a_params.probHi += a_params.domainLength;


}


int setGrids(Vector<DisjointBoxLayout>& vectGrids,
             VCPoissonParameters&         a_params)
{
  Vector<ProblemDomain>     vectDomain;
  Vector<Real>              vectDx;
  getDomainsAndDxes(vectDomain, vectDx, a_params);

  int numlevels = a_params.numLevels;

  ParmParse pp;
  // grid generation parameters

  vectGrids.resize(numlevels);
  bool readInGrids = false;
  if (pp.contains("read_in_grids"))
    {
      pp.get("read_in_grids", readInGrids);
    }

  if (readInGrids)
    {

      ProblemDomain levDomain = a_params.coarsestDomain;
      for (int ilev = 0; ilev < a_params.numLevels; ilev++)
        {
          Vector<Box>   boxes;
          char boxCountVar[100];
          int boxCount;
          sprintf(boxCountVar, "level_%d_box_count", ilev);
          pp.get(boxCountVar, boxCount);
          boxes.resize(boxCount);
          for (int ibox = 0; ibox < boxCount; ibox++)
            {
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
              if (!levDomain.contains(boxes[ibox]))
                {
                  MayDay::Error("box outside of domain");
                }
            }
          //check to see if level 0 domain is covered
          if (ilev == 0)
            {
              IntVectSet ivDom(levDomain.domainBox());
              for (int ibox = 0; ibox < boxes.size(); ibox++)
                {
                  ivDom -= boxes[ibox];
                }
              if (!ivDom.isEmpty())
                {
                  MayDay::Error("level 0 boxes must cover the domain");
                }
            }
          Vector<int>  proc(a_params.numLevels);
          LoadBalance(proc,boxes);
          vectGrids[ilev] = DisjointBoxLayout(boxes, proc, levDomain);
          levDomain.refine(a_params.refRatio[ilev]);
        }

    }
  else // no read_grids
    {
      pout() << "tagging on gradient of RHS" << endl;
      int maxLevel = numlevels-1;
      Vector<Vector<Box> > newBoxes(numlevels);
      Vector<Vector<Box> > oldBoxes(numlevels);

      // determine grids dynamically, based on grad(RHS)
      // will need temp storage for RHS
      Vector<LevelData<FArrayBox>* > vectRHS(maxLevel+1,NULL);
      int ncomps = 1;

      // define base level first
      Vector< Vector<int> > procAssign(maxLevel+1);
      domainSplit(vectDomain[0], oldBoxes[0], a_params.maxGridSize, a_params.blockFactor);
      procAssign[0].resize(oldBoxes[0].size());
      LoadBalance(procAssign[0],oldBoxes[0]);

      vectGrids[0].define(oldBoxes[0],procAssign[0],vectDomain[0]);

      vectRHS[0] = new LevelData<FArrayBox>(vectGrids[0], ncomps,
                                            IntVect::Zero);

      int topLevel = 0;

      bool moreLevels = (maxLevel > 0);

      int nesting_radius = 2;
      // create grid generation object
      BRMeshRefine meshrefine(vectDomain[0], a_params.refRatio,
                              a_params.fillRatio,
                              a_params.blockFactor, nesting_radius,
                              a_params.maxGridSize);

      while (moreLevels)
        {
          // default is moreLevels = false
          // (only repeat loop in the case where a new level
          // is generated which is still less than maxLevel)
          moreLevels = false;

          int baseLevel = 0;
          int oldTopLevel = topLevel;

          // now initialize RHS for this existing hierarchy

          for (int level=0; level<=topLevel; level++)
            {
              RealVect dxLevel = vectDx[level]*RealVect::Unit;

              LevelData<FArrayBox> *DummyChi; // REMEMBER TO FIX THIS 
              LevelData<FArrayBox> *DummyPsi; // REMEMBER TO FIX THIS 
      
              DummyChi = new LevelData<FArrayBox>(vectGrids[level], ncomps,
                                            IntVect::Unit);
              DummyPsi = new LevelData<FArrayBox>(vectGrids[level], ncomps,
                                            IntVect::Unit);

              set_initial_chi(*DummyChi, dxLevel, a_params);
              set_initial_psi(*DummyPsi, dxLevel, a_params);
              setRHS(*vectRHS[level], *DummyChi, *DummyPsi, dxLevel, a_params);
            }

          Vector<IntVectSet> tagVect(topLevel+1);
          int tags_grow = 1;
          tagCells(vectRHS, tagVect, vectDx, vectDomain,
                   a_params.refineThresh,
                   tags_grow, baseLevel, topLevel+1);

          int new_finest = meshrefine.regrid(newBoxes, tagVect,
                                             baseLevel,
                                             topLevel, oldBoxes);

          if (new_finest > topLevel)
            {
              topLevel++;
            }

          oldBoxes = newBoxes;

          //  no need to do this for the base level (already done)
          for (int lev=1; lev<= topLevel; lev++)
            {
              // do load balancing
              procAssign[lev].resize(newBoxes[lev].size());
              LoadBalance(procAssign[lev], newBoxes[lev]);
              const DisjointBoxLayout newDBL(newBoxes[lev], procAssign[lev],
                                             vectDomain[lev]);
              vectGrids[lev] = newDBL;
              delete vectRHS[lev];
              vectRHS[lev] = new LevelData<FArrayBox>(vectGrids[lev], ncomps,
                                                      IntVect::Zero);
            } // end loop over levels for initialization

          // figure out whether we need another pass through grid generation
          if ((topLevel<maxLevel) && (topLevel > oldTopLevel))
            moreLevels = true;

        } // end while moreLevels loop
      // clean up temp storage
      for (int ilev=0; ilev <vectRHS.size(); ilev++)
        {
          if (vectRHS[ilev] != NULL)
            {
              delete vectRHS[ilev];
              vectRHS[ilev] = NULL;
            }
        }
    }
  return 0;
}

// EUGENE TODO : make it more general
void set_initial_chi(LevelData<FArrayBox>&    a_chi,

                     const RealVect&          a_dx,
                     const VCPoissonParameters& a_params)
{
  CH_assert(a_chi.nComp() == 1);
  int comp=0;
//  const RealVect&    trig = getTrigRV();

  for (DataIterator dit = a_chi.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox& thisCHIbox = a_chi[dit()];
    Box b = thisCHIbox.box ();
    BoxIterator bit (b);
    for (bit.begin (); bit.ok (); ++bit)
    {
      IntVect iv = bit ();
      for (int comp =0 ; comp < a_chi.nComp () ; ++comp)
      { 
        thisCHIbox (iv,comp) = a_params.initial_chi;
      }
    }
  }
} // end set_initial_chi

// EUGENE : temporary
void set_initial_psi(LevelData<FArrayBox>&    a_psi,
                     const RealVect&          a_dx,
                     const VCPoissonParameters& a_params)
{
  CH_assert(a_psi.nComp() == 1);
  int comp=0;
//  const RealVect&    trig = getTrigRV();

  for (DataIterator dit = a_psi.dataIterator(); dit.ok(); ++dit)
  {
    FArrayBox& thisPSIbox = a_psi[dit()];
    Box b = thisPSIbox.box ();
    BoxIterator bit (b);
    for (bit.begin (); bit.ok (); ++bit)
    {
      IntVect iv = bit ();
      for (int comp =0 ; comp < a_psi.nComp () ; ++comp)
      { 
        thisPSIbox (iv,comp) = a_params.initial_psi;
      }
    }
  }
} // end set_initial_psi

/********/
void setRHS(LevelData<FArrayBox>&    a_rhs,
            LevelData<FArrayBox>&    a_chi,
            LevelData<FArrayBox>&    a_psi,
            const RealVect&          a_dx,
            const VCPoissonParameters& a_params)
{
  CH_assert(a_rhs.nComp() == 1);

  int comp = 0;

  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& thisRHS = a_rhs[dit()];
      Box thisBox = thisRHS.box();

      if (a_params.probtype == 0) // original
        {
          const RealVect&  trig = getTrigRV();
          FORT_GETLOFCHI(CHF_FRA1(thisRHS,comp),
                        CHF_CONST_REALVECT(trig),
                        CHF_CONST_REALVECT(a_dx),
                        CHF_CONST_REALVECT(a_params.probLo),
                        CHF_CONST_REALVECT(a_params.probHi),
                        CHF_CONST_REAL(a_params.alpha),
                        CHF_CONST_REAL(a_params.beta),
                        CHF_CONST_REAL(a_params.gamma),
                        CHF_BOX(thisBox));
        } 
      else if (a_params.probtype == 1) // gaussians
        {
          // rhs is cell-centered...
          RealVect ccOffset = 0.5*a_dx*RealVect::Unit;
          int numGaussians = 3;
          Vector<RealVect> center(numGaussians,RealVect::Zero);
          Vector<Real> scale(numGaussians, 1.0);
          Vector<Real> strength(numGaussians, 1.0);
          
          for (int n=0; n<numGaussians; n++)
            {
              if (n==0)
                {
                  strength[0] = 1.0;
                  scale[0] = 1.0e-2;
                  center[0] = 0.25*RealVect::Unit;
                }
              else if (n == 1)
                {
                  strength[1] = 3.0;
                  scale[1] = 1.0e-2;
                  center[1] = RealVect(D_DECL(0.5,0.75, 0.75));
                }
              else if (n == 2)
                {
                  strength[2] = 2.0;
                  scale[2] = 1.0e-2;
                  center[2] = RealVect(D_DECL(0.75,0.5, 0.5));
                }
              else
                {
                  MayDay::Error("too many Gaussian sources attempted");
                }
            }

          thisRHS.setVal(0,0);

          BoxIterator bit(thisRHS.box());
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc(iv);
              loc *= a_dx;
              loc += ccOffset;

              for (int n=0; n<numGaussians; n++)
                {
                  RealVect dist = loc - center[n];
                  Real radSqr = D_TERM(dist[0]*dist[0],
                                        +dist[1]*dist[1],
                                        +dist[2]*dist[2]);

                  Real val = strength[n]*exp(-radSqr/scale[n]);
                  thisRHS(iv,0) += val;
                }
            }
        }
      else if (a_params.probtype == 2) // scalar field
        {
          //
          // first set up the scalar field config
          //
          FArrayBox& thisPsi= a_psi[dit()];
          
          // rhs is cell-centered...
          RealVect ccOffset = 0.5*a_dx*RealVect::Unit;
          
          // put in a gaussian for now
          int numGaussians = 1;
          Vector<RealVect> center(numGaussians,RealVect::Zero);
          Vector<Real> scale(numGaussians, 1.0);
          Vector<Real> strength(numGaussians, 1.0);
          BoxIterator bit(thisRHS.box());

          // one central gaussian 
          strength[0] = 1.0;
          scale[0] = 1.0e-2; // variance
          center[0] = 0.5*RealVect::Unit;

          thisPsi.setVal(0,0);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc(iv);
              loc *= a_dx;
              loc += ccOffset;

              for (int n=0; n<numGaussians; n++)
                {
                  RealVect dist = loc - center[n];
                  Real radSqr = D_TERM(dist[0]*dist[0],
                                        +dist[1]*dist[1],
                                        +dist[2]*dist[2]);

                  Real val = strength[n]*exp(-radSqr/scale[n]);
                  thisPsi(iv,0) += val;
                }
            }

          // end set up gaussian
        } // end probtype 2
      else if (a_params.probtype == 3) // fixed energy density
        {

          // Total equation is
          // nabla^2 chi = (1/8)chi^5 (2/3K^2 - 16pi G rho(x))
          // 
          // linearized chi^n = (1-n)chi_0^n + n chi_0^(n-1)chi
          //
          // linearized equation is then
          //
          // L chi_0 = RHS
          //
          // L = nabla^2 - (5/8)chi_0^4 * M(rho,K)
          // RHS = -(1/2)M(rho,K) chi_0^5
          // M(rho,K) = (2/3K^2 - 16piG rho)
          // K = constant 

          // cell-centered
          RealVect ccOffset = 0.5*a_dx*RealVect::Unit;

          FArrayBox& thisCHI = a_chi[dit()]; 
          Box thisBox = thisRHS.box();

          thisRHS.setVal(0,0); // set everything to zero
         
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              IntVect iv = bit();
              RealVect loc(iv);
              loc *= a_dx;
              loc += ccOffset;

              // set up M(rho, K)

              Real M = M_value(loc, a_params);

              Real chi_0 = thisCHI(iv,0);

              Real val = (-0.5) * M * (chi_0 * chi_0 * chi_0 * chi_0 * chi_0); // -(1/2)M chi_0^5

              thisRHS(iv,0) += val;

            }

        } // end prob_type
     }

} // end set_RHS


extern
AMRLevelOpFactory<LevelData<FArrayBox> >*
defineOperatorFactory(
                      const Vector<DisjointBoxLayout>&             a_grids,
                      const Vector<ProblemDomain>&                 a_vectDomain,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_aCoef,
                      Vector<RefCountedPtr<LevelData<FluxBox> > >& a_bCoef,
                      Vector<RefCountedPtr<LevelData<FArrayBox> > >& a_cCoef,
                      const VCPoissonParameters&                     a_params)
{
  ParmParse pp2;

  VCAMRPoissonOp2Factory* opFactory = new VCAMRPoissonOp2Factory;

  opFactory->define(a_params.coarsestDomain,
                    a_grids,
                    a_params.refRatio,
                    a_params.coarsestDx,
                    &ParseBC,
                    a_params.alpha,
                    a_aCoef,
                    a_params.beta,
                    a_bCoef,
                    a_params.gamma,
                    a_cCoef);

  if (a_params.coefficient_average_type >= 0)
    {
      opFactory->m_coefficient_average_type
        = a_params.coefficient_average_type;
    }

  return (AMRLevelOpFactory<LevelData<FArrayBox> >*) opFactory;
}

/*
  Set grid hierarchy from input file
 */
void getDomainsAndDxes(  Vector<ProblemDomain>&     vectDomain,
                         Vector<Real>&              vectDx,
                         VCPoissonParameters&         a_params)
{

  vectDomain.resize(a_params.numLevels);
  vectDx.resize(    a_params.numLevels);
  vectDx[0] = a_params.coarsestDx;
  for (int ilev = 1; ilev < a_params.numLevels; ilev++)
    {
      vectDx[ilev] = vectDx[ilev-1]/a_params.refRatio[ilev-1];
    }


  vectDomain[0] = a_params.coarsestDomain;
  for (int ilev = 1;ilev < a_params.numLevels; ilev++)
    {
      vectDomain[ilev] = refine(vectDomain[ilev-1],a_params.refRatio[ilev-1]);
    }
}


/*
  tag cells for refinement based on magnitude(RHS)
*/
void
tagCells(Vector<LevelData<FArrayBox>* >& vectRHS,
         Vector<IntVectSet>& tagVect,
         Vector<Real>& vectDx,
         Vector<ProblemDomain>& vectDomain,
         const Real refine_thresh,
         const int tags_grow,
         const int baseLevel,
         int numLevels)
{
  for (int lev=baseLevel; lev!= numLevels; lev++)
    {
      IntVectSet local_tags;
      LevelData<FArrayBox> & levelRhs = *vectRHS[lev];
      DisjointBoxLayout level_domain = levelRhs.getBoxes();
      DataIterator dit = levelRhs.dataIterator();

      Real maxRHS = 0;

      maxRHS = norm(levelRhs, levelRhs.interval(), 0);

      Real tagVal = maxRHS * refine_thresh;

      // now loop through grids and tag cells where RHS > tagVal
      for (dit.reset(); dit.ok(); ++dit)
        {
          const Box thisBox = level_domain.get(dit());
          const FArrayBox& thisRhs = levelRhs[dit()];
          BoxIterator bit(thisBox);
          for (bit.begin(); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (abs(thisRhs(iv)) >= tagVal)
                local_tags |= iv;
            }
        } // end loop over grids on this level

      local_tags.grow(tags_grow);
      const Box& domainBox = vectDomain[lev].domainBox();
      local_tags &= domainBox;

      tagVect[lev] = local_tags;

    } // end loop over levels
}



void TrigValueNeum(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  RealVect gradChi;
  FORT_GETGRADCHIPOINT(CHF_REALVECT(gradChi),
                       CHF_CONST_REALVECT(trig),
                       CHF_CONST_REALVECT(xval));

  a_values[0] = gradChi[*dir];
}
void TrigValueDiri(Real* pos,
                   int* dir,
                   Side::LoHiSide* side,
                   Real* a_values)
{
  RealVect& trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      xval[idir] = pos[idir];
    }
  Real value;
  FORT_GETCHIPOINT(CHF_REAL(value),
                   CHF_CONST_REALVECT(trig),
                   CHF_CONST_REALVECT(xval));
  a_values[0] = value;
}

// M(K, rho) = 2/3K^2 - 16piG rho
// for now we put rho to be a gaussian with center (0.5,0.5,0.5)
Real M_value(RealVect&          loc, 
             const VCPoissonParameters& a_params)
{
  Real rho;

  // Gaussian for dirichlet BC
  if (a_params.rho_type == 0)  // central gaussian rho = strength * Exp[-r^2/scale]
  {
    RealVect center;
    center = RealVect(D_DECL(a_params.rho_center_1[0], a_params.rho_center_1[1], a_params.rho_center_1[2])); // put it in the middle for now


    RealVect dist = loc - center;
    Real radSqr = D_TERM(dist[0]*dist[0],
                        +dist[1]*dist[1],
                        +dist[2]*dist[2]);
    
    rho = a_params.rho_strength * exp(-radSqr/a_params.rho_scale);
  }
  else if (a_params.rho_type == 1)  // wave modes (for periodic) rho = Sum_i a_i sin{2pi k_i x_i) + baseline (to prevent negative rho)
   
  {
    
    rho = a_params.rho_amplitude[0] * sin(a_params.rho_k[0] * 2.0 * PI * loc[0]) 
        + a_params.rho_amplitude[1] * sin(a_params.rho_k[1] * 2.0 * PI * loc[1]) 
        + a_params.rho_amplitude[2] * sin(a_params.rho_k[2] * 2.0 * PI * loc[2]) 
        + a_params.rho_baseline;
  }
  else if (a_params.rho_type == 2) // 3 gaussians
  {
    RealVect center;

    center = RealVect(D_DECL(a_params.rho_center_1[0], a_params.rho_center_1[1], a_params.rho_center_1[2])); // put it in the middle for now

    RealVect dist = loc - center;
    Real radSqr = D_TERM(dist[0]*dist[0],
                        +dist[1]*dist[1],
                        +dist[2]*dist[2]);
    
    rho = a_params.rho_strength * exp(-radSqr/a_params.rho_scale);

    center = RealVect(D_DECL(a_params.rho_center_2[0], a_params.rho_center_2[1], a_params.rho_center_2[2])); // put it in the middle for now

    dist = loc - center;
    radSqr = D_TERM(dist[0]*dist[0],
                    +dist[1]*dist[1],
                    +dist[2]*dist[2]);

    rho += a_params.rho_strength * exp(-radSqr/a_params.rho_scale);

    center = RealVect(D_DECL(a_params.rho_center_3[0], a_params.rho_center_3[1], a_params.rho_center_3[2])); // put it in the middle for now

    dist = loc - center;
    radSqr = D_TERM(dist[0]*dist[0],
                    +dist[1]*dist[1],
                    +dist[2]*dist[2]);

    rho += a_params.rho_strength * exp(-radSqr/a_params.rho_scale);
  }
  else
  {
    rho = 0.0;
    // probably should put some error handling code here
  }

  return ((2.0/3.0)*(a_params.constant_K * a_params.constant_K)- 2.0*(a_params.kappa_sq)*rho);

}







