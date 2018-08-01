#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VClocalFuncs.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "ChiFunctionsF_F.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

std::vector<bool> GlobalBCRS::s_printedThatLo =
    std::vector<bool>(SpaceDim, false);
std::vector<bool> GlobalBCRS::s_printedThatHi =
    std::vector<bool>(SpaceDim, false);
std::vector<int> GlobalBCRS::s_bcLo = std::vector<int>();
std::vector<int> GlobalBCRS::s_bcHi = std::vector<int>();
RealVect GlobalBCRS::s_trigvec = RealVect::Zero;
bool GlobalBCRS::s_areBCsParsed = false;
bool GlobalBCRS::s_valueParsed = false;
bool GlobalBCRS::s_trigParsed = false;

// BCValueHolder class, which is a pointer to a void-type function with the 4
// arguements given pos [x,y,z] position on center of cell edge int dir
// direction, x being 0 int side -1 for low, +1 = high, fill in the a_values
// array

void ParseValue(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values) {
  ParmParse pp;
  Real bcVal;
  pp.get("bc_value", bcVal);
  a_values[0] = bcVal;
}

void ParseBC(FArrayBox &a_state, const Box &a_valid,
             const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous) {
  if (!a_domain.domainBox().contains(a_state.box())) {

    if (!GlobalBCRS::s_areBCsParsed) {
      ParmParse pp;
      pp.getarr("bc_lo", GlobalBCRS::s_bcLo, 0, SpaceDim);
      pp.getarr("bc_hi", GlobalBCRS::s_bcHi, 0, SpaceDim);
      GlobalBCRS::s_areBCsParsed = true;
    }

    Box valid = a_valid;

    for (int i = 0; i < CH_SPACEDIM; ++i) {

      // periodic?
      if (!a_domain.isPeriodic(i)) {
        // periodic?

        Box ghostBoxLo = adjCellBox(valid, i, Side::Lo, 1);
        Box ghostBoxHi = adjCellBox(valid, i, Side::Hi, 1);
        if (!a_domain.domainBox().contains(ghostBoxLo)) {
          if (GlobalBCRS::s_bcLo[i] == 1) {
            if (!GlobalBCRS::s_printedThatLo[i]) {
              GlobalBCRS::s_printedThatLo[i] = true;
              pout() << "const neum bcs lo for direction " << i << endl;
            }
            NeumBC(a_state, valid, a_dx, a_homogeneous,
                   ParseValue, // BCValueHolder class
                   i, Side::Lo);
          } else if (GlobalBCRS::s_bcLo[i] == 2) {
            if (!GlobalBCRS::s_printedThatLo[i]) {
              GlobalBCRS::s_printedThatLo[i] = true;
              pout() << "trig neum bcs lo for direction " << i << endl;
            }
            NeumBC(a_state, valid, a_dx, a_homogeneous, TrigValueNeum, i,
                   Side::Lo);
          } else if (GlobalBCRS::s_bcLo[i] == 0) {
            if (!GlobalBCRS::s_printedThatLo[i]) {
              GlobalBCRS::s_printedThatLo[i] = true;
              pout() << "const diri bcs lo for direction " << i << endl;
            }
            DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValue, i,
                   Side::Lo);

          } else if (GlobalBCRS::s_bcLo[i] == 3) {
            if (!GlobalBCRS::s_printedThatLo[i]) {
              GlobalBCRS::s_printedThatLo[i] = true;
              pout() << "trig diri bcs lo for direction " << i << endl;
            }
            DiriBC(a_state, valid, a_dx, a_homogeneous, TrigValueDiri, i,
                   Side::Lo);

          } else if (GlobalBCRS::s_bcLo[i] == 4) {
            if (!GlobalBCRS::s_printedThatLo[i]) {
              GlobalBCRS::s_printedThatLo[i] = true;
              pout() << "periodic bcs lo for direction " << i << endl;
            }

          } else {
            MayDay::Error("bogus bc flag lo");
          }
        }

        if (!a_domain.domainBox().contains(ghostBoxHi)) {
          if (GlobalBCRS::s_bcHi[i] == 1) {
            if (!GlobalBCRS::s_printedThatHi[i]) {
              GlobalBCRS::s_printedThatHi[i] = true;
              pout() << "const neum bcs hi for direction " << i << endl;
            }
            NeumBC(a_state, valid, a_dx, a_homogeneous, ParseValue, i,
                   Side::Hi);
          } else if (GlobalBCRS::s_bcHi[i] == 2) {
            if (!GlobalBCRS::s_printedThatHi[i]) {
              GlobalBCRS::s_printedThatHi[i] = true;
              pout() << "trig neum bcs hi for direction " << i << endl;
            }
            NeumBC(a_state, valid, a_dx, a_homogeneous, TrigValueNeum, i,
                   Side::Hi);
          } else if (GlobalBCRS::s_bcHi[i] == 0) {
            if (!GlobalBCRS::s_printedThatHi[i]) {
              GlobalBCRS::s_printedThatHi[i] = true;
              pout() << "const diri bcs hi for direction " << i << endl;
            }
            DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValue, i,
                   Side::Hi);
          } else if (GlobalBCRS::s_bcHi[i] == 3) {
            if (!GlobalBCRS::s_printedThatHi[i]) {
              GlobalBCRS::s_printedThatHi[i] = true;
              pout() << "trig diri bcs hi for direction " << i << endl;
            }
            DiriBC(a_state, valid, a_dx, a_homogeneous, TrigValueDiri, i,
                   Side::Hi);
          } else if (GlobalBCRS::s_bcHi[i] == 4) {
            if (!GlobalBCRS::s_printedThatHi[i]) {
              GlobalBCRS::s_printedThatHi[i] = true;
              pout() << "periodic bcs hi for direction " << i << endl;
            }
          } else {
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
RealVect &getTrigRV() {
  Real pi = 4. * atan(1.0);
  if (!GlobalBCRS::s_trigParsed) {
    GlobalBCRS::s_trigParsed = true;
    ParmParse pp;
    std::vector<Real> trigvec(SpaceDim);
    pp.getarr("trig", trigvec, 0, SpaceDim);

    for (int idir = 0; idir < SpaceDim; idir++) {
      GlobalBCRS::s_trigvec[idir] = pi * trigvec[idir];
    }
  }
  return GlobalBCRS::s_trigvec;
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

void TrigValueNeum(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values) {
  RealVect &trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++) {
    xval[idir] = pos[idir];
  }
  RealVect gradChi;
  FORT_GETGRADCHIPOINT(CHF_REALVECT(gradChi), CHF_CONST_REALVECT(trig),
                       CHF_CONST_REALVECT(xval));

  a_values[0] = gradChi[*dir];
}

void TrigValueDiri(Real *pos, int *dir, Side::LoHiSide *side, Real *a_values) {
  RealVect &trig = getTrigRV();
  RealVect xval;
  for (int idir = 0; idir < SpaceDim; idir++) {
    xval[idir] = pos[idir];
  }
  Real value;
  FORT_GETCHIPOINT(CHF_REAL(value), CHF_CONST_REALVECT(trig),
                   CHF_CONST_REALVECT(xval));
  a_values[0] = value;
}
