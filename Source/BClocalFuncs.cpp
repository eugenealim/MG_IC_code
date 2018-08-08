#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BClocalFuncs.H"
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

// Global BCRS definitions
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
            //} else if (GlobalBCRS::s_bcLo[i] == 2) {
            // if (!GlobalBCRS::s_printedThatLo[i]) {
            //  GlobalBCRS::s_printedThatLo[i] = true;
            //  pout() << "trig neum bcs lo for direction " << i << endl;
            //}
            // NeumBC(a_state, valid, a_dx, a_homogeneous, TrigValueNeum, i,
            //       Side::Lo);
          } else if (GlobalBCRS::s_bcLo[i] == 0) {
            if (!GlobalBCRS::s_printedThatLo[i]) {
              GlobalBCRS::s_printedThatLo[i] = true;
              pout() << "const diri bcs lo for direction " << i << endl;
            }
            DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValue, i,
                   Side::Lo);

            //} else if (GlobalBCRS::s_bcLo[i] == 3) {
            // if (!GlobalBCRS::s_printedThatLo[i]) {
            //  GlobalBCRS::s_printedThatLo[i] = true;
            //  pout() << "trig diri bcs lo for direction " << i << endl;
            //}
            // DiriBC(a_state, valid, a_dx, a_homogeneous, TrigValueDiri, i,
            //       Side::Lo);

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
            //} else if (GlobalBCRS::s_bcHi[i] == 2) {
            //  if (!GlobalBCRS::s_printedThatHi[i]) {
            //    GlobalBCRS::s_printedThatHi[i] = true;
            //    pout() << "trig neum bcs hi for direction " << i << endl;
            //  }
            //  NeumBC(a_state, valid, a_dx, a_homogeneous, TrigValueNeum, i,
            //         Side::Hi);
          } else if (GlobalBCRS::s_bcHi[i] == 0) {
            if (!GlobalBCRS::s_printedThatHi[i]) {
              GlobalBCRS::s_printedThatHi[i] = true;
              pout() << "const diri bcs hi for direction " << i << endl;
            }
            DiriBC(a_state, valid, a_dx, a_homogeneous, ParseValue, i,
                   Side::Hi);
            //} else if (GlobalBCRS::s_bcHi[i] == 3) {
            //  if (!GlobalBCRS::s_printedThatHi[i]) {
            //    GlobalBCRS::s_printedThatHi[i] = true;
            //    pout() << "trig diri bcs hi for direction " << i << endl;
            //  }
            //  DiriBC(a_state, valid, a_dx, a_homogeneous, TrigValueDiri, i,
            //         Side::Hi);
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
