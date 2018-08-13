/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MULTIGRIDUSERVARIABLES_HPP
#define MULTIGRIDUSERVARIABLES_HPP

// assign an enum to each variable
enum {
  c_psi,

  c_A11_t0,
  c_A12_t0,
  c_A13_t0,
  c_A22_t0,
  c_A23_t0,
  c_A33_t0,

  c_phi_t0, // matter field

  NUM_MULTIGRID_VARS
};

namespace MultiGridUserVariables {
static constexpr char const *variable_names[NUM_MULTIGRID_VARS] = {
    "psi",

    "A11_0", "A12_0", "A13_0", "A22_0", "A23_0", "A33_0",

    "phi_0"};
}

#endif /* MULTIGRIDUSERVARIABLES_HPP */
