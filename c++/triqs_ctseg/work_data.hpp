// Copyright (c) 2022-2024 Simons Foundation
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You may obtain a copy of the License at
//     https://www.gnu.org/licenses/gpl-3.0.txt
//
// Authors: Nikita Kavokine, Olivier Parcollet, Nils Wentzell

#pragma once
#include <mpi/mpi.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/det_manip.hpp>

#include "params.hpp"
#include "inputs.hpp"
#include "util.hpp"
#include "dets.hpp"

namespace triqs_ctseg {

  struct work_data_t {

    work_data_t(params_t const &p, inputs_t const &inputs, mpi::communicator c);

    nda::matrix<double> U;  // Density-density interaction: U_ab n_a n_b
    gf_struct_t gf_struct;  // gf_struct of the Green's function (input copied)
    int n_color;            // Number of colors
    nda::vector<double> mu; // Chemical potential per color

    bool has_Delta     = false; // There is a non-zero hybridization term
    bool has_Dt        = false; // There is a non-zero dynamical nn interaction
    bool has_Jperp     = false; // There is a non-zero Jperp interaction
    bool rot_inv       = true;  // The spin-spin interaction is rotationally invariant (matters for F(tau) measure)
    bool minus_sign    = false; // Has a move ever produced a negative sign?
    bool offdiag_Delta = false; // Does Delta(tau) have blocks of size larger than 1?

    // Dynamical and spin-spin interaction kernels
    gf<imtime> D0t, K, Kprime, Jperp, Kprime_spin;

    // Hybridization function
    block_gf<imtime, matrix_real_valued> Delta;

    // The determinants
    // Vector of the det_manip objects, one per block of the input Delta(tau). See dets.hpp
    std::vector<det_t> dets;

    // Color to (block, idx) conversion tables
    std::vector<long> block_number;   // block numbers corresponding to colors
    std::vector<long> index_in_block; // index in block of a given color

    // Find color corresponding to (block, idx)
    int block_to_color(int block, int idx) const;

    // Find block of color
    long find_block_number(int color) const;

    // Find index of color in its block
    long find_index_in_block(int color) const;
  };

  // Additional sign of the trace (computed from dets).
  double trace_sign(work_data_t const &wdata);

  // Functions for checking if a time is already in det.
  bool c_in_det(tau_t const &tau, det_t const &D);

  bool cdag_in_det(tau_t const &tau, det_t const &D);

} // namespace triqs_ctseg
