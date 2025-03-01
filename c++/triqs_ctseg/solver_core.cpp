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

#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

#include "solver_core.hpp"
#include "work_data.hpp"
#include "configuration.hpp"
#include "measures.hpp"
#include "moves.hpp"
#include "logs.hpp"

namespace triqs_ctseg {

  // ---------------------------------------------------------------------------

  solver_core::solver_core(constr_params_t const &p) : constr_params(p) {

    beta = p.beta;
    tau_t::set_beta(beta);

    inputs.Delta  = block_gf<imtime>({beta, Fermion, p.n_tau}, p.gf_struct);
    inputs.D0t    = make_block2_gf<imtime>({beta, Boson, p.n_tau_bosonic}, p.gf_struct);
    inputs.Jperpt = gf<imtime>({beta, Boson, p.n_tau_bosonic}, {1, 1});

    inputs.Delta()  = 0;
    inputs.D0t()    = 0;
    inputs.Jperpt() = 0;
  };

  // ---------------------------------------------------------------------------

  void solver_core::solve(solve_params_t const &solve_params_input) {

    // http://patorjk.com/software/taag/#p=display&f=Calvin%20S&t=TRIQS%20ctseg
    if (c.rank() == 0)
      std::cout << "\n"
                   "╔╦╗╦═╗╦╔═╗ ╔═╗  ┌─┐┌┬┐┌─┐┌─┐┌─┐\n"
                   " ║ ╠╦╝║║═╬╗╚═╗  │   │ └─┐├┤ │ ┬\n"
                   " ╩ ╩╚═╩╚═╝╚╚═╝  └─┘ ┴ └─┘└─┘└─┘\n";

    // ................ Parameters .................
    // Store the solve_params
    solve_params = solve_params_input;
    // Set tau mesh parameters for results to default if not supplied
    if (solve_params_input.n_tau_G == 0) solve_params.n_tau_G = constr_params.n_tau;
    if (solve_params_input.n_tau_G == 0) solve_params.n_tau_chi2 = constr_params.n_tau_bosonic;
    // Merge constr_params and solve_params
    params_t p(constr_params, solve_params);

    // ................   Work data & Configuration  ...................

    // Initialize work data
    work_data_t wdata{p, inputs, c};
    // Initialize configuration
    configuration_t config{wdata.n_color};
    // Start from a non-empty configuration when Delta(tau) = 0
    if (not wdata.has_Delta) { config.seglists[0].push_back(segment_t::full_line()); }

    // ................   QMC  ...................

    auto CTQMC = triqs::mc_tools::mc_generic<double>(p.random_name, p.random_seed, p.verbosity);

    // Initialize moves
    if (wdata.has_Delta) {
      if (p.move_insert_segment) CTQMC.add_move(moves::insert_segment{wdata, config, CTQMC.get_rng()}, "insert");
      if (p.move_remove_segment) CTQMC.add_move(moves::remove_segment{wdata, config, CTQMC.get_rng()}, "remove");
      if (p.move_double_insert_segment) CTQMC.add_move(moves::double_insert_segment{wdata, config, CTQMC.get_rng()}, "double insert");
      if (p.move_double_remove_segment) CTQMC.add_move(moves::double_remove_segment{wdata, config, CTQMC.get_rng()}, "double remove");
      if (p.move_move_segment) CTQMC.add_move(moves::move_segment{wdata, config, CTQMC.get_rng()}, "move");
      if (p.move_split_segment) CTQMC.add_move(moves::split_segment{wdata, config, CTQMC.get_rng()}, "split");
      if (p.move_regroup_segment) CTQMC.add_move(moves::regroup_segment{wdata, config, CTQMC.get_rng()}, "regroup");
    }

    if (wdata.has_Jperp) {
      if (p.move_insert_spin_segment)
        CTQMC.add_move(moves::insert_spin_segment{wdata, config, CTQMC.get_rng()}, "spin insert");

      if (p.move_remove_spin_segment)
        CTQMC.add_move(moves::remove_spin_segment{wdata, config, CTQMC.get_rng()}, "spin remove");
    }

    if (wdata.has_Jperp and wdata.has_Delta) {
      if (p.move_split_spin_segment)
        CTQMC.add_move(moves::split_spin_segment{wdata, config, CTQMC.get_rng()}, "spin split");

      if (p.move_regroup_spin_segment)
        CTQMC.add_move(moves::regroup_spin_segment{wdata, config, CTQMC.get_rng()}, "spin regroup");
    }

    if (wdata.has_Jperp) {
      if (p.move_swap_spin_lines) CTQMC.add_move(moves::swap_spin_lines{wdata, config, CTQMC.get_rng()}, "spin swap");
    }

    // Initialize measurements
    if (p.measure_G_tau) CTQMC.add_measure(measures::G_F_tau{p, wdata, config, results}, "G(tau)/F(tau)");
    if (p.measure_densities) CTQMC.add_measure(measures::densities{p, wdata, config, results}, "Densities");
    if (p.measure_average_sign) CTQMC.add_measure(measures::average_sign{p, wdata, config, results}, "Average Sign");
    if (p.measure_nn_static) CTQMC.add_measure(measures::nn_static{p, wdata, config, results}, "<nn>");
    if (p.measure_nn_tau) CTQMC.add_measure(measures::nn_tau{p, wdata, config, results}, "<n(tau)n(0)>");
    if (p.measure_Sperp_tau) CTQMC.add_measure(measures::Sperp_tau{p, wdata, config, results}, "<S_x(tau)S_x(0)>");
    if (p.measure_pert_order) {
      if (wdata.has_Delta) {
        CTQMC.add_measure(measures::pert_order{[&]() { return config.Delta_order(); }, results.pert_order_Delta,
                                               results.average_order_Delta},
                          "Perturbation order Delta");
      }
      if (wdata.has_Jperp) {
        CTQMC.add_measure(measures::pert_order{[&]() { return config.Jperp_order(); }, results.pert_order_Jperp,
                                               results.average_order_Jperp},
                          "Perturbation order Jperp");
      }
    }
    if (p.measure_state_hist) CTQMC.add_measure(measures::state_hist{p, wdata, config, results}, "State histograms");

    // Run and collect results
    CTQMC.warmup_and_accumulate(p.n_warmup_cycles, p.n_cycles, p.length_cycle,
                                triqs::utility::clock_callback(p.max_time));
    CTQMC.collect_results(c);

    // Report sign and average order
    if (c.rank() == 0) {
      spdlog::info("Average sign: {}", results.average_sign);
      if (results.average_order_Delta)
        spdlog::info("Average perturbation order in Delta: {:.3f}", results.average_order_Delta.value());
      if (results.average_order_Jperp)
        spdlog::info("Average perturbation order in Jperp: {:.3f}", results.average_order_Jperp.value());
    }

  } // solve

  // ----------------- Save to h5 file -----------------------

#define STR(x) #x
#define STRINGIZE(x) STR(x)

  // Function that writes the solver_core to hdf5 file
  void h5_write(h5::group h5group, std::string subgroup_name, solver_core const &s) {
    auto grp = h5group.create_group(subgroup_name);
    h5_write_attribute(grp, "Format", solver_core::hdf5_format());
    h5_write_attribute(grp, "TRIQS_GIT_HASH", std::string(STRINGIZE(TRIQS_GIT_HASH)));
    h5_write_attribute(grp, "CTSEG_GIT_HASH", std::string(STRINGIZE(CTSEG_GIT_HASH)));
    h5_write(grp, "constr_params", s.constr_params);
    h5_write(grp, "solve_params", s.solve_params);
    h5_write(grp, "inputs", s.inputs);
    h5_write(grp, "results", s.results);
  }

  // Function that reads all containers in hdf5 file
  solver_core solver_core::h5_read_construct(h5::group h5group, std::string subgroup_name) {
    auto grp           = h5group.open_group(subgroup_name);
    auto constr_params = h5_read<constr_params_t>(grp, "constr_params");
    auto s             = solver_core{constr_params};
    h5_read(grp, "solve_params", s.solve_params);
    h5_read(grp, "inputs", s.inputs);
    h5_read(grp, "results", s.results);
    return s;
  }

} // namespace triqs_ctseg
