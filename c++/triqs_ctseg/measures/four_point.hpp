// Copyright (c) 2024 Simons Foundation
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
// Authors: Hao Lu

#pragma once
#include "../configuration.hpp"
#include "../work_data.hpp"
#include "../results.hpp"

namespace triqs_ctseg::measures {

  struct four_point {

    work_data_t const &wdata;
    configuration_t const &config;
    results_t &results;
    double beta;
    double w_ini, w_inc;
    bool measure_g3w;
    bool measure_f3w;
    int n_w_fermionic;
    int n_w_bosonic;
    std::vector<std::string> block_names;
    std::vector<array<dcomplex, 4>> Mw_vector, nMw_vector;
    vector<dcomplex> y_exp_ini, y_exp_inc, x_exp_ini, x_exp_inc;
    vector<int> y_inner_index, x_inner_index;

    block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> g3w;
    block_gf<prod<imfreq, imfreq, imfreq>, tensor_valued<4>> f3w;

    double Z = 0;

    four_point(params_t const &params, work_data_t const &wdata, configuration_t const &config, results_t &results);

    void accumulate(double s);
    void collect_results(mpi::communicator const &c);
    double fprefactor(long const &block, std::pair<tau_t, long> const &y);
    
    std::vector<array<dcomplex, 4>> compute_Mw(bool is_nMw);

    dcomplex Mw(long const &block, int const &i, int const &j, int const &n1, int const &n2) {
      return Mw_vector[block](i, j, n1 + n_w_fermionic + n_w_bosonic - 1, n2 + n_w_fermionic + n_w_bosonic - 1);
    }

    dcomplex nMw(long const &block, int const &i, int const &j, int const &n1, int const &n2) {
      return nMw_vector[block](i, j, n1 + n_w_fermionic + n_w_bosonic - 1, n2 + n_w_fermionic + n_w_bosonic - 1);
    }

  };

} // namespace triqs_ctseg::measures