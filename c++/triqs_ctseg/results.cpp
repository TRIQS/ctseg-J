#include "./results.hpp"

namespace triqs_ctseg {

  void h5_write(h5::group h5group, std::string subgroup_name, results_t const &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.create_group(subgroup_name);

    h5_write(grp, "G_tau", c.G_tau);
    h5_write(grp, "sign", c.sign);
    h5_write(grp, "F_tau", c.F_tau);
    h5_write(grp, "nn_tau", c.nn_tau);
    h5_write(grp, "sperp_tau", c.sperp_tau);
    h5_write(grp, "nn_static", c.nn_static);
    h5_write(grp, "densities", c.densities);
    h5_write(grp, "pert_order_histo_Delta", c.pert_order_histo_Delta);
    h5_write(grp, "pert_order_histo_Jperp", c.pert_order_histo_Jperp);
    h5_write(grp, "state_hist", c.state_hist);
  }

  //------------------------------------

  void h5_read(h5::group h5group, std::string subgroup_name, results_t &c) {

    h5::group grp = subgroup_name.empty() ? h5group : h5group.open_group(subgroup_name);

    h5_read(grp, "G_tau", c.G_tau);
    h5_read(grp, "sign", c.sign);
    h5_read(grp, "F_tau", c.F_tau);
    h5_read(grp, "nn_tau", c.nn_tau);
    h5_read(grp, "sperp_tau", c.sperp_tau);
    h5_read(grp, "nn_static", c.nn_static);
    h5_read(grp, "densities", c.densities);
    h5_read(grp, "pert_order_histo_Delta", c.pert_order_histo_Delta);
    h5_read(grp, "pert_order_histo_Jperp", c.pert_order_histo_Jperp);
    h5_read(grp, "state_hist", c.state_hist);
  }

} // namespace triqs_ctseg
