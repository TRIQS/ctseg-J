#include "regroup_spin_segment.hpp"
#include "../logs.hpp"
#include "triqs_ctseg/configuration.hpp"
#include <cmath>

namespace moves {

  double regroup_spin_segment::attempt() {

    prop_failed = false;

    LOG("\n =================== ATTEMPT REGROUP SPIN ================ \n");

    ln_trace_ratio = 0;
    det_ratio      = 1;
    prop_ratio     = 1;

    // ----------- Propose move in each color ----------
    std::tie(idx_c_up, idx_cdag_down, tau_up, prop_failed) = propose(0); // spin up
    if (prop_failed) return 0;
    std::tie(idx_c_down, idx_cdag_up, tau_down, prop_failed) = propose(1); // spin down
    if (prop_failed) return 0;

    // Correct for the dynamical interaction between the two operators that have been moved
    if (wdata.has_Dt) {
      auto &sl_up   = config.seglists[0];
      auto &sl_down = config.seglists[1];
      ln_trace_ratio -= real(wdata.K(double(tau_up - sl_down[idx_c_down].tau_c))(0, 1));
      ln_trace_ratio -= real(wdata.K(double(tau_down - sl_up[idx_c_up].tau_c))(0, 1));
      ln_trace_ratio += real(wdata.K(double(tau_down - tau_up))(0, 1));
    }

    double trace_ratio = std::exp(ln_trace_ratio);
    trace_ratio *= -real(wdata.Jperp(double(tau_up - tau_down))(0, 0)) / 2;
    prop_ratio /= double(config.Jperp_list.size()) + 1;

    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    double prod = trace_ratio * det_ratio * prop_ratio;
    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    return (std::isfinite(prod) ? prod : det_sign);
  }

  //--------------------------------------------------

  double regroup_spin_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    double initial_sign = config_sign(config, wdata.dets);

    // Update dets
    wdata.dets[0].complete_operation();
    wdata.dets[1].complete_operation();

    // Update the segments
    auto &sl_up               = config.seglists[0];
    auto &sl_down             = config.seglists[1];
    sl_up[idx_c_up].tau_c     = tau_up;
    sl_down[idx_c_down].tau_c = tau_down;

    // Update spin tags
    sl_up[idx_c_up].J_c           = true;
    sl_down[idx_cdag_down].J_cdag = true;
    sl_down[idx_c_down].J_c       = true;
    sl_up[idx_cdag_up].J_cdag     = true;

    // Add spin line
    config.Jperp_list.push_back(jperp_line_t{tau_up, tau_down});

    // Compute sign
    double final_sign = config_sign(config, wdata.dets);
    double sign_ratio = final_sign / initial_sign;
    LOG("Sign ratio is {}", sign_ratio);

    // Check invariant
    if constexpr (check_invariants or ctseg_debug) check_invariant(config, wdata.dets);
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);
    LOG("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------

  void regroup_spin_segment::reject() {

    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[0].reject_last_try();
    wdata.dets[1].reject_last_try();
  }

  //--------------------------------------------------

  std::tuple<long, long, tau_t, bool> regroup_spin_segment::propose(int color) {

    auto &sl         = config.seglists[color];
    int other_color  = 1 - color;
    auto &dsl        = config.seglists[other_color];
    std::string spin = (color == 0) ? "up" : "down";
    auto result      = std::tuple<long, long, tau_t, bool>{0, 0, tau_t::zero(), true};

    // --------- Eliminate cases where move is impossible ---------
    if (sl.empty()) {
      LOG("Spin up is empty, cannot regroup spin.");
      return result;
    }
    if (is_full_line(sl[0])) {
      LOG("Spin up is full line, cannot regroup spin.");
      return result;
    }
    if (dsl.empty()) {
      LOG("Spin down is empty, cannot regroup spin.");
      return result;
    }
    if (is_full_line(dsl[0])) {
      LOG("Spin down is full line, cannot regroup spin.");
      return result;
    }

    // --------- Randomly choose a c operator -----------
    auto idx_c = rng(sl.size());
    LOG("Spin {}: regrouping c at position {}.", spin, idx_c);
    if (sl[idx_c].J_c) {
      LOG("Spin {}: cannot regroup because c is connected to a spin line.", spin);
      return result;
    }
    tau_t tau_c = sl[idx_c].tau_c;

    // -------- Propose new position for the c ---------

    // Determine window in which the c can be moved
    auto idx_left       = (idx_c == 0) ? sl.size() - 1 : idx_c - 1;
    tau_t wtau_left     = sl[idx_left].tau_cdag;
    tau_t wtau_right    = sl[idx_c].tau_cdag;
    tau_t window_length = (idx_left == idx_c) ? tau_t::beta() : wtau_left - wtau_right;

    // Find the cdag in opposite spin that are within the window
    auto cdag_list = cdag_in_window(wtau_left, wtau_right, dsl);
    if (cdag_list.empty()) {
      LOG("Spin {}: cannot regroup because there are no suitable cdag operators.", spin);
      return result;
    }
    // Choose one of them randomly
    auto idx_cdag = cdag_list[rng(cdag_list.size())];
    if (dsl[idx_cdag].J_cdag) {
      LOG("Spin {}: cannot regroup because chosen cdag is connected to a spin line.", spin);
      return result;
    }
    tau_t tau_c_new = dsl[idx_cdag].tau_cdag;
    auto new_seg    = segment_t{tau_c_new, sl[idx_c].tau_cdag};
    LOG("Spin {}: moving c from {} to {}.", spin, tau_c, tau_c_new);

    // -------- Trace ratio ---------
    for (auto const &[c, slc] : itertools::enumerate(config.seglists)) {
      if (c != color) {
        ln_trace_ratio += -wdata.U(c, color) * overlap(slc, new_seg);
        ln_trace_ratio -= -wdata.U(c, color) * overlap(slc, sl[idx_c]);
      }
      ln_trace_ratio +=
         K_overlap(slc, tau_c_new, true, wdata.K, c, color) - K_overlap(slc, tau_c, true, wdata.K, c, color);
    }

    // -------- Det ratio ----------
    auto &D          = wdata.dets[color];
    auto det_c_time  = [&](long i) { return D.get_y(i).first; };
    long det_index_c = lower_bound(det_c_time, D.size(), sl[idx_c].tau_c);
    det_ratio *= D.try_change_col(det_index_c, {tau_c_new, 0});

    // --------- Prop ratio ---------
    prop_ratio *= (double(sl.size()) * cdag_list.size()) / window_length;

    std::get<0>(result) = idx_c;
    std::get<1>(result) = idx_cdag;
    std::get<2>(result) = tau_c_new;
    std::get<3>(result) = false;
    return result;
  }

}; // namespace moves
