#include "remove_segment.hpp"
#include "../logs.hpp"

namespace moves {

  double remove_segment::attempt() {

    LOG("\n =================== ATTEMPT REMOVE ================ \n");

    // ------------ Choice of segment --------------

    // Select removal color
    color    = rng(wdata.n_color);
    auto &sl = config.seglists[color];
    LOG("Removing at color {}", color);

    // If color is empty, nothing to remove
    if (sl.empty()) {
      LOG("Line is empty.");
      return 0;
    }

    // Select segment to remove
    proposed_segment_index = rng(sl.size());
    proposed_segment       = sl[proposed_segment_index];
    if (proposed_segment.tau_c == wdata.qmc_beta and proposed_segment.tau_cdag == wdata.qmc_zero) {
      LOG("Line is full.");
      return 0; // If segment is a full line do not remove
    }
    auto proposed_segment_length = double(proposed_segment.tau_c - proposed_segment.tau_cdag);

    LOG("Removing c at {}, cdag at {}", double(proposed_segment.tau_c), double(proposed_segment.tau_cdag));

    // ------------  Trace ratio  -------------
    double ln_trace_ratio = 0;
    for (auto c : range(wdata.n_color)) {
      if (c != color) {
        ln_trace_ratio -= -wdata.U(color, c) * overlap(config.seglists[c], proposed_segment, fac);
        ln_trace_ratio -= wdata.mu(c) * proposed_segment_length;
        if (wdata.has_Dt)
          ln_trace_ratio -= K_overlap(config.seglists[c], proposed_segment, slice_target_to_scalar(wdata.K, color, c));
      }
    }
    double trace_ratio = std::exp(ln_trace_ratio);

    // ------------  Det ratio  ---------------

    auto det_ratio = wdata.dets[color].try_remove(proposed_segment_index, proposed_segment_index);

    // ------------  Proposition ratio ------------

    double current_number_segments = sl.size();
    double future_number_intervals = std::max(1.0, sl.size() - 1.0);
    // Limits of insertion interval for reverse move, initialise at (beta,0)
    qmc_time_t tau_left  = wdata.qmc_beta;
    qmc_time_t tau_right = wdata.qmc_zero;
    if (current_number_segments != 1) {
      bool is_last_segment  = proposed_segment_index == sl.size() - 1;
      bool is_first_segment = proposed_segment_index == 0;
      // Look at segments on either side of proposed_segment, accounting for cyclicity
      tau_right = sl[is_last_segment ? 0 : proposed_segment_index + 1].tau_c;
      tau_left  = sl[is_first_segment ? sl.size() - 1 : proposed_segment_index - 1].tau_cdag;
    }
    qmc_time_t l = (current_number_segments == 1) ? wdata.qmc_beta : tau_left - tau_right;
    double prop_ratio =
       (future_number_intervals * l * l / (current_number_segments == 1 ? 1 : 2)) / current_number_segments;
    LOG("trace_ratio  = {}, prop_ratio = {}, det_ratio = {}", trace_ratio, prop_ratio, det_ratio);

    det_sign    = (det_ratio > 0) ? 1.0 : -1.0;
    double prod = trace_ratio * det_ratio * prop_ratio;

    return (std::isfinite(prod)) ? prod : det_sign;
  }

  //--------------------------------------------------

  double remove_segment::accept() {

    LOG("\n - - - - - ====> ACCEPT - - - - - - - - - - -\n");

    // Update the dets
    wdata.dets[color].complete_operation();

    auto &sl = config.seglists[color];
    // Compute the sign ratio
    double sign_ratio = is_cyclic(sl[proposed_segment_index]) ? -1 : 1;
    LOG("Sign ratio is {}", sign_ratio);
    ALWAYS_EXPECTS((sign_ratio * det_sign == 1.0),
                   "Error: move has produced negative sign! Det sign is {} and additional sign is {}.", det_sign,
                   sign_ratio);

    // Remove the segment
    sl.erase(sl.begin() + proposed_segment_index);

    // Check invariant
    check_invariant(config);
    SPDLOG_TRACE("Configuration is {}", config);

    return sign_ratio;
  }

  //--------------------------------------------------
  void remove_segment::reject() {
    LOG("\n - - - - - ====> REJECT - - - - - - - - - - -\n");
    wdata.dets[color].reject_last_try();
  }
} // namespace moves
