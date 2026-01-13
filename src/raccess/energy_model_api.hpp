/*
 * Minimal adapter around ScoreModelEnergy for beam-search DP.
 */
#ifndef RACCESS__ENERGY_MODEL_API_HPP
#define RACCESS__ENERGY_MODEL_API_HPP
#include "raccess/score_model_energy.hpp"
namespace Raccess {
class EnergyModelApi {
public:
  typedef ScoreModelEnergy SM;
  typedef SM::IntT IntT;
  typedef SM::Seq Seq;
  explicit EnergyModelApi(SM& sm) : _sm(&sm) {}
  void initialize() { _sm->initialize(); }
  // NOTE: ScoreModelEnergy pads sequence at both ends (1-based indexing).
  // Indices i,j are interpreted in Raccess coordinates and are not raw 0-based.
  void set_seq(const Seq& seq) { _sm->set_seq(seq); }
  IntT seqlen() const { return _sm->seqlen(); }
  IntT max_loop() const { return SM::MAXLOOP; }
  IntT min_hairpin() const { return SM::MINHPIN; }
  ScoreT rt_kcal_mol() const { return SM::RT_KCAL_MOL(); }
  ScoreT energy_to_score(ScoreT energy) const { return _sm->energy_to_score(energy); }
  ScoreT score_to_energy(ScoreT score) const { return _sm->score_to_energy(score); }
  // Log Boltzmann factors.
  // Stack: in DP coords this is the stack between pairs (i, j) and (i+1, j-1).
  // In padded seq coordinates it corresponds to (i+1, j) and (i+2, j-1).
  ScoreT log_boltz_stack(IntT i, IntT j) const { return _sm->score_stack(i, j); }
  ScoreT log_boltz_stem_close(IntT i, IntT j) const { return _sm->score_stem_close(i, j); }

  // Hairpin: i and j are inner loop bounds in DP coordinates.
  // Closing pair is (i-1, j+1) in DP, which maps to seq(i) and seq(j+1) after padding.
  ScoreT log_boltz_hairpin(IntT i, IntT j) const { return _sm->score_hairpin(i, j); }
  
  // Interior/bulge: (i-1, j+1) closes the outer pair, (ip, jp) closes the inner pair in DP.
  // Example (1-origin, DP coords): outer (1,10), inner (2,8) bulge length 1:
  // log_boltz_interior(2, 9, 2, 8).
  ScoreT log_boltz_interior(IntT i, IntT j, IntT ip, IntT jp) const {
    return _sm->score_interior(i, j, ip, jp);
  }
  ScoreT log_boltz_loop(IntT i, IntT j, IntT p, IntT q) const {
    // LinearCapR energy_loop equivalent (stack/bulge/internal).
    // stack: no unpaired bases between (i,j) and (p,q)
    if ((p == (i + 1)) && (q == (j - 1))) {
      return _sm->score_stack(i, j);
    }
    return _sm->score_interior(i, j, p, q);
  }
  ScoreT log_boltz_multi_close(IntT i, IntT j) const { return _sm->score_multi_close(i, j); }
  ScoreT log_boltz_multi_open(IntT i, IntT j) const { return _sm->score_multi_open(i, j); }
  ScoreT log_boltz_multi_extend(IntT i, IntT j) const { return _sm->score_multi_extend(i, j); }
  ScoreT log_boltz_outer_extend(IntT i, IntT j) const { return _sm->score_outer_extend(i, j); }
  ScoreT log_boltz_outer_branch(IntT i, IntT j) const { return _sm->score_outer_branch(i, j); }
  // Boltzmann factors.
  ScoreT boltz_stack(IntT i, IntT j) const { return EXP(log_boltz_stack(i, j)); }
  ScoreT boltz_stem_close(IntT i, IntT j) const { return EXP(log_boltz_stem_close(i, j)); }
  ScoreT boltz_hairpin(IntT i, IntT j) const { return EXP(log_boltz_hairpin(i, j)); }
  ScoreT boltz_interior(IntT i, IntT j, IntT ip, IntT jp) const {
    return EXP(log_boltz_interior(i, j, ip, jp));
  }
  ScoreT boltz_loop(IntT i, IntT j, IntT p, IntT q) const {
    return EXP(log_boltz_loop(i, j, p, q));
  }
  ScoreT boltz_multi_close(IntT i, IntT j) const { return EXP(log_boltz_multi_close(i, j)); }
  ScoreT boltz_multi_open(IntT i, IntT j) const { return EXP(log_boltz_multi_open(i, j)); }
  ScoreT boltz_multi_extend(IntT i, IntT j) const { return EXP(log_boltz_multi_extend(i, j)); }
  ScoreT boltz_outer_extend(IntT i, IntT j) const { return EXP(log_boltz_outer_extend(i, j)); }
  ScoreT boltz_outer_branch(IntT i, IntT j) const { return EXP(log_boltz_outer_branch(i, j)); }

  // Closed-pair wrappers (accept closing pairs in "1-origin" "closed coordinates").
  // These avoid manual padding/half-open adjustments in callers.

  // Hairpin: (a, b) closes the hairpin
  ScoreT log_boltz_hairpin_closed(IntT a, IntT b) const {
    return log_boltz_hairpin(a + 1, b - 1);
  }

  // 
  ScoreT boltz_hairpin_closed(IntT a, IntT b) const {
    return EXP(log_boltz_hairpin_closed(a, b));
  }
  ScoreT log_boltz_stack_closed(IntT a, IntT b) const {
    return log_boltz_stack(a - 1, b);
  }
  ScoreT boltz_stack_closed(IntT a, IntT b) const {
    return EXP(log_boltz_stack_closed(a, b));
  }
  ScoreT log_boltz_interior_closed(IntT a, IntT b, IntT c, IntT d) const {
    // outer pair (a,b), inner pair (c,d) with a < c < d < b
    // Maps to dp(i,j,ip,jp)=(a, b-1, c-1, d) to align with score_interior_nuc().
    return log_boltz_interior(a, b - 1, c - 1, d);
  }
  ScoreT boltz_interior_closed(IntT a, IntT b, IntT c, IntT d) const {
    return EXP(log_boltz_interior_closed(a, b, c, d));
  }
  ScoreT log_boltz_loop_closed(IntT a, IntT b, IntT c, IntT d) const {
    // LinearCapR energy_loop equivalent in closed coordinates.
    if ((c == (a + 1)) && (d == (b - 1))) {
      return log_boltz_stack_closed(a, b);
    }
    return log_boltz_interior_closed(a, b, c, d);
  }
  ScoreT boltz_loop_closed(IntT a, IntT b, IntT c, IntT d) const {
    return EXP(log_boltz_loop_closed(a, b, c, d));
  }

  // Multiloop/external wrappers (accept closing pairs in "1-origin" closed coordinates).
  // - multi_close uses inner bounds, like hairpin/interior.
  // - multi_open/outer_branch use paired coordinates (pair = (i+1, j) in padded seq).
  ScoreT log_boltz_multi_close_closed(IntT a, IntT b) const {
    return log_boltz_multi_close(a + 1, b - 1);
  }
  ScoreT boltz_multi_close_closed(IntT a, IntT b) const {
    return EXP(log_boltz_multi_close_closed(a, b));
  }
  ScoreT log_boltz_multi_open_closed(IntT a, IntT b) const {
    return log_boltz_multi_open(a - 1, b);
  }
  ScoreT boltz_multi_open_closed(IntT a, IntT b) const {
    return EXP(log_boltz_multi_open_closed(a, b));
  }
  ScoreT log_boltz_outer_branch_closed(IntT a, IntT b) const {
    return log_boltz_outer_branch(a - 1, b);
  }
  ScoreT boltz_outer_branch_closed(IntT a, IntT b) const {
    return EXP(log_boltz_outer_branch_closed(a, b));
  }
private:
  SM* _sm;
};
} // namespace Raccess
#endif
