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
  void set_seq(const Seq& seq) { _sm->set_seq(seq); }
  IntT seqlen() const { return _sm->seqlen(); }
  IntT max_loop() const { return SM::MAXLOOP; }
  IntT min_hairpin() const { return SM::MINHPIN; }
  ScoreT rt_kcal_mol() const { return SM::RT_KCAL_MOL(); }
  ScoreT energy_to_score(ScoreT energy) const { return _sm->energy_to_score(energy); }
  ScoreT score_to_energy(ScoreT score) const { return _sm->score_to_energy(score); }
  // Log Boltzmann factors.
  ScoreT log_boltz_stack(IntT i, IntT j) const { return _sm->score_stack(i, j); }
  ScoreT log_boltz_stem_close(IntT i, IntT j) const { return _sm->score_stem_close(i, j); }
  ScoreT log_boltz_hairpin(IntT i, IntT j) const { return _sm->score_hairpin(i, j); }
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
private:
  SM* _sm;
};
} // namespace Raccess
#endif
