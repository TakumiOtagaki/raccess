/*
 * GrammarApi: thin wrapper for biological constraints and log-space utilities.
 */
#ifndef RACCESS__GRAMMAR_API_HPP
#define RACCESS__GRAMMAR_API_HPP
#include "raccess/score.hpp"
#include "raccess/prob_model.hpp"
namespace Raccess {
template <typename ProbModelT> class GrammarApi {
public:
  typedef ProbModelT PM;
  typedef typename PM::IntT IntT;
  explicit GrammarApi(PM& pm) : _pm(&pm) {}
  // Biological constraints.
  bool allow_pair(IntT i, IntT j) const { return _pm->allow_pair(i, j); }
  bool allow_inner_loop(IntT i, IntT j) const { return _pm->allow_inner_loop(i, j); }
  bool allow_outer_loop(IntT i, IntT j) const { return _pm->allow_outer_loop(i, j); }
  // Log-space utilities (re-export from score.hpp).
  static ScoreT neg_inf() { return NEG_INF(); }
  static bool impossible(ScoreT sc) { return ::impossible(sc); }
  static void logadd(ScoreT& a, const ScoreT& b) { LOGADD(a, b); }
private:
  PM* _pm;
};
} // namespace Raccess
#endif
