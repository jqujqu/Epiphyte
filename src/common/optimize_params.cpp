/*****************************************************************************
 *  Copyright (C) 2016 University of Southern California and
 *                     Jenny Qu and Andrew D Smith
 *
 *  Authors: Jenny Qu and Andrew D Smith
 *
 *  This program is free software: you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see
 *  <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "param_set.hpp"
#include "sufficient_statistics_helpers.hpp"
#include "epiphy_utils.hpp"

#include <iostream>
#include <vector>
#include <string>

using std::string;
using std::vector;
using std::endl;
using std::cerr;

using std::max;
using std::min;
using std::abs;
using std::log;

using std::pair;


double
log_likelihood(const vector<size_t> &subtree_sizes, const param_set &ps,
               const pair<double, double> &root_start_counts,
               const pair_state &root_counts,
               const vector<pair_state> &start_counts, // dim=[treesize x 2 x 2]
               const vector<triple_state> &triad_counts) { // dim=[treesize x 2 x 2 x 2]

  vector<pair_state> P;
  vector<triple_state> GP;
  get_transition_matrices(ps, P, GP);

  double llk =
    (root_counts(0, 0)*log(ps.g0) + root_counts(0, 1)*log(1.0 - ps.g0) +
     root_counts(1, 0)*log(1.0 - ps.g1) + root_counts(1, 1)*log(ps.g1));

  for (size_t node = 1; node < subtree_sizes.size(); ++node)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        llk += start_counts[node](j, k)*log(P[node](j, k));
        for (size_t i = 0; i < 2; ++i)
          llk += triad_counts[node](i, j, k)*log(GP[node](i, j, k));
      }
  return llk;
}

static void
objective_branch(const param_set &ps,
                 const vector<pair_state> &start_counts,
                 const vector<triple_state> &triad_counts,
                 const vector<pair_state> &P,
                 const vector<triple_state> &GP,
                 const vector<triple_state> &GP_dT,
                 const size_t node_id,
                 double &F, double &deriv) {

  F = 0.0;
  deriv = 0.0;
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
        deriv += triad_counts[node_id](i, j, k)*
          (GP_dT[node_id](i, j, k)/GP[node_id](i, j, k));
      }

  double rate0 = ps.rate0;
  pair_state P_dT(-rate0, rate0, 1.0 - rate0, rate0 - 1.0);

  for (size_t j = 0; j < 2; ++j)
    for (size_t k = 0; k < 2; ++k) {
      F += start_counts[node_id](j, k)*log(P[node_id](j, k));
      deriv += start_counts[node_id](j, k)*(P_dT(j, k)/P[node_id](j, k));
    }
}


template <class T> static bool
btwn_01_eps(const T x, const double &epsilon) {
  return x > (0.0 + epsilon) && x < (1.0 - epsilon);
}


static double
find_next_branch(const param_set &ps, const double deriv,
                 const size_t node_id, double step_size, param_set &next_ps) {

  while (step_size > param_set::tolerance &&
         !btwn_01_eps(ps.T[node_id] + step_size*sign(deriv), param_set::tolerance))
    step_size /= 2.0;

  next_ps = ps;
  next_ps.T[node_id] = ps.T[node_id] + step_size*sign(deriv);

  return step_size;
}


void
max_likelihood_branch(const bool VERBOSE, const vector<size_t> &subtree_sizes,
                      const size_t node_id,
                      const vector<pair_state> &start_counts,
                      const vector<triple_state> &triad_counts,
                      param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  double F = 0.0, deriv = 0.0;
  objective_branch(ps, start_counts, triad_counts, P, GP, GP_dT, node_id, F, deriv);

  double step_size = 1.0;
  while (step_size > param_set::tolerance) {

    param_set next_ps(ps);
    step_size = find_next_branch(ps, deriv, node_id, step_size, next_ps);

    get_transition_matrices_deriv(next_ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

    double next_F = 0.0, next_deriv = 0.0;
    objective_branch(next_ps, start_counts, triad_counts, P, GP,
                     GP_dT, node_id, next_F, next_deriv);

    if (next_F > F) {
      if (VERBOSE)
        cerr << "[max_likelihood_branch: "
             << "delta=" << next_F - F << ", "
             << "step_size=" << step_size << ", "
             << "branch[" << node_id << "]="
             << next_ps.T[node_id] << ']' << endl;
      F = next_F;
      deriv = next_deriv;
      ps.T[node_id] = next_ps.T[node_id];
    }
    else step_size /= 2.0;
  }
}


static void
objective_rate(const vector<size_t> &subtree_sizes,
               const param_set &ps,
               const vector<pair_state> &start_counts,
               const vector<triple_state> &triad_counts,
               const vector<pair_state> &P,
               const vector<triple_state> &GP,
               const vector<triple_state> &GP_drate,
               double &F, double &deriv_rate) {

  const size_t n_nodes = subtree_sizes.size();

  F = 0.0;
  deriv_rate = 0.0;
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k) {
          F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
          deriv_rate += triad_counts[node_id](i, j, k)*
            (GP_drate[node_id](i, j, k)/GP[node_id](i, j, k));
        }

    const double T_val = ps.T[node_id];
    const pair_state P_drate(-T_val, T_val, -T_val, T_val);
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        F += start_counts[node_id](j, k)*log(P[node_id](j, k));
        deriv_rate += start_counts[node_id](j, k)*
          P_drate(j, k)/P[node_id](j, k);
      }
  }
}


static double
find_next_rate(const param_set &ps, const double deriv,
               double step_size, param_set &next_ps) {

  const double denom = abs(deriv);
  while (step_size > param_set::tolerance &&
         !btwn_01_eps(ps.rate0 + step_size*(deriv/denom), param_set::tolerance))
    step_size /= 2.0;

  next_ps = ps;
  next_ps.rate0 = ps.rate0 + step_size*(deriv/denom);

  return step_size;
}


void
max_likelihood_rate(const bool VERBOSE, const vector<size_t> &subtree_sizes,
                    const vector<pair_state> &start_counts,
                    const vector<triple_state> &triad_counts,
                    param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  double F = 0.0;
  double deriv = 0.0;
  objective_rate(subtree_sizes, ps, start_counts, triad_counts,
                 P, GP, GP_drate, F, deriv);

  double step_size = 1.0;
  while (step_size > param_set::tolerance) {

    param_set next_ps(ps);
    step_size = find_next_rate(ps, deriv, step_size, next_ps);

    get_transition_matrices_deriv(next_ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

    double next_F = 0.0;
    double next_deriv = 0.0;
    objective_rate(subtree_sizes, next_ps, start_counts, triad_counts,
                   P, GP, GP_drate, next_F, next_deriv);

    if (next_F > F) {
      if (VERBOSE)
        cerr << "[max_likelihood_rate: "
             << "delta=" << next_F - F << ", "
             << "step_size=" << step_size << ", "
             << "rate0=" << next_ps.rate0 << ']' << endl;
      F = next_F;
      deriv = next_deriv;
      ps.rate0 = next_ps.rate0;
    }
    else step_size /= 2.0;
  }
}


static void
objective_horiz(const vector<size_t> &subtree_sizes, const param_set &ps,
                const pair_state &root_counts,
                const vector<triple_state> &triad_counts,
                const vector<triple_state> &GP,
                const vector<triple_state> &GP_dg0,
                const vector<triple_state> &GP_dg1,
                double &F, pair<double, double> &deriv_G) { // first=0, second=1

  const size_t n_nodes = subtree_sizes.size();

  F = 0.0;
  deriv_G = std::make_pair(0.0, 0.0);
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

    for (size_t i = 0; i < 2; ++i) // for previous
      for (size_t j = 0; j < 2; ++j) // for parent
        for (size_t k = 0; k < 2; ++k) // for current
          F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));

    for (size_t j = 0; j < 2; ++j) // for parent
      for (size_t k = 0; k < 2; ++k) { // for current
        deriv_G.first += triad_counts[node_id](0, j, k)*
          (GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
        deriv_G.second += triad_counts[node_id](1, j, k)*
          (GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
      }
  }

  const pair_state G(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);

  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      F += root_counts(i, j)*log(G(i, j));

  deriv_G.first += root_counts(0, 0)/G(0, 0) - root_counts(0, 1)/G(0, 1);
  deriv_G.second += -1.0*root_counts(1, 0)/G(1, 0) + root_counts(1, 1)/G(1, 1);
}


static double
find_next_horiz(const param_set &ps, const pair<double, double> &deriv,
                double step_size, param_set &next_ps) {

  const double denom = abs(deriv.first) + abs(deriv.second);
  while (step_size > param_set::tolerance &&
         !(btwn_01_eps(ps.g0 + step_size*(deriv.first/denom),
                       param_set::tolerance) &&
           btwn_01_eps(ps.g1 + step_size*(deriv.second/denom),
                       param_set::tolerance)))
    step_size /= 2.0;

  next_ps = ps;
  next_ps.g0 = ps.g0 + step_size*(deriv.first/denom);
  next_ps.g1 = ps.g1 + step_size*(deriv.second/denom);

  return step_size;
}


void
max_likelihood_horiz(const bool VERBOSE,
                     const vector<size_t> &subtree_sizes,
                     const pair_state &root_counts,
                     const vector<triple_state> &triad_counts,
                     param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  double F = 0.0;
  pair<double, double> deriv(0.0, 0.0);
  objective_horiz(subtree_sizes, ps, root_counts, triad_counts,
                  GP, GP_dg0, GP_dg1, F, deriv);

  double step_size = 1.0;
  while (step_size > param_set::tolerance) {

    param_set next_ps;
    step_size = find_next_horiz(ps, deriv, step_size, next_ps);

    get_transition_matrices_deriv(next_ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

    double next_F = 0.0;
    pair<double, double> next_deriv;
    objective_horiz(subtree_sizes, next_ps, root_counts, triad_counts,
                    GP, GP_dg0, GP_dg1, next_F, next_deriv);

    // update if we have improved, otherwise reduce step size
    if (next_F > F) {
      if (VERBOSE)
        cerr << "[update_G: "
             << "delta=" << next_F - F << ", "
             << "step_size=" << step_size << ", "
             << "G=(" << next_ps.g0 << ", " << next_ps.g1 << ")]" << endl;
      F = next_F;
      deriv.swap(next_deriv);
      ps.g0 = next_ps.g0;
      ps.g1 = next_ps.g1;
    }
    else step_size /= 2.0;
  }
}


void
max_likelihood_pi0(const bool VERBOSE,
                   const pair<double, double> &root_start_counts,
                   param_set &ps) {
  ps.pi0 = root_start_counts.first/(root_start_counts.first +
                                    root_start_counts.second);
  if (VERBOSE)
    cerr << "[max_likelihood_pi0: pi0=" << ps.pi0 << ']' << endl;
}


void
optimize_params(const bool VERBOSE, const vector<size_t> &subtree_sizes,
                const vector<size_t> &reset_points,
                const pair<double, double> &root_start_counts,
                const pair_state &root_counts,
                const vector<pair_state> &start_counts,
                const vector<triple_state> &triad_counts, param_set &ps) {

  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id)
    max_likelihood_branch(VERBOSE, subtree_sizes, node_id,
                          start_counts, triad_counts, ps);

  max_likelihood_rate(VERBOSE, subtree_sizes, start_counts, triad_counts, ps);

  max_likelihood_pi0(VERBOSE, root_start_counts, ps);

  max_likelihood_horiz(VERBOSE, subtree_sizes, root_counts, triad_counts, ps);
}
