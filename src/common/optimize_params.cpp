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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
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

static const double TOL = 1e-10;

double
log_likelihood(const vector<size_t> &subtree_sizes, const param_set &ps,
               const pair<double, double> &root_start_counts,
               const pair_state &root_counts,
               const vector<pair_state> &start_counts, // dim=[treesize x 2 x 2]
               const vector<triple_state> &triad_counts) { // dim=[treesize x 2 x 2 x 2]
  vector<pair_state> P;
  vector<triple_state> GP;
  get_transition_matrices(ps, P, GP);
  pair_state logG(log(ps.g0[0]), log(1.0 - ps.g0[0]),
                  log(1.0 - ps.g1[0]), log(ps.g1[0]));
  double llk =
    (root_start_counts.first*log(ps.pi0) +
     root_start_counts.second*log(1.0 - ps.pi0)) + 
    (root_counts(0, 0)*logG(0, 0) + root_counts(0, 1)*logG(0, 1) +
     root_counts(1, 0)*logG(1, 0) + root_counts(1, 1)*logG(1, 1));
  for (size_t node = 1; node < subtree_sizes.size(); ++node)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        llk += start_counts[node](j, k)*log(P[node](j, k));
        for (size_t i = 0; i < 2; ++i)
          llk += triad_counts[node](i, j, k)*log(GP[node](i, j, k));
      }
  return llk;
}


// Treat sites as independent
double
log_likelihood(const vector<size_t> &subtree_sizes, const param_set &ps,
               const pair<double, double> &root_start_counts,
               const pair_state &root_counts,
               const vector<pair_state> &start_counts) { // dim=[treesize x 2 x 2]
  vector<pair_state> P;
  vector<triple_state> GP;
  get_transition_matrices(ps, P, GP);
  pair_state logG(log(ps.g0[0]), log(1.0 - ps.g0[0]),
                  log(1.0 - ps.g1[0]), log(ps.g1[0]));
  double llk =
    (root_start_counts.first*log(ps.pi0) +
     root_start_counts.second*log(1.0 - ps.pi0)) + // ADS: is this right?
    (root_counts(0, 0)*logG(0, 0) + root_counts(0, 1)*logG(0, 1) +
     root_counts(1, 0)*logG(1, 0) + root_counts(1, 1)*logG(1, 1));
  for (size_t node = 1; node < subtree_sizes.size(); ++node)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k)
        llk += start_counts[node](j, k)*log(P[node](j, k));
  return llk;
}


// where is this used?
double
log_likelihood(const vector<size_t> &subtree_sizes,
               const pair<double, double> &root_start_counts,
               const pair_state &root_counts,
               const vector<pair_state> &start_counts, // dim=[treesize x 2 x 2]
               const vector<triple_state> &triad_counts) { // dim=[treesize x 2 x 2 x 2]

  const double denom = root_start_counts.first + root_start_counts.second;
  std::pair<double, double> log_pi =
    std::make_pair(log(root_start_counts.first/denom),
                   log(root_start_counts.second/denom));
  pair_state logG(root_counts);
  logG.to_probabilities();
  logG.make_logs();

  std::vector<pair_state> logP;
  std::vector<triple_state> logGP;
  for (size_t i = 0; i < start_counts.size(); ++i) {
    logP.push_back(start_counts[i]);
    logP[i].to_probabilities();
    logP[i].make_logs();
    logGP.push_back(triad_counts[i]);
    logGP[i].to_probabilities();
    logGP[i].make_logs();
  }

  double llk =
    (root_start_counts.first*log_pi.first +
     root_start_counts.second*log_pi.second) + // ADS: is this right?
    (root_counts(0, 0)*logG(0, 0) + root_counts(0, 1)*logG(0, 1) +
     root_counts(1, 0)*logG(1, 0) + root_counts(1, 1)*logG(1, 1));

  for (size_t node = 1; node < subtree_sizes.size(); ++node)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        llk += start_counts[node](j, k)*logP[node](j, k);
        for (size_t i = 0; i < 2; ++i)
          llk += triad_counts[node](i, j, k)*logGP[node](i, j, k);
      }
  return llk;
}

////////////////////////////////////////////////////////////////////////////////
//////// Optimize inividual parameters to maximize log-likelihood  /////////////
////////////////////////////////////////////////////////////////////////////////
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


static bool
btwn_01_eps(const param_set &ps, const double &epsilon,
            const param_set &deriv, double step_size) {
  bool valid = true;
  valid = valid & btwn_01_eps(ps.rate0 + step_size*deriv.rate0, epsilon);
  for (size_t i = 0; i < ps.T.size(); ++i) {
    valid = valid & btwn_01_eps(ps.g0[i] + step_size*deriv.g0[i], epsilon);
    valid = valid & btwn_01_eps(ps.g1[i] + step_size*deriv.g1[i], epsilon);
    if (i > 0)
      valid = valid & btwn_01_eps(ps.T[i] + step_size*deriv.T[i], epsilon);
  }
  return valid;
}


static double
bound_01_eps(const double x, const double epsilon) {
  return min(max(x, epsilon), 1.0-epsilon);
}


static double
find_next_branch(const param_set &ps, const double deriv,
                 const size_t node_id, double step_size, param_set &next_ps) {
  const double sgn = sign(deriv);
  while (step_size > TOL &&
         !btwn_01_eps(ps.T[node_id] + step_size*sgn, TOL))
    step_size /= 2.0;

  next_ps = ps;
  next_ps.T[node_id] = bound_01_eps(ps.T[node_id] + step_size*sgn, TOL);

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
  objective_branch(ps, start_counts, triad_counts,
                   P, GP, GP_dT, node_id, F, deriv);
  double step_size = 1.0;
  while (step_size > TOL) {
    param_set next_ps(ps);
    step_size = find_next_branch(ps, deriv, node_id, step_size, next_ps);
    get_transition_matrices_deriv(next_ps, P, GP, GP_drate,
                                  GP_dg0, GP_dg1, GP_dT);
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
    } else {
      step_size /= 2.0;
    }
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
  const double sgn = sign(deriv);
  while (step_size > TOL && !btwn_01_eps(ps.rate0 + step_size*sgn, TOL))
    step_size /= 2.0;
  next_ps = ps;
  next_ps.rate0 = bound_01_eps(ps.rate0 + step_size*sgn, TOL);
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
  while (step_size > TOL) {
    param_set next_ps(ps);
    step_size = find_next_rate(ps, deriv, step_size, next_ps);
    get_transition_matrices_deriv(next_ps, P, GP, GP_drate,
                                  GP_dg0, GP_dg1, GP_dT);
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
    } else {
      step_size /= 2.0;
    }
  }
}


/* horizontal rates: same for all nodes*/
static void
objective_horiz_one_rate(const vector<size_t> &subtree_sizes, const param_set &ps,
                         const pair_state &root_counts,
                         const vector<triple_state> &triad_counts,
                         const vector<triple_state> &GP,
                         const vector<triple_state> &GP_dg0,
                         const vector<triple_state> &GP_dg1,
                         double &F, vector<pair<double, double> > &deriv_G) { // first=0, second=1
  const size_t n_nodes = subtree_sizes.size();
  F = 0.0;
  deriv_G = vector<pair<double, double> >(n_nodes, std::make_pair(0.0, 0.0));
  pair<double, double> dG = std::make_pair(0.0, 0.0);
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    if (node_id == 0) { // root
      pair_state G(ps.g0[0], 1.0 - ps.g0[0], 1.0 - ps.g1[0], ps.g1[0]);
      for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 2; ++j)
          F += root_counts(i, j)*log(G(i, j));
      dG.first += root_counts(0, 0)/G(0, 0) - root_counts(0, 1)/G(0, 1);
      dG.second += -1.0*root_counts(1, 0)/G(1, 0) + root_counts(1, 1)/G(1, 1);
    } else { // non-root
      for (size_t i = 0; i < 2; ++i) // for previous
        for (size_t j = 0; j < 2; ++j) // for parent
          for (size_t k = 0; k < 2; ++k) // for current
            F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
      for (size_t j = 0; j < 2; ++j) // for parent
        for (size_t k = 0; k < 2; ++k) { // for current
          dG.first += triad_counts[node_id](0, j, k)*
            (GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
          dG.second += triad_counts[node_id](1, j, k)*
            (GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
        }
    }
  }
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    // non-root nodes share horizontal rates and derivatives
    deriv_G[node_id].first = dG.first;
    deriv_G[node_id].second = dG.second;
  }
}


/* horizontal rates: root vs non-root*/
static void
objective_horiz_root_nonroot(const vector<size_t> &subtree_sizes,
                             const param_set &ps,
                             const pair_state &root_counts,
                             const vector<triple_state> &triad_counts,
                             const vector<triple_state> &GP,
                             const vector<triple_state> &GP_dg0,
                             const vector<triple_state> &GP_dg1,
                             double &F,
                             vector<pair<double, double> > &deriv_G) { // first=0, second=1
  const size_t n_nodes = subtree_sizes.size();
  F = 0.0;
  deriv_G = vector<pair<double, double> >(n_nodes, std::make_pair(0.0, 0.0));
  pair<double, double> dG = std::make_pair(0.0, 0.0);
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    if (node_id == 0) { // root
      pair_state G(ps.g0[0], 1.0 - ps.g0[0], 1.0 - ps.g1[0], ps.g1[0]);
      for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 2; ++j)
          F += root_counts(i, j)*log(G(i, j));
      deriv_G[0].first += root_counts(0, 0)/G(0, 0) - root_counts(0, 1)/G(0, 1);
      deriv_G[0].second += (-1.0*root_counts(1, 0)/G(1, 0) +
                            root_counts(1, 1)/G(1, 1));
    } else { // non-root
      for (size_t i = 0; i < 2; ++i) // for previous
        for (size_t j = 0; j < 2; ++j) // for parent
          for (size_t k = 0; k < 2; ++k) // for current
            F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
      for (size_t j = 0; j < 2; ++j) // for parent
        for (size_t k = 0; k < 2; ++k) { // for current
          dG.first += triad_counts[node_id](0, j, k)*
            (GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
          dG.second += triad_counts[node_id](1, j, k)*
            (GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
        }
    }
  }
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    // non-root nodes share horizontal rates and derivatives
    deriv_G[node_id].first = dG.first;
    deriv_G[node_id].second = dG.second;
  }
}

/* horizontal rates: by node */
static void
objective_horiz_by_node(const vector<size_t> &subtree_sizes,
                        const param_set &ps,
                        const pair_state &root_counts,
                        const vector<triple_state> &triad_counts,
                        const vector<triple_state> &GP,
                        const vector<triple_state> &GP_dg0,
                        const vector<triple_state> &GP_dg1,
                        double &F, vector<pair<double, double> > &deriv_G) { // first=0, second=1
  const size_t n_nodes = subtree_sizes.size();
  F = 0.0;
  deriv_G = vector<pair<double, double> >(n_nodes, std::make_pair(0.0, 0.0));
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    if (node_id == 0) { // root
      pair_state G(ps.g0[0], 1.0 - ps.g0[0], 1.0 - ps.g1[0], ps.g1[0]);
      for (size_t i = 0; i < 2; ++i)
        for (size_t j = 0; j < 2; ++j)
          F += root_counts(i, j)*log(G(i, j));
      deriv_G[0].first += root_counts(0, 0)/G(0, 0) - root_counts(0, 1)/G(0, 1);
      deriv_G[0].second += (-1.0*root_counts(1, 0)/G(1, 0) +
                            root_counts(1, 1)/G(1, 1));
    } else { // non-root
      for (size_t i = 0; i < 2; ++i) // for previous
        for (size_t j = 0; j < 2; ++j) // for parent
          for (size_t k = 0; k < 2; ++k) // for current
            F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
      for (size_t j = 0; j < 2; ++j) // for parent
        for (size_t k = 0; k < 2; ++k) { // for current
          deriv_G[node_id].first += triad_counts[node_id](0, j, k)*
            (GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
          deriv_G[node_id].second += triad_counts[node_id](1, j, k)*
            (GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
        }
    }
  }
}


static double
find_next_horiz(const param_set &ps, const vector<pair<double, double> > &deriv,
                double step_size, param_set &next_ps) {
  double denom = 0.0;
  for (size_t i = 0; i < deriv.size(); ++i)
    denom = max(denom, abs(deriv[i].first) + abs(deriv[i].second));
  bool bound = false;
  step_size *= 2.0;
  while (step_size > TOL && !bound) {
    step_size /= 2.0;
    bound = true;
    for (size_t i = 0; i < deriv.size(); ++i)
      bound = bound & ( btwn_01_eps(ps.g0[i] +
                                    step_size*(deriv[i].first/denom), TOL) &&
                        btwn_01_eps(ps.g1[i]+
                                    step_size*(deriv[i].second/denom), TOL));
  }
  next_ps = ps;
  for (size_t i = 0; i < deriv.size(); ++i) {
    next_ps.g0[i] = bound_01_eps(ps.g0[i] +
                                 step_size*(deriv[i].first/denom), TOL);
    next_ps.g1[i] = bound_01_eps(ps.g1[i] +
                                 step_size*(deriv[i].second/denom), TOL);
  }
  return step_size;
}


void
max_likelihood_horiz(const bool VERBOSE,
                     const size_t HORIZ_MODE,                     
                     const vector<size_t> &subtree_sizes,
                     const pair_state &root_counts,
                     const vector<triple_state> &triad_counts,
                     param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  double F = 0.0;
  vector<pair<double, double> > deriv;
  switch (HORIZ_MODE) {
  case 1:
    objective_horiz_one_rate(subtree_sizes, ps, root_counts, triad_counts,
                             GP, GP_dg0, GP_dg1, F, deriv);
    assert (deriv[0].first == deriv[1].first);
    break;
  case 2:
    objective_horiz_root_nonroot(subtree_sizes, ps, root_counts, triad_counts,
                                 GP, GP_dg0, GP_dg1, F, deriv);
    break;
  case 3:
    objective_horiz_by_node(subtree_sizes, ps, root_counts, triad_counts,
                            GP, GP_dg0, GP_dg1, F, deriv);
    break;
  }
  double step_size = 1.0;
  while (step_size > TOL) {
    param_set next_ps;
    step_size = find_next_horiz(ps, deriv, step_size, next_ps);
    get_transition_matrices_deriv(next_ps, P, GP, GP_drate,
                                  GP_dg0, GP_dg1, GP_dT);
    double next_F = 0.0;
    vector<pair<double, double> > next_deriv;
    switch (HORIZ_MODE) {
    case 1:
      objective_horiz_one_rate(subtree_sizes, next_ps, root_counts,
                               triad_counts, GP, GP_dg0, GP_dg1, next_F,
                               next_deriv);
      assert (deriv[0].first == deriv[1].first);
      break;
    case 2:
      objective_horiz_root_nonroot(subtree_sizes, next_ps, root_counts,
                                   triad_counts, GP, GP_dg0, GP_dg1, next_F,
                                   next_deriv);
      break;
    case 3:
      objective_horiz_by_node(subtree_sizes, next_ps, root_counts, triad_counts,
                              GP, GP_dg0, GP_dg1, next_F, next_deriv);
      break;
    }
    // update if we have improved, otherwise reduce step size
    if (next_F > F) {
      if (VERBOSE) {
        cerr << "[update_G: "
             << "delta=" << next_F - F << ", "
             << "step_size=" << step_size << ", "
             << "G={";
        for (size_t i = 0; i < subtree_sizes.size();++i)
          cerr << "("<< next_ps.g0[i] << ", " << next_ps.g1[i] << ")";
        cerr << "}" << endl;
      }
      F = next_F;
      deriv.swap(next_deriv);
      ps.g0 = next_ps.g0;
      ps.g1 = next_ps.g1;
    } else {
      step_size /= 2.0;
    }
  }
}


void
max_likelihood_pi0(const bool VERBOSE,
                   const pair<double, double> &root_start_counts,
                   param_set &ps) {
  ps.pi0 = root_start_counts.first/(root_start_counts.first +
                                    root_start_counts.second);

  if (ps.pi0 == 0.0 || ps.pi0 == 1.0) ps.pi0 = 0.5;
  // ADS: needs to be improved to use more information
  if (VERBOSE)
    cerr << "[max_likelihood_pi0: pi0=" << ps.pi0 << ']' << endl;
}


// may get rid of subtree_sizes
// Optimize branch and rate parameters separately
void
optimize_params(const bool VERBOSE, const size_t HORIZ_MODE,
                const vector<size_t> &subtree_sizes,
                const pair<double, double> &root_start_counts,
                const pair_state &root_counts,
                const vector<pair_state> &start_counts,
                const vector<triple_state> &triad_counts, param_set &ps) {
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id)
    max_likelihood_branch(VERBOSE, subtree_sizes, node_id,
                          start_counts, triad_counts, ps);
  max_likelihood_rate(VERBOSE, subtree_sizes, start_counts, triad_counts, ps);
  max_likelihood_pi0(VERBOSE, root_start_counts, ps);
  max_likelihood_horiz(VERBOSE, HORIZ_MODE, subtree_sizes, root_counts,
                       triad_counts, ps);
}


////////////////////////////////////////////////////////////////////////////////
//////// Optimize by gradient ascent to maximize log-likelihood   //////////////
////////////////////////////////////////////////////////////////////////////////
static void
objective_params_one_rate(const param_set &ps,
                          const pair<double, double> &root_start_counts,
                          const pair_state &root_counts,
                          const vector<pair_state> &start_counts,
                          const vector<triple_state> &triad_counts,
                          const pair_state &root_G,
                          const vector<pair_state> &P,
                          const vector<triple_state> &GP,
                          const vector<pair_state> &P_drate,
                          const pair_state &P_dT, //share by all T
                          const vector<triple_state> &GP_drate,
                          const vector<triple_state> &GP_dg0,
                          const vector<triple_state> &GP_dg1,
                          const vector<triple_state> &GP_dT,
                          double &F, param_set &deriv) {
  //set all values to 0
  F = 0.0;
  deriv = ps;
  std::fill(deriv.T.begin(), deriv.T.end(), 0.0);
  deriv.pi0 = 0.0;
  deriv.rate0 = 0.0;
  deriv.g0 = vector<double>(ps.T.size(), 0.0);
  deriv.g1 = vector<double>(ps.T.size(), 0.0);

  // root_start_counts contribute to pi0
  deriv.pi0 = (root_start_counts.first/ps.pi0 -
               root_start_counts.second/(1.0 - ps.pi0));
  // root_counts countribute to g0, g1
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j) 
      F += root_counts(i, j)*log(root_G(i, j));

  double dg0 = 0.0, dg1 = 0.0;
  dg0 += root_counts(0, 0)/root_G(0, 0) - root_counts(0, 1)/root_G(0, 1);
  dg1 += (-1.0*root_counts(1, 0)/root_G(1, 0) +
                  root_counts(1, 1)/root_G(1, 1));

  const size_t n_nodes = ps.T.size();
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    // triad_counts contribute to rate, g0, g1, and T
    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k) {
          F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
          deriv.rate0 += (triad_counts[node_id](i, j, k)*
                          GP_drate[node_id](i, j, k)/GP[node_id](i, j, k));
          deriv.T[node_id] += (triad_counts[node_id](i, j, k)*
                               GP_dT[node_id](i, j, k)/GP[node_id](i, j, k));
        }
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        dg0 += (triad_counts[node_id](0, j, k)*
                     GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
        dg1 += (triad_counts[node_id](1, j, k)*
                     GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
      }
    // start_counts contribute to rate and T
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        F += start_counts[node_id](j, k)*log(P[node_id](j, k));
        deriv.rate0 += (start_counts[node_id](j, k)*
                        P_drate[node_id](j, k)/P[node_id](j, k));
        deriv.T[node_id] += (start_counts[node_id](j, k)*
                             (P_dT(j, k)/P[node_id](j, k)));
      }
  }
  // copy from the dg0, dg1
  for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
    deriv.g0[node_id] = dg0; 
    deriv.g1[node_id] = dg1;
  }
}


static void
objective_params_root_nonroot(const param_set &ps,
                              const pair<double, double> &root_start_counts,
                              const pair_state &root_counts,
                              const vector<pair_state> &start_counts,
                              const vector<triple_state> &triad_counts,
                              const pair_state &root_G,
                              const vector<pair_state> &P,
                              const vector<triple_state> &GP,
                              const vector<pair_state> &P_drate,
                              const pair_state &P_dT, //share by all T
                              const vector<triple_state> &GP_drate,
                              const vector<triple_state> &GP_dg0,
                              const vector<triple_state> &GP_dg1,
                              const vector<triple_state> &GP_dT,
                              double &F, param_set &deriv) {
  //set all values to 0
  F = 0.0;
  deriv = ps;
  std::fill(deriv.T.begin(), deriv.T.end(), 0.0);
  deriv.pi0 = 0.0;
  deriv.rate0 = 0.0;
  deriv.g0 = vector<double>(ps.T.size(), 0.0);
  deriv.g1 = vector<double>(ps.T.size(), 0.0);

  // root_start_counts contribute to pi0
  deriv.pi0 = (root_start_counts.first/ps.pi0 -
               root_start_counts.second/(1.0 - ps.pi0));

  // root_counts countribute to f0, f1
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j) 
      F += root_counts(i, j)*log(root_G(i, j));

  deriv.g0[0] += (root_counts(0, 0)/root_G(0, 0) -
                  root_counts(0, 1)/root_G(0, 1));
  deriv.g1[0] += (-1.0*root_counts(1, 0)/root_G(1, 0) +
                  root_counts(1, 1)/root_G(1, 1));
  
  double dg0 = 0, dg1 = 0;
  const size_t n_nodes = ps.T.size();
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    // triad_counts contribute to rate, g0, g1, and T
    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k) {
          F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
          deriv.rate0 += (triad_counts[node_id](i, j, k)*
                          GP_drate[node_id](i, j, k)/GP[node_id](i, j, k));
          deriv.T[node_id] += (triad_counts[node_id](i, j, k)*
                               GP_dT[node_id](i, j, k)/GP[node_id](i, j, k));
        }
    // g0 g1 at the first nonroot node
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        dg0 += (triad_counts[node_id](0, j, k)*
                     GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
        dg1 += (triad_counts[node_id](1, j, k)*
                     GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
      }
    // start_counts contribute to rate and T
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        F += start_counts[node_id](j, k)*log(P[node_id](j, k));
        deriv.rate0 += (start_counts[node_id](j, k)*
                        P_drate[node_id](j, k)/P[node_id](j, k));
        deriv.T[node_id] += (start_counts[node_id](j, k)*
                             (P_dT(j, k)/P[node_id](j, k)));
      }
  }
  // nonroot nodes copy from the dg0, dg1
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    deriv.g0[node_id] = dg0; 
    deriv.g1[node_id] = dg1;
  }
}


static void
objective_params_by_node(const param_set &ps,
                         const pair<double, double> &root_start_counts,
                         const pair_state &root_counts,
                         const vector<pair_state> &start_counts,
                         const vector<triple_state> &triad_counts,
                         const pair_state &root_G,
                         const vector<pair_state> &P,
                         const vector<triple_state> &GP,
                         const vector<pair_state> &P_drate,
                         const pair_state &P_dT, //share by all T
                         const vector<triple_state> &GP_drate,
                         const vector<triple_state> &GP_dg0,
                         const vector<triple_state> &GP_dg1,
                         const vector<triple_state> &GP_dT,
                         double &F, param_set &deriv) {
  //set all values to 0
  F = 0.0;
  deriv = ps;
  std::fill(deriv.T.begin(), deriv.T.end(), 0.0);
  deriv.pi0 = 0.0;
  deriv.rate0 = 0.0;
  deriv.g0 = vector<double>(ps.T.size(), 0.0);
  deriv.g1 = vector<double>(ps.T.size(), 0.0);
  // root_start_counts contribute to pi0
  deriv.pi0 = (root_start_counts.first/ps.pi0 -
               root_start_counts.second/(1.0 - ps.pi0));
  // root_counts countribute to g0[0], g1[0]
  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j) 
      F += root_counts(i, j)*log(root_G(i, j));
  
  deriv.g0[0] += (root_counts(0, 0)/root_G(0, 0) -
                  root_counts(0, 1)/root_G(0, 1));
  deriv.g1[0] += (-1.0*root_counts(1, 0)/root_G(1, 0) +
                  root_counts(1, 1)/root_G(1, 1));

  const size_t n_nodes = ps.T.size();
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    // triad_counts contribute to rate, g0, g1, and T
    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k) {
          F += triad_counts[node_id](i, j, k)*log(GP[node_id](i, j, k));
          deriv.rate0 += (triad_counts[node_id](i, j, k)*
                          GP_drate[node_id](i, j, k)/GP[node_id](i, j, k));
          deriv.T[node_id] += (triad_counts[node_id](i, j, k)*
                               GP_dT[node_id](i, j, k)/GP[node_id](i, j, k));
        }   
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        deriv.g0[node_id] += (triad_counts[node_id](0, j, k)*
                              GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
        deriv.g1[node_id] += (triad_counts[node_id](1, j, k)*
                     GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
      }
    // start_counts contribute to rate and T
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        F += start_counts[node_id](j, k)*log(P[node_id](j, k));
        deriv.rate0 += (start_counts[node_id](j, k)*
                        P_drate[node_id](j, k)/P[node_id](j, k));
        deriv.T[node_id] += (start_counts[node_id](j, k)*
                             (P_dT(j, k)/P[node_id](j, k)));
      }
  }
}


static double
find_next_ps(const param_set &ps,
             const param_set &deriv,
             double &step_size, param_set &next_ps) {
  const size_t n_nodes = ps.T.size();
  // find the max gradient
  double denom = max(abs(deriv.pi0), abs(deriv.rate0));
  for (size_t i = 0; i < n_nodes; ++i)
    denom = max(max(abs(deriv.g0[i]), abs(deriv.g1[i])), denom);
  for (size_t i = 1; i < n_nodes; ++i)
    denom = max(denom, abs(deriv.T[i]));

  while (step_size > TOL && !btwn_01_eps(ps, TOL, deriv, step_size))
    step_size /= 2.0;

  next_ps = ps;
  next_ps.pi0 = bound_01_eps(ps.pi0 + step_size*(deriv.pi0/denom), TOL);
  next_ps.rate0 = bound_01_eps(ps.rate0 + step_size*(deriv.rate0/denom), TOL);
  for (size_t i = 0; i < n_nodes; ++i) {
    next_ps.g0[i] = bound_01_eps(ps.g0[i] + step_size*(deriv.g0[i]/denom), TOL);
    next_ps.g1[i] = bound_01_eps(ps.g1[i] + step_size*(deriv.g1[i]/denom), TOL);
  }
  for (size_t i = 1; i < n_nodes; ++i)
    next_ps.T[i] = bound_01_eps(ps.T[i] + step_size*(deriv.T[i]/denom), TOL);

  return step_size;
}


void
max_likelihood_params(const bool VERBOSE,
                      const size_t HORIZ_MODE,
                      const vector<size_t> &subtree_sizes,
                      const pair<double, double> &root_start_counts,
                      const pair_state &root_counts,
                      const vector<pair_state> &start_counts,
                      const vector<triple_state> &triad_counts,
                      param_set &ps) {
  pair_state root_G(ps.g0[0], 1.0 - ps.g0[0], 1.0 - ps.g1[0], ps.g1[0]);
  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  const size_t n_nodes = ps.T.size();
  vector<pair_state>P_drate(n_nodes);
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
    const double T_val = ps.T[node_id];
    P_drate[node_id] = pair_state(-T_val, T_val, -T_val, T_val);
  }
  const double rate0 = ps.rate0;
  pair_state P_dT(-rate0, rate0, 1.0 - rate0, rate0 - 1.0);
  double F; //log-likelihood
  param_set deriv;
  switch (HORIZ_MODE) {
  case 1:
    objective_params_one_rate(ps, root_start_counts,root_counts, start_counts,
                              triad_counts, root_G, P, GP, P_drate, P_dT,
                              GP_drate, GP_dg0, GP_dg1, GP_dT, F, deriv);
    assert (deriv.g0[0] == deriv.g0[1]);
    break;
  case 2:
    objective_params_root_nonroot(ps, root_start_counts,root_counts,
                                  start_counts, triad_counts, root_G, P, GP,
                                  P_drate, P_dT, GP_drate, GP_dg0, GP_dg1,
                                  GP_dT, F, deriv);
  case 3:
    objective_params_by_node(ps, root_start_counts,root_counts, start_counts,
                             triad_counts, root_G, P, GP, P_drate, P_dT,
                             GP_drate, GP_dg0, GP_dg1, GP_dT, F, deriv);
    break;
  } 
  double step_size = 1.0;
  while (step_size > TOL) {
    // find next point
    param_set next_ps(ps);
    step_size = find_next_ps(ps, deriv, step_size, next_ps);
    // evaluate auxiliary quantities at next_ps
    root_G = pair_state(next_ps.g0[0], 1.0 - next_ps.g0[0],
                        1.0 - next_ps.g1[0], next_ps.g1[0]);
    get_transition_matrices_deriv(next_ps, P, GP, GP_drate,
                                  GP_dg0, GP_dg1, GP_dT);
    for (size_t node_id = 1; node_id < n_nodes; ++node_id) {
      const double T_val = next_ps.T[node_id];
      P_drate[node_id] = pair_state(-T_val, T_val, -T_val, T_val);
    }
    const double rate0 = ps.rate0;
    P_dT = pair_state(-rate0, rate0, 1.0 - rate0, rate0 - 1.0);
    // get log-likelihood and gradients at next_ps
    double next_F = 0.0;
    param_set next_deriv;
    switch (HORIZ_MODE) {
    case 1:
      objective_params_one_rate(next_ps, root_start_counts, root_counts,
                                    start_counts, triad_counts, root_G, P, GP,
                                    P_drate, P_dT, GP_drate, GP_dg0, GP_dg1,
                                    GP_dT, next_F, next_deriv);
      assert (deriv.g0[0] == deriv.g0[1]);
      break;
    case 2:
      objective_params_root_nonroot(next_ps, root_start_counts, root_counts,
                                    start_counts, triad_counts, root_G, P, GP,
                                    P_drate, P_dT, GP_drate, GP_dg0, GP_dg1,
                                    GP_dT, next_F, next_deriv);
      break;
    case 3:
      objective_params_by_node(next_ps, root_start_counts, root_counts,
                               start_counts, triad_counts, root_G, P, GP,
                               P_drate, P_dT, GP_drate, GP_dg0, GP_dg1,
                               GP_dT, next_F, next_deriv);
      break;
    }
    if (next_F > F) {
      if (VERBOSE)
        cerr << "[Maximization step]\tlog_lik=" << next_F
             << "\tparam:" << next_ps << endl;
      F = next_F;
      deriv = next_deriv;
      ps = next_ps;
    } else {
      step_size /= 2.0;
    }
  }
}


// may get rid of  subtree_sizes
void
optimize_all_params(const bool VERBOSE,
                    const size_t HORIZ_MODE,                      
                    const vector<size_t> &subtree_sizes,
                    const pair<double, double> &root_start_counts,
                    const pair_state &root_counts,
                    const vector<pair_state> &start_counts,
                    const vector<triple_state> &triad_counts, param_set &ps) {
  max_likelihood_params(VERBOSE, HORIZ_MODE, subtree_sizes, root_start_counts, root_counts,
                        start_counts, triad_counts, ps);
}
