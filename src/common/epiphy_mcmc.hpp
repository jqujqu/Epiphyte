/* Copyright (C) 2015-16 University of Southern California and
 *                       Andrew D. Smith and Jenny Qu
 *
 * Authors: Jenny Qu and Andrew Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#ifndef EPIPHY_MCMC_HPP
#define EPIPHY_MCMC_HPP

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////   MCMC sampling using Markov Blancket    ////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <epiphy_utils.hpp>
#include <param_set.hpp>
#include <sufficient_statistics_helpers.hpp>

#include <vector>
#include <random>
#include <iostream>

class epiphy_mcmc {
public:
  epiphy_mcmc(const size_t mi, const size_t s) :
    max_iterations(mi), the_seed(s) {
    dis = std::uniform_real_distribution<>(0, 1);
    gen = std::mt19937_64(the_seed);
  }

  //JQU: this functions is not defined/used anywhere yet.
  template <class T> void
  estimate_posterior(const bool VERBOSE,
                     const std::vector<size_t> &subtree_sizes,
                     const std::vector<size_t> &parent_ids,
                     const std::vector<std::vector<double> > &marginals,
                     const std::vector<std::pair<size_t, size_t> > &blocks,
                     const param_set &params,
                     std::vector<std::vector<T> > &states,
                     std::vector<std::vector<double> > &posteriors) const;

  template <class T> void
  sample_states(const std::vector<size_t> &subtree_sizes,
                const std::vector<size_t> &parent_ids,
                const param_set &params,
                const std::vector<std::vector<double> > &marginals,
                const std::vector<std::pair<size_t, size_t> > &blocks,
                std::vector<std::vector<T> > &states) const;

  template <class T> void
  sample_states_sym(const std::vector<size_t> &subtree_sizes,
                    const std::vector<size_t> &parent_ids,
                    const param_set &params,
                    const std::vector<std::vector<double> > &marginals,
                    const std::vector<std::pair<size_t, size_t> > &blocks,
                    std::vector<std::vector<T> > &states) const;

  template <class T> void
  sample_states(const std::vector<size_t> &subtree_sizes,
                const std::vector<size_t> &parent_ids,
                const std::pair<double, double> &root_start_distr,
                const pair_state &root_distr,
                const std::vector<pair_state> &start_distr,
                const std::vector<triple_state> &triad_distr,
                const std::vector<std::vector<double> > &marginals,
                const std::vector<std::pair<size_t, size_t> > &blocks,
                std::vector<std::vector<T> > &states) const;

  //for general case
  template <class T> void
  sample_states(const std::vector<size_t> &subtree_sizes,
                const std::vector<size_t> &parent_ids,
                const std::vector<MSite> &sites,
                const size_t desert_size,
                const param_set &params,
                const std::vector<std::vector<double> > &marginals,
                const std::vector<std::vector<bool> > &marks,
                std::vector<std::vector<T> > &states) const;
  
private:

  template <class T> void
  accept_or_reject_proposal(const double log_ratio, const size_t idx,
                            std::vector<T> &states) const;
  template <class T> void
  sample_states_middle(const std::vector<size_t> &subtree_sizes,
                       const std::vector<size_t> &parent_ids,
                       const pair_state &logG,
                       const std::vector<triple_state> &logGP,
                       const std::vector<double> &marginals,
                       const std::vector<T> &states_prev,
                       std::vector<T> &states_curr,
                       const std::vector<T> &states_next,
                       const size_t node_id) const;

  template <class T> void
  sample_states_end(const std::vector<size_t> &subtree_sizes,
                    const std::vector<size_t> &parent_ids,
                    const pair_state &logG,
                    const std::vector<triple_state> &logGP,
                    const std::vector<double> &marginals,
                    const std::vector<T> &states_prev,
                    std::vector<T> &states_curr,
                    const size_t node_id) const;

  template <class T> void
  sample_states_start(const std::vector<size_t> &subtree_sizes,
                      const std::vector<size_t> &parent_ids,
                      const std::pair<double, double> &log_pi,
                      const pair_state &logG,
                      const std::vector<pair_state> &logP,
                      const std::vector<triple_state> &logGP,
                      const std::vector<double> &marginals,
                      std::vector<T> &states_curr,
                      const std::vector<T> &states_next,
                      const size_t node_id) const;

  template <class T> void
  sample_states_middle_sym(const std::vector<size_t> &subtree_sizes,
                           const std::vector<size_t> &parent_ids,
                           const pair_state &logG,
                           const std::vector<triple_state> &logGP,
                           const std::vector<double> &marginals,
                           const std::vector<T> &states_prev,
                           std::vector<T> &states_curr,
                           const std::vector<T> &states_next,
                           const size_t node_id) const;

  template <class T> void
  sample_states_end_sym(const std::vector<size_t> &subtree_sizes,
                        const std::vector<size_t> &parent_ids,
                        const std::pair<double, double> &log_pi,
                        const pair_state &logG,
                        const std::vector<pair_state> &logP,
                        const std::vector<triple_state> &logGP,
                        const std::vector<double> &marginals,
                        const std::vector<T> &states_prev,
                        std::vector<T> &states_curr,
                        const size_t node_id) const;
  
  template <class T> void
  sample_states_start_sym(const std::vector<size_t> &subtree_sizes,
                          const std::vector<size_t> &parent_ids,
                          const std::pair<double, double> &log_pi,
                          const pair_state &logG,
                          const std::vector<pair_state> &logP,
                          const std::vector<triple_state> &logGP,
                          const std::vector<double> &marginals,
                          std::vector<T> &states_curr,
                          const std::vector<T> &states_next,
                          const size_t node_id) const;
  
  template <class T> void
  sample_states_indep(const std::vector<size_t> &subtree_sizes,
                      const std::vector<size_t> &parent_ids,
                      const std::pair<double, double> &log_pi,
                      const std::vector<pair_state> &logP,
                      const std::vector<double> &marginals,
                      std::vector<T> &states_curr,
                      const size_t node_id) const;

  mutable std::uniform_real_distribution<> dis;
  mutable std::mt19937_64 gen;
  size_t max_iterations;
  size_t the_seed;
};


/*
  prev/curr/next denote site positions in the genome
  parent/self/child denote species (nodes) in the phylogenetic tree
*/
template <class T>
static double
MB_state_middle(const std::vector<size_t> &subtree_sizes,
                const std::vector<size_t> &parent_ids,
                const pair_state &logG,
                const std::vector<triple_state> &logGP,
                const std::vector<T> &meth_prev,
                const std::vector<T> &meth_curr,
                const std::vector<T> &meth_next, const size_t node_id) {

  const size_t parent = parent_ids[node_id];
  double lp0 = 0.0, lp1 = 0.0;

  // contribution of parent states in the network
  if (is_root(node_id)) {
    const size_t prev_self_state = meth_prev[node_id];
    lp0 += logG(prev_self_state, 0);
    lp1 += logG(prev_self_state, 1);
  } else {
    const T curr_parent_state = meth_curr[parent];
    const T prev_self_state = meth_prev[node_id];
    lp0 += logGP[node_id](prev_self_state, curr_parent_state, 0);
    lp1 += logGP[node_id](prev_self_state, curr_parent_state, 1);
  }

  // contribution of states in child species (includes network
  // children and parents of network children)
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const T prev_child_state = meth_prev[child_id];
    const T curr_child_state = meth_curr[child_id];
    lp0 += logGP[child_id](prev_child_state, 0, curr_child_state);
    lp1 += logGP[child_id](prev_child_state, 1, curr_child_state);
    count += subtree_sizes[child_id];
  }

  // contribution of states at next site (includes network children
  // and parents of network children)
  if (is_root(node_id)) {
    const T next_self_state = meth_next[node_id];
    lp0 += logG(0, next_self_state);
    lp1 += logG(1, next_self_state);
  } else {
    const T next_parent_state = meth_next[parent];
    const T next_self_state = meth_next[node_id];
    lp0 += logGP[node_id](0, next_parent_state, next_self_state);
    lp1 += logGP[node_id](1, next_parent_state, next_self_state);
  }
  return lp0 - lp1;
}


template <class T>
static double
MB_state_end(const std::vector<size_t> &subtree_sizes,
             const std::vector<size_t> &parent_ids,
             const pair_state &logG,
             const std::vector<triple_state> &logGP,
             const std::vector<T> &meth_prev,
             const std::vector<T> &meth_curr, const size_t node_id) {

  double lp0 = 0.0, lp1 = 0.0;
  // parents in network
  if (is_root(node_id)) {
    const size_t prev_self_state = meth_prev[node_id];
    lp0 += logG(prev_self_state, 0);
    lp1 += logG(prev_self_state, 1);
  } else {
    const T curr_parent_state = meth_curr[parent_ids[node_id]];
    const T prev_self_state = meth_prev[node_id];
    lp0 += logGP[node_id](prev_self_state, curr_parent_state, 0);
    lp1 += logGP[node_id](prev_self_state, curr_parent_state, 1);
  }

  // children in network
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const T prev_child_state = meth_prev[child_id];
    const T curr_child_state = meth_curr[child_id];
    lp0 += logGP[child_id](prev_child_state, 0, curr_child_state);
    lp1 += logGP[child_id](prev_child_state, 1, curr_child_state);
    count += subtree_sizes[child_id];
  }

  return lp0 - lp1;
}


template <class T>
static double
MB_state_start(const std::vector<size_t> &subtree_sizes,
               const std::vector<size_t> &parent_ids,
               const std::pair<double, double> &log_pi,
               const pair_state &logG,
               const std::vector<pair_state> &logP,
               const std::vector<triple_state> &logGP,
               const std::vector<T> &meth_curr,
               const std::vector<T> &meth_next, const size_t node_id) {

  double lp0 = 0.0, lp1 = 0.0;
  const size_t parent = parent_ids[node_id];
  // transition from parents in network [same here as parents in tree]
  // (or root distribution)
  if (is_root(node_id)) {
    lp0 += log_pi.first;
    lp1 += log_pi.second;
  } else {
    const T curr_parent_state = meth_curr[parent];
    lp0 += logP[node_id](curr_parent_state, 0);
    lp1 += logP[node_id](curr_parent_state, 1);
  }

  // children in network (only one parent; current self)
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const T curr_child_state = meth_curr[child_id];
    lp0 += logP[child_id](0, curr_child_state);
    lp1 += logP[child_id](1, curr_child_state);
    count += subtree_sizes[child_id];
  }

  // horizontal children in network and their other parents
  if (is_root(node_id)) {
    const T next_self_state = meth_next[node_id];
    lp0 += logG(0, next_self_state);
    lp1 += logG(1, next_self_state);
  } else {
    const T next_parent_state = meth_next[parent];
    const T next_self_state = meth_next[node_id];
    lp0 += logGP[node_id](0, next_parent_state, next_self_state);
    lp1 += logGP[node_id](1, next_parent_state, next_self_state);
  }
  return lp0 - lp1;
}


template <class T>
static double
MB_state_indep(const std::vector<size_t> &subtree_sizes,
               const std::vector<size_t> &parent_ids,
               const std::pair<double, double> &log_pi,
               const std::vector<pair_state> &logP,
               const std::vector<T> &meth_curr,
               const size_t node_id) {

  double lp0 = 0.0, lp1 = 0.0;
  const size_t parent = parent_ids[node_id];
  // transition from parents in network [same here as parents in tree]
  // (or root distribution)
  if (is_root(node_id)) {
    lp0 += log_pi.first;
    lp1 += log_pi.second;
  } else {
    const T curr_parent_state = meth_curr[parent];
    lp0 += logP[node_id](curr_parent_state, 0);
    lp1 += logP[node_id](curr_parent_state, 1);
  }

  // children in network (only one parent; current self)
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const T curr_child_state = meth_curr[child_id];
    lp0 += logP[child_id](0, curr_child_state);
    lp1 += logP[child_id](1, curr_child_state);
    count += subtree_sizes[child_id];
  }

  return lp0 - lp1;
}

// general case
// pi0 and G are shared by all nodes
template <class T>
static double
MB_state(const std::vector<size_t> &subtree_sizes,
         const std::vector<size_t> &parent_ids,
         const std::pair<double, double> log_pi,
         const pair_state &logG,
         const std::vector<pair_state> &logP,
         const std::vector<triple_state> &logGP,
         const std::vector<T> &meth_prev,
         const std::vector<T> &meth_curr,
         const std::vector<T> &meth_next,
         const std::vector<bool> &mark_prev,
         const std::vector<bool> &mark_curr,
         const std::vector<bool> &mark_next,
         const size_t node_id) {

  double lp0 = 0;
  double lp1 = 0;

  if (node_id == 0) { // root
    if (mark_prev[node_id]) {
      lp0 += logG(meth_prev[node_id], 0);
      lp1 += logG(meth_prev[node_id], 1);
    } else {
      lp0 += log_pi.first;
      lp1 += log_pi.second;
    }

    if (mark_next[node_id]) {
      lp0 += logG(0, meth_next[node_id]);
      lp1 += logG(1, meth_next[node_id]);
    }
  } else { // not root
    // between par-level and self-level
    if (mark_prev[node_id] && mark_curr[parent_ids[node_id]]) {
      lp0 += logGP[node_id](meth_prev[node_id],
                            meth_curr[parent_ids[node_id]], 0);
      lp1 += logGP[node_id](meth_prev[node_id],
                            meth_curr[parent_ids[node_id]], 1);
    } else if (mark_prev[node_id]) {
      lp0 += logG(meth_prev[node_id], 0);
      lp1 += logG(meth_prev[node_id], 1);
    } else if (mark_curr[parent_ids[node_id]]) {
      lp0 += logP[node_id](meth_curr[parent_ids[node_id]], 0);
      lp1 += logP[node_id](meth_curr[parent_ids[node_id]], 1);
    } else {
      lp0 += log_pi.first;
      lp1 += log_pi.second;
    }

    if (mark_next[node_id] && mark_next[parent_ids[node_id]]) {
      lp0 += logGP[node_id](0, meth_next[parent_ids[node_id]],
                            meth_next[node_id]);
      lp1 += logGP[node_id](1, meth_next[parent_ids[node_id]],
                            meth_next[node_id]);
    } else if (mark_next[node_id]) {
      lp0 += logG(0, meth_next[node_id]);
      lp1 += logG(1, meth_next[node_id]);
    }
  }

  if (!is_leaf(subtree_sizes[node_id])) {  // between self-level and child-level
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      if (mark_prev[child_id]) {
        lp0 += logGP[child_id](meth_prev[child_id], 0, meth_curr[child_id]);
        lp1 += logGP[child_id](meth_prev[child_id], 1, meth_curr[child_id]);
      } else {
        lp0 += logP[child_id](0, meth_curr[child_id]);
        lp1 += logP[child_id](1, meth_curr[child_id]);
      }
      count += subtree_sizes[child_id];
    }
  }

  return lp0 - lp1;
}


template <class T>
void
epiphy_mcmc::accept_or_reject_proposal(const double log_ratio, const size_t idx,
                                       std::vector<T> &states) const {
  const T tmp_state = states[idx];
  const double ratio = (!tmp_state) ? exp(-log_ratio) : exp(log_ratio);
  if (dis(gen) < ratio)  {
    // ADS: need to make sure the "!" is defined the type T
    states[idx] = !tmp_state; // = (1 - tmp_state)
  }
}


template <class T>
void
epiphy_mcmc::sample_states_middle(const std::vector<size_t> &subtree_sizes,
                                  const std::vector<size_t> &parent_ids,
                                  const pair_state &logG,
                                  const std::vector<triple_state> &logGP,
                                  const std::vector<double> &marginals,
                                  const std::vector<T> &states_prev,
                                  std::vector<T> &states_curr,
                                  const std::vector<T> &states_next,
                                  const size_t node_id) const {

  // ADS: this test might not need to check whether or not the node is
  // a leaf, because we only have data at leaf nodes. So the check for
  // missing data might be enough.
  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else { // this happens if internal node or missing data at leaf
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      sample_states_middle(subtree_sizes, parent_ids, logG, logGP, marginals,
                           states_prev, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_middle(subtree_sizes, parent_ids, logG, logGP,
                      states_prev, states_curr, states_next, node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


template <class T>
void
epiphy_mcmc::sample_states_end(const std::vector<size_t> &subtree_sizes,
                               const std::vector<size_t> &parent_ids,
                               const pair_state &logG,
                               const std::vector<triple_state> &logGP,
                               const std::vector<double> &marginals,
                               const std::vector<T> &states_prev,
                               std::vector<T> &states_curr,
                               const size_t node_id) const {

  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      sample_states_end(subtree_sizes, parent_ids, logG, logGP,
                        marginals, states_prev, states_curr, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_end(subtree_sizes, parent_ids, logG, logGP,
                   states_prev, states_curr, node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


template <class T>
void
epiphy_mcmc::sample_states_start(const std::vector<size_t> &subtree_sizes,
                                 const std::vector<size_t> &parent_ids,
                                 const std::pair<double, double> &log_pi,
                                 const pair_state &logG,
                                 const std::vector<pair_state> &logP,
                                 const std::vector<triple_state> &logGP,
                                 const std::vector<double> &marginals,
                                 std::vector<T> &states_curr,
                                 const std::vector<T> &states_next,
                                 const size_t node_id) const {

  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      sample_states_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                          marginals, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                     states_curr, states_next, node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  symmetric sampling  /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T>
void
epiphy_mcmc::sample_states_middle_sym(const std::vector<size_t> &subtree_sizes,
                                      const std::vector<size_t> &parent_ids,
                                      const pair_state &logG,
                                      const std::vector<triple_state> &logGP,
                                      const std::vector<double> &marginals,
                                      const std::vector<T> &states_prev,
                                      std::vector<T> &states_curr,
                                      const std::vector<T> &states_next,
                                      const size_t node_id) const {

  // ADS: this test might not need to check whether or not the node is
  // a leaf, because we only have data at leaf nodes. So the check for
  // missing data might be enough.
  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else { // this happens if internal node or missing data at leaf
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      sample_states_middle(subtree_sizes, parent_ids, logG, logGP, marginals,
                           states_prev, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_middle(subtree_sizes, parent_ids, logG, logGP,
                      states_prev, states_curr, states_next, node_id);
    const double log_ratio_sym = log_ratio +
      MB_state_middle(subtree_sizes, parent_ids, logG, logGP,
                      states_next, states_curr, states_prev, node_id);

    accept_or_reject_proposal(log_ratio_sym, node_id, states_curr);
  }
}

template <class T>
void
epiphy_mcmc::sample_states_end_sym(const std::vector<size_t> &subtree_sizes,
                                   const std::vector<size_t> &parent_ids,
                                   const std::pair<double, double> &log_pi,
                                   const pair_state &logG,
                                   const std::vector<pair_state> &logP,
                                   const std::vector<triple_state> &logGP,
                                   const std::vector<double> &marginals,
                                   const std::vector<T> &states_prev,
                                   std::vector<T> &states_curr,
                                   const size_t node_id) const {
  
  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      sample_states_end(subtree_sizes, parent_ids, logG, logGP,
                        marginals, states_prev, states_curr, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_end(subtree_sizes, parent_ids, logG, logGP,
                   states_prev, states_curr, node_id);

    const double log_ratio_sym = log_ratio + 
      MB_state_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                     states_curr, states_prev, node_id);
    accept_or_reject_proposal(log_ratio_sym, node_id, states_curr);
  }
}


template <class T>
void
epiphy_mcmc::sample_states_start_sym(const std::vector<size_t> &subtree_sizes,
                                     const std::vector<size_t> &parent_ids,
                                     const std::pair<double, double> &log_pi,
                                     const pair_state &logG,
                                     const std::vector<pair_state> &logP,
                                     const std::vector<triple_state> &logGP,
                                     const std::vector<double> &marginals,
                                     std::vector<T> &states_curr,
                                     const std::vector<T> &states_next,
                                     const size_t node_id) const {
  
  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      sample_states_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                          marginals, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                     states_curr, states_next, node_id);
    const double log_ratio_sym = log_ratio +
      MB_state_end(subtree_sizes, parent_ids, logG, logGP,
                   states_next, states_curr, node_id);
  
    accept_or_reject_proposal(log_ratio_sym, node_id, states_curr);
  }
}



////////////////////////////////////////////////////////////////////////////////
///////////////////////////// independent sampling /////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T>
void
epiphy_mcmc::sample_states_indep(const std::vector<size_t> &subtree_sizes,
                                 const std::vector<size_t> &parent_ids,
                                 const std::pair<double, double> &log_pi,
                                 const std::vector<pair_state> &logP,
                                 const std::vector<double> &marginals,
                                 std::vector<T> &states_curr,
                                 const size_t node_id) const {

  if (is_leaf(subtree_sizes[node_id]) &&
      !missing_meth_value(marginals[node_id])) {
    states_curr[node_id] = dis(gen) > marginals[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      sample_states_indep(subtree_sizes, parent_ids, log_pi, logP,
                          marginals, states_curr, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_state_indep(subtree_sizes, parent_ids, log_pi, logP,
                     states_curr, node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


template <class T>
void
epiphy_mcmc::sample_states(const std::vector<size_t> &subtree_sizes,
                           const std::vector<size_t> &parent_ids,
                           const param_set &params,
                           const std::vector<std::vector<double> > &marginals,
                           const std::vector<std::pair<size_t, size_t> > &blocks,
                           std::vector<std::vector<T> > &states) const {

  std::pair<double, double> log_pi(std::make_pair(log(params.pi0),
                                                  log(1.0 - params.pi0)));
  pair_state logG(params.g0, 1.0 - params.g0, 1.0 - params.g1, params.g1);
  logG.make_logs();

  std::vector<pair_state> logP;
  std::vector<triple_state> logGP;
  get_transition_matrices(params, logP, logGP);
  for (size_t i = 0; i < logP.size(); ++i) {
    logP[i].make_logs();
    logGP[i].make_logs();
  }

  for (size_t i = 0; i < blocks.size(); ++i) {
    const size_t start = blocks[i].first;
    const size_t end = blocks[i].second;
    if (start < end) {
      sample_states_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                          marginals[start], states[start], states[start + 1], 0);
      for (size_t pos = start + 1; pos < end; ++pos)
        sample_states_middle(subtree_sizes, parent_ids, logG, logGP,
                             marginals[pos], states[pos - 1], states[pos],
                             states[pos + 1], 0);
      sample_states_end(subtree_sizes, parent_ids, logG, logGP, marginals[end],
                        states[end - 1], states[end], 0);
    } else {
      sample_states_indep(subtree_sizes, parent_ids, log_pi, logP,
                          marginals[start], states[start], 0);
    }
  }
}


template <class T>
void
epiphy_mcmc::sample_states_sym(const std::vector<size_t> &subtree_sizes,
                               const std::vector<size_t> &parent_ids,
                               const param_set &params,
                               const std::vector<std::vector<double> > &marginals,
                               const std::vector<std::pair<size_t, size_t> > &blocks,
                               std::vector<std::vector<T> > &states) const {
  
  std::pair<double, double> log_pi(std::make_pair(log(params.pi0),
                                                  log(1.0 - params.pi0)));
  pair_state logG(params.g0, 1.0 - params.g0, 1.0 - params.g1, params.g1);
  logG.make_logs();

  std::vector<pair_state> logP;
  std::vector<triple_state> logGP;
  get_transition_matrices(params, logP, logGP);
  for (size_t i = 0; i < logP.size(); ++i) {
    logP[i].make_logs();
    logGP[i].make_logs();
  }

  for (size_t i = 0; i < blocks.size(); ++i) {
    const size_t start = blocks[i].first;
    const size_t end = blocks[i].second;
    if (start < end) {
      sample_states_start_sym(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                              marginals[start], states[start], states[start + 1], 0);
      for (size_t pos = start + 1; pos < end; ++pos)
        sample_states_middle_sym(subtree_sizes, parent_ids, logG, logGP,
                                 marginals[pos], states[pos - 1], states[pos],
                                 states[pos + 1], 0);
      sample_states_end_sym(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                            marginals[end], states[end - 1], states[end], 0);
    } else {
      sample_states_indep(subtree_sizes, parent_ids, log_pi, logP,
                          marginals[start], states[start], 0);
    }
  }
}

// construct transition probabilities directly from sufficient statistics
template <class T>
void
epiphy_mcmc::sample_states(const std::vector<size_t> &subtree_sizes,
                           const std::vector<size_t> &parent_ids,
                           const std::pair<double, double> &root_start_distr,
                           const pair_state &root_distr,
                           const std::vector<pair_state> &start_distr,
                           const std::vector<triple_state> &triad_distr,
                           const std::vector<std::vector<double> > &marginals,
                           const std::vector<std::pair<size_t, size_t> > &blocks,
                           std::vector<std::vector<T> > &states) const {

  const double denom = root_start_distr.first + root_start_distr.second;
  std::pair<double, double> log_pi =
    std::make_pair(log(root_start_distr.first/denom),
                   log(root_start_distr.second/denom));
  pair_state logG(root_distr);
  logG.to_probabilities();
  logG.make_logs();

  std::vector<pair_state> logP;
  std::vector<triple_state> logGP;
  for (size_t i = 0; i < start_distr.size(); ++i) {
    logP.push_back(start_distr[i]);
    logP[i].to_probabilities();
    logP[i].make_logs();
    logGP.push_back(triad_distr[i]);
    logGP[i].to_probabilities();
    logGP[i].make_logs();
  }

  for (size_t i = 0; i < blocks.size(); ++i) {
    const size_t start = blocks[i].first;
    const size_t end = blocks[i].second;
    if (start < end) {
      sample_states_start(subtree_sizes, parent_ids, log_pi, logG, logP, logGP,
                          marginals[start], states[start], states[start + 1], 0);
      for (size_t pos = start + 1; pos < end; ++pos)
        sample_states_middle(subtree_sizes, parent_ids, logG, logGP,
                             marginals[pos], states[pos - 1], states[pos],
                             states[pos + 1], 0);
      sample_states_end(subtree_sizes, parent_ids, logG, logGP, marginals[end],
                        states[end - 1], states[end], 0);
    } else {
      sample_states_indep(subtree_sizes, parent_ids, log_pi, logP,
                          marginals[end], states[start], 0);
    }
  }
}


//for general case
template <class T> void
epiphy_mcmc::sample_states(const std::vector<size_t> &subtree_sizes,
                           const std::vector<size_t> &parent_ids,
                           const std::vector<MSite> &sites,
                           const size_t desert_size,
                           const param_set &params,
                           const std::vector<std::vector<double> > &marginals,
                           const std::vector<std::vector<bool> > &marks,
                           std::vector<std::vector<T> > &states) const {

  std::pair<double, double> log_pi(std::make_pair(log(params.pi0),
                                                  log(1.0 - params.pi0)));
  pair_state logG(params.g0, 1.0 - params.g0, 1.0 - params.g1, params.g1);
  logG.make_logs();

  std::vector<pair_state> logP;
  std::vector<triple_state> logGP;
  get_transition_matrices(params, logP, logGP);
  for (size_t i = 0; i < logP.size(); ++i) {
    logP[i].make_logs();
    logGP[i].make_logs();
  }

  std::vector<bool> mark_prev;
  std::vector<bool> mark_next;
  std::vector<T> meth_prev;
  std::vector<T> meth_next;
  const size_t n_sites = marginals.size();
  const size_t n_nodes = subtree_sizes.size();
  for (size_t i = 0; i < n_sites; ++i) { // across all sites
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size) {
      // has no prev neighbor
      mark_prev = std::vector<bool>(n_nodes, false);
      meth_prev = std::vector<T>(n_nodes);
    } else {
      mark_prev = marks[i - 1];
      meth_prev = states[i - 1];
    }

    if (i == n_sites - 1 || distance(sites[i], sites[i+1]) > desert_size) {
      mark_next = std::vector<bool>(n_nodes, false);
      meth_next = std::vector<T>(n_nodes);
    } else {
      mark_next = marks[i + 1];
      meth_next = states[i + 1];
    }

    for (size_t node_id = 0; node_id < n_nodes; ++node_id) { // across all nodes
      if (!missing_meth_value(marginals[i][node_id])) {
        states[i][node_id] = dis(gen) > marginals[i][node_id] ? 0 : 1;
      } else {
        const double log_ratio =
          MB_state(subtree_sizes, parent_ids,
                   log_pi, logG, logP, logGP,
                   meth_prev, states[i], meth_next,
                   mark_prev, marks[i], mark_next, node_id);

        accept_or_reject_proposal(log_ratio, node_id, states[i]);
      }
    }
  }
}



struct mcmc_stat{
  std::pair<double, double> root_start_distr;
  pair_state root_distr;
  std::vector<pair_state> start_distr;
  std::vector<triple_state> triad_distr;

  mcmc_stat() {}

  mcmc_stat(const std::pair<double, double> &root_start,
            const pair_state &root,
            const std::vector<pair_state> &start,
            const std::vector<triple_state> &triad) :
    root_start_distr(root_start), root_distr(root),
    start_distr(start), triad_distr(triad) {}

  void scale();
};

void
sum(const std::vector<mcmc_stat> &mcmcstats,
    mcmc_stat &ave_mcmc_stat);


// measure MCMC convergence
// outter vector for parallel chains
// inner vector for within-chain sample stats
void
EPSR(std::vector<std::vector<mcmc_stat> > &mcmcstats,
     vector<double> &epsr);

#endif
