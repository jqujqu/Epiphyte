/*    Copyright (C) 2015 University of Southern California and
 *                       Andrew D. Smith and Jenny Qu
 *
 *    Authors: Jenny Qu and Andrew Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/******************************************************************************/

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>    //std::accumulate
#include <random>
#include <algorithm>  //std::max, min
#include <cmath>      //std::abs
#include <limits>     //std::numeric_limits
#include <iterator>   //std::distance
#include <unistd.h>

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

/* from methpipe */
#include "MethpipeFiles.hpp"
#include "MethpipeSite.hpp"

/* headers for epigenomic evolution */
#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"
#include "param_set.hpp"
#include "epiphy_utils.hpp"
#include "sufficient_statistics_helpers.hpp"
#include "optimize_params.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::pair;
using std::make_pair;
using std::abs;
using std::numeric_limits;
using std::accumulate;
using std::inner_product;
using std::min;
using std::max;
using std::istringstream;
using std::to_string;

using std::ostream_iterator;

#include <functional>
using std::placeholders::_1;
using std::bind;
using std::plus;

static const double PROBABILITY_GUARD = 1e-10;
static const double POST_CONV_TOL = 1e-6;  //MAGIC
static const double KL_CONV_TOL = 1e-8;    //MAGIC

#include <random>
using std::uniform_real_distribution;
std::random_device rd; //seed generator
//std::mt19937_64 gen(rd()); //generator initialized with seed from rd
std::mt19937_64 gen(0); //generator initialized with seed from rd


static size_t
distance(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom ? max(a.pos, b.pos) - min(a.pos, b.pos) :
    std::numeric_limits<size_t>::max();
}


static void
separate_regions(const size_t desert_size, vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &reset_points) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      reset_points.push_back(std::make_pair(i, i));
    else reset_points.back().second = i;
}


// putting this function here because it is related to the parameters
static double
estimate_pi0(const vector<vector<double> > &meth) {
  double total = 0.0;
  size_t N = 0;
  for (auto i(meth.begin()); i != meth.end(); ++i)
    for (auto j(i->begin()); j != i->end(); ++j)
      if (!missing_meth_value(*j)) {
        total += *j;
        ++N;
      }
  return (N - total)/N;
}

static void
estimate_g0_g1(const vector<vector<double> > &meth,
               double &g0_est, double &g1_est) {
  const size_t n_leaves = meth[0].size();
  double g00 = 0.0, g01 = 0.0, g10 = 0.0, g11 = 0.0;
  for (size_t i = 0; i < n_leaves; ++i)
    for (size_t j = 0; j < meth.size() - 1; ++j)
      if (!missing_meth_value(meth[j][i]) && !missing_meth_value(meth[j + 1][i])) {
        g00 += (1.0 - meth[j][i])*(1.0 - meth[j + 1][i]);
        g01 += (1.0 - meth[j][i])*meth[j + 1][i];
        g10 += meth[j][i]*(1.0 - meth[j + 1][i]);
        g11 += meth[j][i]*meth[j + 1][i];
      }

  g0_est = (g00 + 1.0)/(g01 + g00 + 2.0);
  g1_est = (g11 + 1.0)/(g11 + g10 + 2.0);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////    FOR LOADING DATA     //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template <class T> bool
parse_line(const string &line,
           vector<MSite> &sites, vector<vector<T> > &states) {

  std::istringstream iss(line);

  sites.push_back(MSite());
  if (!(iss >> sites.back().chrom >> sites.back().pos))
    return false;

  states.push_back(vector<T>(std::istream_iterator<T>(iss),
                             std::istream_iterator<T>()));

  return true;
}


static void
read_meth_table(const string &table_file,
                vector<MSite> &sites,
                vector<string> &species_names,
                vector<vector<double> > &states) {

  std::ifstream table_in(table_file.c_str());
  if (!table_in)
    throw std::runtime_error("bad table file: " + table_file);

  string line;
  getline(table_in, line);
  std::istringstream iss(line);
  copy(std::istream_iterator<string>(iss), std::istream_iterator<string>(),
       std::back_inserter(species_names));

  while (getline(table_in, line))
    if (line[0] != '#')
      if (!parse_line(line, sites, states) ||
          states.back().size() != species_names.size())
        throw std::runtime_error("bad table file line: " + line);
}


static void
kl_divergence(const vector<triple_state> &P, const vector<triple_state> &Q,
              vector<double> &kld) {
  assert(P.size() == Q.size());
  kld.resize(P.size(), 0.0); // clearing not needed; values not used here
  for (size_t i = 0; i < P.size(); ++i) {
    vector<double> p, q;
    P[i].flatten(p);
    Q[i].flatten(q);
    transform(p.begin(), p.end(), p.begin(),
              bind(plus<double>(), _1, PROBABILITY_GUARD));
    transform(q.begin(), q.end(), q.begin(),
              bind(plus<double>(), _1, PROBABILITY_GUARD));
    kld[i] = kl_divergence(p, q);
  }
}


static void
leaf_to_tree_prob(const vector<size_t> &subtree_sizes,
                  const vector<vector<double> > &leaf_prob_table,
                  vector<vector<double> > &tree_prob_table) {
  // ADS: should the size of

  // copy leaves first
  vector<size_t> leaves_preorder;
  subtree_sizes_to_leaves_preorder(subtree_sizes, leaves_preorder);

  for (size_t i = 0; i < leaf_prob_table.size(); ++i)
    for (size_t j = 0; j < leaf_prob_table[i].size(); ++j)
      if (valid_probability(leaf_prob_table[i][j]))
        tree_prob_table[i][leaves_preorder[j]] =
          away_from_extremes(leaf_prob_table[i][j], PROBABILITY_GUARD);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////   Markov blanket probabilities    ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <class T> static double
MB_middle_site(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
               const pair_state &logG, const vector<triple_state> &logGP,
               const vector<T> &meth_prev,
               const vector<T> &meth_curr,
               const vector<T> &meth_next, const size_t node_id) {

  const size_t parent = parent_ids[node_id];
  double lp0 = 0.0, lp1 = 0.0;

  // contribution of parent states in the network
  if (is_root(node_id)) {
    const size_t prev_self_state = meth_prev[node_id];
    lp0 += logG(prev_self_state, 0);
    lp1 += logG(prev_self_state, 1);
  }
  else {
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
    const T next_curr_state = meth_next[node_id];
    lp0 += logG(0, next_curr_state);
    lp1 += logG(1, next_curr_state);
  }
  else {
    const T next_parent_state = meth_next[parent];
    const T next_curr_state = meth_next[node_id];
    lp0 += logGP[node_id](0, next_parent_state, next_curr_state);
    lp1 += logGP[node_id](1, next_parent_state, next_curr_state);
  }
  return lp0 - lp1;
}


// ADS: why is pi0 only used for start nodes?
template <class T> static double
MB_end_site(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
            const pair_state &logG, const vector<triple_state> &logGP,
            const vector<T> &meth_prev,
            const vector<T> &meth_curr, const size_t node_id) {

  double lp0 = 0.0, lp1 = 0.0;
  // parents in network
  if (is_root(node_id)) {
    const size_t prev_self_state = meth_prev[node_id];
    lp0 += logG(prev_self_state, 0);
    lp1 += logG(prev_self_state, 1);
  }
  else {
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


template <class T> static double
MB_start_site(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
              const double pi0, const pair_state &logG,
              const vector<pair_state> &logP, const vector<triple_state> &logGP,
              const vector<T> &meth_curr,
              const vector<T> &meth_next, const size_t node_id) {

  double lp0 = 0.0, lp1 = 0.0;
  const size_t parent = parent_ids[node_id];
  // transition from parents in network [same here as parents in tree]
  // (or root distribution)
  if (is_root(node_id)) {
    lp0 += log(pi0);
    lp1 += log(1.0 - pi0);
  }
  else {
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
  }
  else {
    const T next_parent_state = meth_next[parent];
    const T next_curr_state = meth_next[node_id];
    lp0 += logGP[node_id](0, next_parent_state, next_curr_state);
    lp1 += logGP[node_id](1, next_parent_state, next_curr_state);
  }
  return lp0 - lp1;
}


template <class T>
static void
accept_or_reject_proposal(const double log_ratio, const size_t idx,
                          vector<T> &states) {
  const T tmp_state = states[idx];
  const double ratio = (!tmp_state) ? exp(-log_ratio) : exp(log_ratio);
  uniform_real_distribution<> dis(0, 1);  // interval [0, 1)
  if (dis(gen) < ratio)
    states[idx] = !tmp_state; // 1 - tmp_state;
}


template <class T>
static void
MH_update_site_middle(const vector<size_t> &subtree_sizes,
                      const vector<size_t> &parent_ids,
                      const pair_state &logG,
                      const vector<triple_state> &logGP,
                      const vector<double> &probs_table,
                      const vector<T> &states_prev,
                      vector<T> &states_curr,
                      const vector<T> &states_next,
                      const size_t node_id) {

  /* if leaf observed, sample from the observed probability*/
  if (is_leaf(subtree_sizes[node_id]) && probs_table[node_id] >= 0.0) {
    // ADS: fix the way leaf states are updated when their value is
    // missing
    uniform_real_distribution<> dis(0, 1); // ADS: remember to fix this
    states_curr[node_id] = dis(gen) > probs_table[node_id] ? 0 : 1;
  }
  else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      MH_update_site_middle(subtree_sizes, parent_ids, logG, logGP, probs_table,
                            states_prev, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio = MB_middle_site(subtree_sizes, parent_ids, logG, logGP,
                                            states_prev, states_curr, states_next,
                                            node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


template <class T>
static void
MH_update_site_end(const vector<size_t> &subtree_sizes,
                   const vector<size_t> &parent_ids,
                   const pair_state &logG,
                   const vector<triple_state> &logGP,
                   const vector<double> &probs_table,
                   const vector<T> &states_prev,
                   vector<T> &states_curr,
                   const size_t node_id) {

  /* if leaf observed, sample from the observed probability*/
  if (is_leaf(subtree_sizes[node_id]) && probs_table[node_id] >= 0) {
    uniform_real_distribution<> dis(0, 1); // ADS: remember to fix this
    states_curr[node_id] = dis(gen) > probs_table[node_id] ? 0 : 1;
  }
  else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      MH_update_site_end(subtree_sizes, parent_ids, logG, logGP,
                         probs_table, states_prev, states_curr, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio = MB_end_site(subtree_sizes, parent_ids, logG, logGP,
                                         states_prev, states_curr, node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


template <class T> static void
MH_update_site_start(const vector<size_t> &subtree_sizes,
                     const vector<size_t> &parent_ids,
                     const double pi0,
                     const pair_state &logG,
                     const vector<pair_state> &logP,
                     const vector<triple_state> &logGP,
                     const vector<double> &probs_table,
                     vector<T> &states_curr,
                     const vector<T> &states_next,
                     const size_t node_id) {

  /* if leaf observed, sample from the observed probability*/
  if (is_leaf(subtree_sizes[node_id]) && probs_table[node_id] >= 0) {
    uniform_real_distribution<> dis(0, 1); // ADS: remember to fix this
    states_curr[node_id] = dis(gen) > probs_table[node_id] ? 0 : 1;
  }
  else {
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      MH_update_site_start(subtree_sizes, parent_ids, pi0, logG, logP, logGP,
                           probs_table, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_start_site(subtree_sizes, parent_ids, pi0, logG, logP, logGP,
                    states_curr, states_next, node_id);

    accept_or_reject_proposal(log_ratio, node_id, states_curr);
  }
}


template <class T> static void
MH_update(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
          const param_set &ps, const vector<pair<size_t, size_t> > &reset_points,
          const vector<vector<double> > &probs_table, // only used at leaves
          vector<vector<T> > &states_table) {

  pair_state logG(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);
  logG.make_logs();

  vector<pair_state> logP;
  vector<triple_state> logGP;
  get_transition_matrices(ps, logP, logGP);
  for (size_t i = 0; i < logP.size(); ++i) {
    logP[i].make_logs();
    logGP[i].make_logs();
  }

  for (size_t i = 0; i < reset_points.size(); ++i) {
    const size_t start = reset_points[i].first;
    const size_t end = reset_points[i].second; // !!!!
    if (start < end) {
      MH_update_site_start(subtree_sizes, parent_ids, ps.pi0, logG, logP, logGP,
                           probs_table[start], states_table[start],
                           states_table[start + 1], 0);
      for (size_t pos = start + 1; pos < end; ++pos)
        MH_update_site_middle(subtree_sizes, parent_ids, logG, logGP,
                              probs_table[pos], states_table[pos - 1],
                              states_table[pos], states_table[pos + 1], 0);
      MH_update_site_end(subtree_sizes, parent_ids, logG, logGP, probs_table[end],
                         states_table[end - 1], states_table[end], 0);
    }
  }
}


static void
count_triads(const vector<size_t> &subtree_sizes,
             const vector<size_t> &parent_ids,
             const vector<vector<bool> > &tree_state_table,
             const vector<pair<size_t, size_t> > &reset_points,
             pair<double, double> &root_start_counts,
             pair_state &root_counts,
             vector<pair_state> &start_counts,
             vector<triple_state> &triad_counts) {

  root_start_counts = make_pair(0.0, 0.0);
  triad_counts = vector<triple_state>(subtree_sizes.size());
  start_counts = vector<pair_state>(subtree_sizes.size());
  root_counts = pair_state();

  for (size_t i = 0; i < reset_points.size(); ++i) {

    const size_t start = reset_points[i].first;
    const size_t end = reset_points[i].second;

    root_start_counts.first += !tree_state_table[start][0];
    root_start_counts.second += tree_state_table[start][0];

    for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
      const size_t parent = tree_state_table[start][parent_ids[node_id]];
      const size_t curr = tree_state_table[start][node_id];
      start_counts[node_id](parent, curr) += 1.0;
    }

    for (size_t pos = start + 1; pos <= end; ++pos) {

      root_counts(tree_state_table[pos - 1][0], tree_state_table[pos][0])++;

      for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
        const size_t parent = tree_state_table[pos][parent_ids[node_id]];
        const size_t prev = tree_state_table[pos - 1][node_id];
        const size_t curr = tree_state_table[pos][node_id];
        triad_counts[node_id](prev, parent, curr) += 1.0;
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////// EXPECTATION MAXIMIZATION ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
maximization_step(const bool VERBOSE, const size_t MAXITER,
                  const vector<size_t> &subtree_sizes,
                  const pair<double, double> &root_start_counts,
                  const pair_state &root_counts,
                  const vector<pair_state> &start_counts,
                  const vector<triple_state> &triad_counts,
                  param_set &ps, PhyloTreePreorder &t) {

  // one M-step: optimize parameters
  double diff = std::numeric_limits<double>::max();
  for (size_t iter = 0; iter < MAXITER && diff > param_set::tolerance; ++iter) {
    if (VERBOSE)
      cerr << "[inside maximization: iter=" << iter << "]" << endl
           << ps << endl;

    param_set prev_ps(ps);
    optimize_params(VERBOSE, subtree_sizes, root_start_counts,
                    root_counts, start_counts, triad_counts, ps);

    diff = param_set::absolute_difference(ps, prev_ps); // convergence criteria
  }
}


template <class T>
static void
expectation_step(const bool VERBOSE, const size_t mh_max_iterations,
                 const size_t max_app_iter,
                 const vector<size_t> &subtree_sizes,
                 const vector<size_t> &parent_ids,
                 const vector<pair<size_t, size_t> > &reset_points,
                 const param_set &ps,
                 pair<double, double> &root_start_counts,
                 pair_state &root_counts,
                 vector<pair_state> &start_counts,
                 vector<triple_state> &triad_counts,
                 vector<vector<double> > &tree_prob_table, // only use leaf vals?
                 vector<vector<T> > &sampled_states) {

  const size_t n_nodes = subtree_sizes.size();

  root_start_counts = make_pair(0.0, 0.0);
  triad_counts = vector<triple_state>(n_nodes);
  start_counts = vector<pair_state>(n_nodes);
  vector<triple_state> triad_counts_prev(n_nodes);

  bool converged = false;
  size_t mh_iter = 0;
  for (mh_iter = 0; mh_iter < mh_max_iterations && !converged; ++mh_iter) {
    if (VERBOSE)
      cerr << "\r[inside expectation: M-H (iter=" << mh_iter << ")]";

    pair<double, double> root_start_counts_samp;
    pair_state root_counts_samp;
    vector<pair_state> start_counts_samp(n_nodes);
    vector<triple_state> triad_counts_samp(n_nodes);
    count_triads(subtree_sizes, parent_ids, sampled_states, reset_points,
                 root_start_counts_samp, root_counts_samp, start_counts_samp,
                 triad_counts_samp);

    // update the counts with current iteration samples
    root_start_counts.first += root_start_counts_samp.first;
    root_start_counts.second += root_start_counts_samp.second;
    root_counts = root_counts + root_counts_samp;
    for (size_t i = 0; i < triad_counts.size(); ++i) {
      triad_counts[i] = triad_counts[i] + triad_counts_samp[i];
      start_counts[i] = start_counts[i] + start_counts_samp[i];
    }

    // check how far the proportions have moved
    // only the triad_counts are used for convergence testing; the
    // start and root counts do not contribute enough to check
    vector<double> divergence;
    kl_divergence(triad_counts_prev, triad_counts, divergence);

    // determine convergence based on how far proportions have moved
    converged =
      *max_element(divergence.begin(), divergence.end()) < KL_CONV_TOL;

    if (!converged) // take next sample (all sites)
      // tree_prob_table only used at leaf nodes (to sample leaf states)
      MH_update(subtree_sizes, parent_ids, ps, reset_points,
                tree_prob_table, sampled_states);

    // retain previous values for comparison with next iteration
    triad_counts_prev = triad_counts; // ADS: should this be a swap?
  }
  if (VERBOSE)
    cerr << endl;

  root_start_counts.first /= mh_iter;
  root_start_counts.second /= mh_iter;
  root_counts.div(mh_iter);
  for (size_t i = 0; i < triad_counts.size(); ++i) {
    start_counts[i].div(mh_iter);
    triad_counts[i].div(mh_iter);
  }

  if (VERBOSE) {
    cerr << "[MH iterations=" << mh_iter << ']' << endl
         << "triad_counts:" << endl;
    copy(triad_counts.begin(), triad_counts.end(),
         ostream_iterator<triple_state>(cerr, "\n"));
    // ADS: should we print more summary statistics here?
  }
}


template <class T> static void
mcmc_estimate_posterior(const bool VERBOSE, const size_t mh_max_iterations,
                        const vector<size_t> &subtree_sizes,
                        const vector<size_t> &parent_ids,
                        const vector<pair<size_t, size_t> > &reset_points,
                        const param_set &ps,
                        vector<vector<T> > &sampled_states, // has value to start
                        vector<vector<double> > &probs_table, // ADS: values for each leaf ???????
                        vector<vector<double> > &posteriors) {

  const size_t n_nodes = subtree_sizes.size();

  posteriors = vector<vector<double> >(probs_table.size(),
                                       vector<double>(n_nodes, 0.0));

  for (size_t mh_iter = 0; mh_iter < mh_max_iterations; ++mh_iter) {
    if (VERBOSE)
      cerr << "\r[inside mcmc_estimate_posterior (iter=" << mh_iter << ")]";

    MH_update(subtree_sizes, parent_ids, ps, reset_points,
              probs_table, sampled_states);

    for (size_t i = 0; i < sampled_states.size(); ++i)
      for (size_t j = 0; j < sampled_states[i].size(); ++j)
        posteriors[i][j] += sampled_states[i][j];
  }
  if (VERBOSE)
    cerr << endl;

  for (size_t i = 0; i < posteriors.size(); ++i)
    for (size_t j = 0; j < posteriors[i].size(); ++j)
      posteriors[i][j] /= mh_max_iterations;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////// FORMATTING AND WRITING OUTPUT /////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


static void
local_parsimony(const vector<size_t> &subtree_sizes,
                const vector<size_t> &parent_ids,
                const size_t node_id, string &state) {
  if (!is_leaf(subtree_sizes[node_id])) {

    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      local_parsimony(subtree_sizes, parent_ids, child_id, state);
      count += subtree_sizes[child_id];
    }

    char s = state[node_id + 1];
    bool SAME = true;
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      if (state[child_id] != s) {
        SAME = false;
        break; // ADS: "break" should never be used
      }
      count += subtree_sizes[child_id];
    }
    if (SAME) {
      state[node_id] = s;
    }
    else {
      s = state[parent_ids[node_id]];
      state[node_id] = s;
    }
  }
}

static void
tree_prob_to_states(const vector<size_t> &subtree_sizes,
                    const vector<size_t> &parent_ids,
                    const vector<vector<double> > &tree_prob_table,
                    const double cutoff, vector<string> &states) {

  const size_t n_nodes = tree_prob_table[0].size();
  const size_t n_sites= tree_prob_table.size();

  for (size_t i = 0; i < n_sites; ++i) {
    string state;
    for (size_t j = 0; j < n_nodes; ++j)
      state += (tree_prob_table[i][j] <= cutoff ? '1' : '0');
    states.push_back(state);
  }

  for (size_t i = 0; i < n_sites; ++i)
    local_parsimony(subtree_sizes, parent_ids, 0, states[i]);
}



template <class T> void
write_treeprob_states(const vector<size_t> &subtree_sizes,
                      const vector<MSite> &sites,
                      const vector<vector<T> > &states,
                      const vector<string> &node_names,
                      const string &outfile) {

  const string post_outfile(outfile + ".posteriors");
  std::ofstream out(post_outfile.c_str());
  if (!out)
    throw std::runtime_error("bad output file: " + post_outfile);

  copy(node_names.begin(), node_names.end(),
       ostream_iterator<string>(out, "\t"));
  out << endl;

  for (size_t i = 0; i < states.size(); ++i) {
    out << sites[i].chrom << '\t' << sites[i].pos;
    for (size_t j = 0; j < states[i].size(); ++j)
      out << '\t' << states[i][j];
    out << endl;
  }
}


static void
build_domain(const size_t minCpG, const size_t desert_size,
             const vector<MSite> &sites, const vector<string> &states,
             vector<GenomicRegion> &domains) {

  // first round collapsing
  GenomicRegion cpg(sites[0].chrom, sites[0].pos, sites[0].pos + 1,
                    states[0], 1.0, '+');
  for (size_t i = 1; i < sites.size(); ++i) {
    const size_t d = distance(sites[i], sites[i - 1]);
    if (d < desert_size && states[i] == states[i - 1]) {
      cpg.set_score(1.0 + cpg.get_score());
      cpg.set_end(sites[i].pos + 1);
    }
    else {
      domains.push_back(cpg);
      cpg.set_chrom(sites[i].chrom);
      cpg.set_start(sites[i].pos);
      cpg.set_end(sites[i].pos + 1);
      cpg.set_name(states[i]);
      cpg.set_score(1.0);
    }
  }
  domains.push_back(cpg);

  // Iteratively merge domains
  for (size_t i = 1; i <= minCpG; ++i) {
    vector<GenomicRegion> merging_domains;
    size_t skip = 0;
    for (size_t j = 0; j < domains.size(); ++j) {
      if (domains[j].get_score() <= i) {
        skip += domains[j].get_score();
      } else {
        if (merging_domains.size() > 0 &&
            domains[j].get_name() == merging_domains.back().get_name() &&
            domains[j].distance(merging_domains.back()) < desert_size) {
          merging_domains.back().set_score(merging_domains.back().get_score() +
                                           domains[j].get_score() + skip);
          merging_domains.back().set_end(domains[j].get_end());
          skip = 0;
        } else {
          merging_domains.push_back(domains[j]);
          skip = 0;
        }
      }
    }
    domains.swap(merging_domains);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


// ADS: not sure if this function is ever really needed???
template <class T> static void
sample_initial_states(const vector<vector<double> > &tree_probs,
                      vector<vector<T> > &sampled_states) {

  uniform_real_distribution<> dis(0, 1);  // interval [0, 1)

  for (size_t i = 0; i < tree_probs.size(); ++i)
    for (size_t j = 0; j < tree_probs[i].size(); ++j)
      sampled_states[i][j] =
        (dis(gen) > (missing_meth_value(tree_probs[i][j]) ? 0.5 :
                     tree_probs[i][j])) ? 0 : 1;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////          MAIN             ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    static const size_t max_app_iter = 100;        //MAGIC
    static const size_t mh_max_iterations = 500;   //MAGIC

    size_t desert_size = 1000;
    size_t minCpG = 10;                            //MAGIC
    size_t min_sites_per_frag = 5;  // ADS: probably should be a post-filter step
    size_t MAXITER = 100;           // iterations inside the M-step
    size_t EMMAXITER = 30;          //rounds of EM iterations

    string outfile;
    string paramfile, outparamfile;

    // run mode flags
    bool VERBOSE = false;
    bool report_posteriors = false;
    bool COMPLETE = false;

    OptionParser opt_parse(strip_path(argv[0]), "Estimate phylogeny shape "
                           "and methylation state transition rates for "
                           "methylome evolution",
                           "<newick> <hypoprob-tab>");
    opt_parse.add_opt("minCpG", 'm', "minimum observed #CpGs in a block "
                      "(default: " + to_string(minCpG) + ")", false, minCpG);
    opt_parse.add_opt("maxiter", 'i', "max number of iterations "
                      "(default: " + to_string(EMMAXITER) + ")", false, EMMAXITER);
    opt_parse.add_opt("complete", 'c', "complete observations",
                      false, COMPLETE);
    opt_parse.add_opt("verbose", 'v', "print more run info "
                      "(default: " + to_string(VERBOSE) + ")", false, VERBOSE);
    opt_parse.add_opt("params", 'p', "given parameters", false, paramfile);
    opt_parse.add_opt("outparams", 'P', "output parameters", false, outparamfile);
    opt_parse.add_opt("output", 'o', "output file name", true, outfile);
    opt_parse.add_opt("min-fragsize", 'f', "min CpG in fragments to output "
                      "(default: " + toa(min_sites_per_frag) + ")", false,
                      min_sites_per_frag);
    opt_parse.add_opt("sites", 's', "write posteriors for each site"
                      "(when -o is used)", false, report_posteriors);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string tree_file = leftover_args.front();
    const string meth_table_file = leftover_args.back();
    /********************* END COMMAND LINE OPTIONS ***************************/

    /******************** LOAD PHYLOGENETIC TREE ******************************/
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in)
      throw SMITHLABException("cannot read: " + tree_file);
    // ADS: is there a need to check if the tree is consistent with
    // the one given in the params file (assuming params file
    // specified)?
    PhyloTreePreorder t;
    tree_in >> t;

    vector<size_t> subtree_sizes, node_degrees;
    vector<double> branches;
    t.get_subtree_sizes(subtree_sizes);
    t.get_branch_lengths(branches);
    get_degrees(subtree_sizes, node_degrees);

    vector<string> node_names;
    t.assign_missing_node_names();
    t.get_node_names(node_names);

    if (!is_semi_binary(node_degrees))
      throw SMITHLABException("invalid tree structure");
    const size_t n_leaves = count_leaves(subtree_sizes);

    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);

    if (VERBOSE)
      cerr << "loaded tree (leaves=" << n_leaves << ")" << endl;

    /******************** INITIALIZE PARAMETERS *******************************/
    double pi0 = 0.5; // MAGIC
    double rate0 = 0.5;// MAGIC
    double g0 = 0.9;   // MAGIC
    double g1 = 0.9;  // MAGIC

    param_set init_ps(pi0, rate0, g0, g1);

    if (!paramfile.empty()) {
      init_ps.read(paramfile, t);
      if (VERBOSE)
        cerr << "[tree:]\n" << t.tostring() << endl
             << "[given params:]\n" << init_ps << endl;
    }
    else {
      // if optimizing params, tree still needed from its own file
      vector<double> branches;
      t.get_branch_lengths(branches);
      // ADS: transformation below should be done in the same place as
      // if tree is read from the param file
      for (size_t i = 0; i < branches.size(); ++i)
        init_ps.T.push_back(1.0 - 1.0/exp(branches[i]));
      if (VERBOSE)
        cerr << "[tree:]\n" << t.tostring() << endl
             << "[starting_params={" << init_ps << "}]" << endl;
    }

    const size_t n_nodes = subtree_sizes.size();

    if (COMPLETE) {
      if (VERBOSE)
        cerr << "============ COMPLETE MODE =============" << endl
             << "[Loading full table]" << endl;

      vector<MSite> sites;
      vector<string> meth_table_file_species;
      vector<vector<double> > tree_prob_table;
      read_meth_table(meth_table_file, sites, meth_table_file_species,
                      tree_prob_table);

      if (meth_table_file_species.size() != n_nodes) {
        cerr << "complete data specified but inconsistent tree sizes:" << endl
             << meth_table_file << endl
             << tree_file << endl;
        return EXIT_SUCCESS;
      }

      if (VERBOSE)
        cerr << "[separating deserts]" << endl;
      vector<pair<size_t, size_t> > reset_points;
      separate_regions(desert_size, sites, reset_points);
      const size_t n_sites = tree_prob_table.size();

      // true counts
      pair<double, double> root_start_counts;
      pair_state root_counts;
      vector<pair_state> start_counts;
      vector<triple_state> triad_counts;

      vector<vector<bool> > tree_states(n_sites, vector<bool>(n_nodes, false));
      sample_initial_states(tree_prob_table, tree_states);
      count_triads(subtree_sizes, parent_ids, tree_states, reset_points,
                   root_start_counts, root_counts, start_counts, triad_counts);


      if (VERBOSE) {
        cerr << "-----triad_counts_true------" << endl;
        for (size_t i = 0; i < triad_counts.size(); ++i)
          cerr << triad_counts[i];
      }

      // One M-step: optimize parameters
      /*M-step: optimize parameters */
      maximization_step(VERBOSE, MAXITER, subtree_sizes,
                        root_start_counts, root_counts, start_counts,
                        triad_counts, init_ps, t);
    }
    else { /* not complete data */
      if (VERBOSE)
        cerr << endl << "========= LEAF (INCOMPLETE DATA) MODE =========" << endl
             << "[Loading LEAF table]" << endl;

      vector<MSite> sites;
      vector<vector<double> > leaf_prob_table; // can have missing leaf values
      vector<string> meth_table_species;
      read_meth_table(meth_table_file, sites, meth_table_species, leaf_prob_table);

      if (VERBOSE)
        cerr << "loaded meth table (species="
             << meth_table_species.size() << ")" << endl;

      // make sure meth data and tree info is in sync
      if (meth_table_species.size() != n_leaves ||
          !has_same_species_order(t, meth_table_species))
        throw SMITHLABException("inconsistent species counts, names or order");

      if (VERBOSE) {
        vector<size_t> species_in_order;
        subtree_sizes_to_leaves_preorder(subtree_sizes, species_in_order);
        for (size_t i = 0; i < species_in_order.size(); ++i)
          cerr << meth_table_species[i] << "\t" << species_in_order[i] << endl;
        cerr << "[total_sites=" << sites.size() << "]" << endl;
      }

      if (VERBOSE)
        cerr << "[separating deserts]" << endl;
      vector<pair<size_t, size_t> > reset_points;
      separate_regions(desert_size, sites, reset_points);
      const size_t n_sites = sites.size();

      // initialize parameters based on leaf values
      double g0_est = 0.0, g1_est = 0.0;
      estimate_g0_g1(leaf_prob_table, g0_est, g1_est);
      const double pi0_est = estimate_pi0(leaf_prob_table);

      if (VERBOSE)
        cerr << "[post-filter=" << n_sites << "; blocks="
             << reset_points.size() << "]" << endl;

      /* heuristic estimates for a few parameters*/
      if (paramfile.empty()) {
        init_ps.pi0 = pi0_est;
        init_ps.g0 = g0_est;
        init_ps.g1 = g1_est;
      }

      /* tree_prob_table_full has missing data filled */
      vector<vector<double> > tree_probs(sites.size(),
                                         vector<double>(n_nodes, -1.0));
      // ADS: this should probably be done as soon as the leaf states
      // are read.
      leaf_to_tree_prob(subtree_sizes, leaf_prob_table, tree_probs);

      // summary statistics
      pair<double, double> root_start_counts;
      pair_state root_counts;
      vector<pair_state> start_counts;
      vector<triple_state> triad_counts;

      bool em_converged = !paramfile.empty(); // params fixed => converged

      // Initialize "sample" states from the heuristic probabilities
      vector<vector<bool> > tree_states(n_sites, vector<bool>(n_nodes));
      sample_initial_states(tree_probs, tree_states);

      param_set curr_ps(init_ps);
      for (size_t iter = 0; iter < EMMAXITER && !em_converged; ++iter) {

        if (VERBOSE)
          cerr << endl << "====================[EM ITERATION=" << iter
               << "]=======================" << endl;

        const param_set prev_ps(curr_ps);
        // next line is optional, just for consistency with older version.
        expectation_step(VERBOSE, mh_max_iterations, max_app_iter,
                         subtree_sizes, parent_ids, reset_points, curr_ps,
                         root_start_counts, root_counts,
                         start_counts, triad_counts,
                         tree_probs, tree_states);

        if (VERBOSE) {
          cerr << "[E-step iter=" << iter << "]\ntriad counts:\n" << endl;
          for (size_t i = 0; i < triad_counts.size(); ++i)
            cerr << triad_counts[i] << "\tnode=" << i << endl;
        }

        /****************** M-step: optimize parameters ************************/
        maximization_step(VERBOSE, MAXITER, subtree_sizes,
                          root_start_counts, root_counts, start_counts,
                          triad_counts, curr_ps, t);
        if (VERBOSE)
          cerr << "[M-step iter=" << iter << ", params:" << endl
               << curr_ps << ']' << endl;

        const double diff = param_set::absolute_difference(prev_ps, curr_ps);

        const double llk = log_likelihood(subtree_sizes, curr_ps, root_start_counts,
                                          root_counts, start_counts, triad_counts);

        em_converged = (diff < param_set::tolerance*curr_ps.T.size());

        if (VERBOSE)
          cerr << "[EM iter=" << iter << ", "
               << "log_lik=" << llk << ", "
               << "delta=" << diff << ", "
               << "conv=" << em_converged << ']' << endl;
      } // end E-M iterations

      if (VERBOSE)
        cerr << "[computing posterior probabilities]" << endl;
      vector<vector<double> > posteriors;
      mcmc_estimate_posterior(VERBOSE, mh_max_iterations, subtree_sizes,
                              parent_ids, reset_points, curr_ps,
                              tree_states, tree_probs, posteriors);

      init_ps = curr_ps; // final parameters

      // update phylotree with branch lenghts
      init_ps.assign_branches(t);

      /************************** Output states ********************************/
      std::ofstream out(outfile.c_str());
      if (!out) // ADS: crazy to check at end...
        throw SMITHLABException("bad output file: " + outfile);

      /*build domain*/
      vector<string> states;
      vector<GenomicRegion> domains;
      double cutoff = 0.5; // MAGIC
      tree_prob_to_states(subtree_sizes, parent_ids, tree_probs, cutoff, states);
      build_domain(min_sites_per_frag, desert_size, sites, states, domains);
      if (VERBOSE)
        cerr << "Built total " << domains.size() << " domains" << endl;

      copy(domains.begin(), domains.end(),
           ostream_iterator<GenomicRegion>(out, "\n"));

      /* output individual sites */
      if (report_posteriors)
        write_treeprob_states(subtree_sizes, sites, posteriors,
                              node_names, outfile);
    }

    if (!outparamfile.empty())
      init_ps.write(t, outparamfile);
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/*
  TODO:

  (1) Make sure that the pi0 is being estimated correctly in the
  maximum likelihood calculations. Right now it is only using the
  first site of each block, and at the root.

  (2) Eliminate the optimize_params function from the
  optimize_params.*pp. This function is actually completely arbitrary,
  and shouldn't be used, due to the arbitrary order of estimates of
  parameters. The contents of that function should be in the
  maximization_step function.

  (3) [DONE] Ensure that the transformation between branch lengths
  is entirely taken care of inside the param_set class. This should
  never be seen by readers of this source file. At the same time,
  likely it should be printed with some notation to remind the user
  that the verbose output prints transformed values.

  (4) Determine whether any states should not have their values
  printed in the final output, due to their distance from any actual
  leaf data making them meaningless.

  (5) Determine whether any states should not be involved in parameter
  estimation, as they are so far from any observed data that they just
  waste time.

  (6) Decide whether we need to have the ternary leaf, or if there is
  enough information to estimate both branches out of it. This would
  mean we have more information than is assumed in the papers that
  show non-identifiability, but I wouldn't be too surprised.

  (7) I'm stull concerned about convergence of the mcmc, and we can
  implement multiple chains without taking too much extra space. We
  can encode many boolean values inside a machine word, and have them
  updated concurrently.

*/
