/*    Copyright (C) 2015-16 University of Southern California and
 *                          Andrew D. Smith and Jenny Qu
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
#include <numeric>
#include <random>
#include <algorithm>  //std::max, min
#include <cmath>      //std::abs
#include <limits>     //std::numeric_limits
#include <iterator>   //std::distance

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

/* from methpipe */
#include "MethpipeSite.hpp"

/* headers for epigenomic evolution */
#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"
#include "param_set.hpp"
#include "epiphy_utils.hpp"
#include "sufficient_statistics_helpers.hpp"
#include "optimize_params.hpp"
#include "epiphy_mcmc.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::pair;
using std::make_pair;
using std::numeric_limits;
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
static const double KL_CONV_TOL = 1e-8;    //MAGIC

#include <random>
using std::uniform_real_distribution;
//std::random_device rd; //seed generator
//std::mt19937_64 gen(rd()); //generator initialized with seed from rd
std::mt19937_64 gen(0); //generator initialized with seed 0


static void
separate_regions(const size_t desert_size, vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &blocks) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      blocks.push_back(std::make_pair(i, i));
    else blocks.back().second = i;
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
estimate_g0_g1(const vector<vector<double> > &meth, double &g0, double &g1) {

  static const double pseudo_count = 1.0;

  double g00 = pseudo_count, g01 = pseudo_count;
  double g10 = pseudo_count, g11 = pseudo_count;
  for (size_t i = 0; i < meth.size() - 1; ++i)
    for (size_t j = 0; j < meth[i].size(); ++j) {
      const double curr = meth[i][j], next = meth[i + 1][j];
      if (!missing_meth_value(curr) && !missing_meth_value(next)) {
        g00 += (1.0 - curr)*(1.0 - next);
        g01 += (1.0 - curr)*next;
        g10 += curr*(1.0 - next);
        g11 += curr*next;
      }
    }

  g0 = g00/(g01 + g00);
  g1 = g11/(g11 + g10);
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
    // ADS: this is not the right transformation to put here; we do
    // not want to add to these probabilities without ensuring they
    // sum to not more than 1.0
    transform(p.begin(), p.end(), p.begin(),
              bind(plus<double>(), _1, PROBABILITY_GUARD));
    transform(q.begin(), q.end(), q.begin(),
              bind(plus<double>(), _1, PROBABILITY_GUARD));
    kld[i] = kl_divergence(p, q);
  }
}


static void
add_internal_node_probs(const vector<size_t> &subtree_sizes,
                        vector<vector<double> > &prob_table) {

  const size_t n_nodes = subtree_sizes.size();
  const size_t n_sites = prob_table.size();

  vector<vector<double> > expanded_probs(n_sites, vector<double>(n_nodes, -1.0));

  vector<size_t> leaves_preorder;
  subtree_sizes_to_leaves_preorder(subtree_sizes, leaves_preorder);

  for (size_t i = 0; i < n_sites; ++i)
    for (size_t j = 0; j < leaves_preorder.size(); ++j)
      // ADS: should these be checked for "extremes"???
      // JQU: checking for extremes should be done when reading in prob_table.
      expanded_probs[i][leaves_preorder[j]] = prob_table[i][j];

  prob_table.swap(expanded_probs);
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
                  const vector<triple_state> &triad_counts, param_set &ps) {

  // one M-step: optimize parameters
  double diff = std::numeric_limits<double>::max();
  for (size_t iter = 0; iter < MAXITER && diff > param_set::tolerance; ++iter) {
    if (VERBOSE)
      cerr << "[inside maximization: iter=" << iter << "]" << endl
           << ps << endl;

    param_set prev_ps(ps);
    for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id)
      max_likelihood_branch(VERBOSE, subtree_sizes, node_id,
                            start_counts, triad_counts, ps);

    max_likelihood_rate(VERBOSE, subtree_sizes, start_counts, triad_counts, ps);

    max_likelihood_pi0(VERBOSE, root_start_counts, ps);

    max_likelihood_horiz(VERBOSE, subtree_sizes, root_counts, triad_counts, ps);

    diff = param_set::absolute_difference(ps, prev_ps); // convergence criteria
  }
}


template <class T>
static void
expectation_step(const bool VERBOSE,
                 const size_t mh_max_iterations,
                 const size_t burnin,
                 const epiphy_mcmc &sampler,
                 const vector<size_t> &subtree_sizes,
                 const vector<size_t> &parent_ids,
                 const vector<vector<double> > &tree_probs,
                 const vector<pair<size_t, size_t> > &blocks,
                 const param_set &ps,
                 pair<double, double> &root_start_counts,
                 pair_state &root_counts,
                 vector<pair_state> &start_counts,
                 vector<triple_state> &triad_counts,
                 vector<vector<T> > &sampled_states) {

  const size_t n_nodes = subtree_sizes.size();

  root_start_counts = make_pair(0.0, 0.0);
  triad_counts = vector<triple_state>(n_nodes);
  start_counts = vector<pair_state>(n_nodes);
  root_counts = pair_state();

  bool converged = false;
  size_t mh_iter = 0;
  size_t burned = 0;
  for (mh_iter = 0; mh_iter < mh_max_iterations && !converged; ++mh_iter) {
    if (VERBOSE)
      cerr << "\r[inside expectation: M-H (iter=" << mh_iter << ")]";

    // take the sample
    sampler.sample_states(subtree_sizes, parent_ids, ps, tree_probs,
                          blocks, sampled_states);

    pair<double, double> root_start_counts_samp;
    pair_state root_counts_samp;
    vector<pair_state> start_counts_samp;
    vector<triple_state> triad_counts_samp;
    count_triads(subtree_sizes, parent_ids, sampled_states, blocks,
                 root_start_counts_samp, root_counts_samp, start_counts_samp,
                 triad_counts_samp);

    // retain previous values for comparison with next iteration
    vector<triple_state> triad_counts_prev(triad_counts);

    if (burned < burnin) {
      ++burned;
      if (VERBOSE)
        cerr << "\r[inside expectation: M-H (iter=" << mh_iter << ") burned]";
    } else {
      // update the counts with current iteration samples
      root_start_counts.first += root_start_counts_samp.first;
      root_start_counts.second += root_start_counts_samp.second;
      root_counts += root_counts_samp;
      for (size_t i = 0; i < triad_counts.size(); ++i) {
        triad_counts[i] += triad_counts_samp[i];
        start_counts[i] += start_counts_samp[i];
      }

      // check how far the proportions have moved; only the triad_counts
      // used for convergence testing; the start and root counts do not
      // contribute enough to check
      vector<double> divergence;
      kl_divergence(triad_counts_prev, triad_counts, divergence);

      // determine convergence based on how far proportions have moved
      converged =
        *max_element(divergence.begin(), divergence.end()) < KL_CONV_TOL;
    }
  }
  if (VERBOSE)
    cerr << endl;

  root_start_counts.first /= mh_iter - burnin;
  root_start_counts.second /= mh_iter - burnin;
  root_counts.div(mh_iter);
  for (size_t i = 0; i < triad_counts.size(); ++i) {
    start_counts[i].div(mh_iter - burnin);
    triad_counts[i].div(mh_iter - burnin);
  }

  if (VERBOSE)
    // ADS: should we print more summary statistics here?
    cerr << "[MH iterations=" << mh_iter << ']' << endl;
}


template <class T>
static void
expectation_maximization(const bool VERBOSE,
                         const size_t em_max_iterations,
                         const size_t opt_max_iterations,
                         const size_t mh_max_iterations,
                         const epiphy_mcmc &sampler,
                         const vector<size_t> &subtree_sizes,
                         const vector<size_t> &parent_ids,
                         const vector<vector<double> > &tree_probs,
                         const vector<pair<size_t, size_t> > &blocks,
                         param_set &params,
                         pair<double, double> &root_start_counts,
                         pair_state &root_counts,
                         vector<pair_state> &start_counts,
                         vector<triple_state> &triad_counts,
                         vector<vector<T> > &sampled_states) {

  bool em_converged = false;
  for (size_t iter = 0; iter < em_max_iterations && !em_converged; ++iter) {

    if (VERBOSE)
      cerr << endl << "====================[EM ITERATION=" << iter
           << "]=======================" << endl;

    const param_set prev_ps(params);
    // next line is optional, just for consistency with older version.
    const size_t burnin = 20;
    expectation_step(VERBOSE, mh_max_iterations, burnin, sampler, subtree_sizes,
                     parent_ids, tree_probs, blocks, params, root_start_counts,
                     root_counts, start_counts, triad_counts, sampled_states);

    if (VERBOSE) {
      cerr << "[E-step iter=" << iter << "]\ntriad counts:\n" << endl;
      for (size_t i = 0; i < triad_counts.size(); ++i)
        cerr << triad_counts[i] << "\tnode=" << i << endl;
    }

    /****************** M-step: optimize parameters ************************/
    maximization_step(VERBOSE, opt_max_iterations, subtree_sizes,
                      root_start_counts, root_counts, start_counts,
                      triad_counts, params);
    if (VERBOSE)
      cerr << "[M-step iter=" << iter << ", params:" << endl
           << params << ']' << endl;

    const double diff = param_set::absolute_difference(prev_ps, params);

    const double llk = log_likelihood(subtree_sizes, params, root_start_counts,
                                      root_counts, start_counts, triad_counts);

    em_converged = (diff < param_set::tolerance*params.T.size());

    if (VERBOSE)
      cerr << "[EM iter=" << iter << ", "
           << "log_lik=" << llk << ", "
           << "delta=" << diff << ", "
           << "conv=" << em_converged << ']' << endl;
  }
}


template <class T>
static void
estimate_posterior(const bool VERBOSE,
                   const size_t max_iterations,
                   const epiphy_mcmc &sampler,
                   const vector<size_t> &subtree_sizes,
                   const vector<size_t> &parent_ids,
                   const vector<vector<double> > &marginals,
                   const vector<pair<size_t, size_t> > &blocks,
                   const param_set &params,
                   vector<vector<T> > &states,
                   vector<vector<double> > &posteriors) {

  const size_t n_nodes = subtree_sizes.size();

  posteriors = vector<vector<double> >(marginals.size(),
                                       vector<double>(n_nodes, 0.0));

  // ADS: need convergence criteria here??
  for (size_t mh_iter = 0; mh_iter < max_iterations; ++mh_iter) {
    if (VERBOSE)
      cerr << "\r[inside estimate_posterior (iter=" << mh_iter << ")]";

    sampler.sample_states(subtree_sizes, parent_ids, params,
                          marginals, blocks, states);

    for (size_t i = 0; i < states.size(); ++i)
      for (size_t j = 0; j < states[i].size(); ++j)
        posteriors[i][j] += states[i][j];
  }
  if (VERBOSE)
    cerr << endl;

  for (size_t i = 0; i < posteriors.size(); ++i)
    for (size_t j = 0; j < posteriors[i].size(); ++j)
      posteriors[i][j] /= max_iterations;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////// FORMATTING AND WRITING OUTPUT /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


static size_t
binary_parsimony_forward(const vector<size_t> &subtree_sizes,
                         const size_t node_id, string &states) {

  size_t pars_score = 0;
  if (!is_leaf(subtree_sizes[node_id])) {

    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      pars_score +=
        binary_parsimony_forward(subtree_sizes, child_id, states);
      count += subtree_sizes[child_id];
    }

    // check intersection
    char s = '-';
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      if (s != states[child_id]) // do nothing if no conflict
        s = (s == '-') ? states[child_id] : 0;
      count += subtree_sizes[child_id];
    }
    states[node_id] = (s == 0) ? '-' : s;
    pars_score += (s == 0);
  }
  return pars_score;
}

static void
binary_parsimony_traceback(const vector<size_t> &subtree_sizes,
                           const size_t node_id, string &states) {

  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    if (states[node_id] == states[child_id] || states[child_id] == '-')
      states[child_id] = states[node_id];
    binary_parsimony_traceback(subtree_sizes, child_id, states);
    count += subtree_sizes[child_id];
  }
}


// static void
// local_parsimony(const vector<size_t> &subtree_sizes,
//                 const vector<size_t> &parent_ids,
//                 const size_t node_id, string &state) {

//   if (!is_leaf(subtree_sizes[node_id])) {

//     for (size_t count = 1; count < subtree_sizes[node_id];) {
//       const size_t child_id = node_id + count;
//       local_parsimony(subtree_sizes, parent_ids, child_id, state);
//       count += subtree_sizes[child_id];
//     }

//     char s = state[node_id + 1];
//     bool SAME = true;
//     for (size_t count = 1; count < subtree_sizes[node_id] && SAME;) {
//       const size_t child_id = node_id + count;
//       if (state[child_id] != s)
//         SAME = false;
//       count += subtree_sizes[child_id];
//     }
//     if (SAME) state[node_id] = s;
//     else state[node_id] = state[parent_ids[node_id]];
//   }
// }


static size_t
tree_prob_to_states(const vector<size_t> &subtree_sizes,
                    const vector<vector<double> > &tree_probs,
                    vector<string> &states) {
  static const double state_cutoff = 0.5; // MAGIC
  uniform_real_distribution<> dis(0, 1);  // interval [0, 1)

  states.clear();
  size_t pars_score = 0;
  for (size_t i = 0; i < tree_probs.size(); ++i) {

    string state;
    for (size_t j = 0; j < tree_probs[i].size(); ++j)
      state += (tree_probs[i][j] < state_cutoff ? '0' : '1');

    pars_score += binary_parsimony_forward(subtree_sizes, 0, state);

    if (state[0] == '-')
      state[0] = dis(gen) > 0.5 ? '0' : '1';

    binary_parsimony_traceback(subtree_sizes, 0, state);
    states.push_back(state);
  }
  return pars_score;
}


template <class T> void
write_treeprob_states(const vector<size_t> &subtree_sizes,
                      const vector<MSite> &sites,
                      const vector<vector<T> > &states,
                      const vector<string> &node_names,
                      const string &outfile) {

  std::ofstream out(outfile.c_str());
  if (!out)
    throw std::runtime_error("bad output file: " + outfile);

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
  domains.push_back(GenomicRegion(sites[0].chrom, sites[0].pos,
                                  sites[0].pos + 1, states[0], 1.0, '+'));
  for (size_t i = 1; i < sites.size(); ++i) {
    const size_t d = distance(sites[i], sites[i - 1]);
    if (d < desert_size && states[i] == states[i - 1]) {
      domains.back().set_score(1.0 + domains.back().get_score());
      domains.back().set_end(sites[i].pos + 1);
    }
    else domains.push_back(GenomicRegion(sites[i].chrom, sites[i].pos,
                                         sites[i].pos + 1, states[i], 1.0, '+'));
  }

  // Iteratively merge domains
  for (size_t i = 1; i <= minCpG; ++i) {
    vector<GenomicRegion> merging;
    size_t skip = 0;
    for (size_t j = 0; j < domains.size(); ++j) {
      if (domains[j].get_score() <= i)
        skip += domains[j].get_score();
      else {
        if (merging.size() > 0 &&
            domains[j].get_name() == merging.back().get_name() &&
            domains[j].distance(merging.back()) < desert_size) {
          merging.back().set_score(merging.back().get_score() +
                                   domains[j].get_score() + skip);
          merging.back().set_end(domains[j].get_end());
        }
        else merging.push_back(domains[j]);
        skip = 0;
      }
    }
    domains.swap(merging);
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template <class T> static void
sample_initial_states(const vector<vector<double> > &tree_probs,
                      vector<vector<T> > &sampled_states) {

  uniform_real_distribution<> dis(0, 1);

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

    static const string param_file_suffix = "epiphy_param";
    static const string posteriors_file_suffix = "epiphy_post";
    static const string intervals_file_suffix = "epiphy_intervals";

    string outprefix;
    string paramfile;

    size_t desert_size = 1000;
    size_t min_sites_per_frag = 5;  // ADS: should this be a
                                    // post-processing filter step in
                                    // a separate program?
    size_t minCpG = 10;                     // MAGIC

    size_t opt_max_iterations = 100;        // iterations inside the M-step
    size_t em_max_iterations = 100;         // rounds of EM iterations
    size_t mh_max_iterations = 500;         // MAGIC

    // run mode flags
    bool VERBOSE = false;
    bool assume_complete_data = false;

    /********************* COMMAND LINE OPTIONS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), "Estimate phylogeny shape "
                           "and methylation state transition rates for "
                           "methylome evolution",
                           "<newick> <hypoprob-tab>");
    opt_parse.add_opt("minCpG", 'm', "minimum observed #CpGs in a block "
                      "(default: " + to_string(minCpG) + ")", false, minCpG);
    opt_parse.add_opt("maxiter", 'i', "max EM iterations (default: " +
                      to_string(em_max_iterations) + ")",
                      false, em_max_iterations);
    opt_parse.add_opt("mcmc-iter", 'h', "max mcmc iterations (default: " +
                      to_string(mh_max_iterations) + ")",
                      false, mh_max_iterations);
    opt_parse.add_opt("complete", 'c', "input is complete observations",
                      false, assume_complete_data);
    opt_parse.add_opt("verbose", 'v', "print more run info "
                      "(default: " + string(VERBOSE ? "true" : "false") + ")",
                      false, VERBOSE);
    opt_parse.add_opt("params", 'p', "given parameters", false, paramfile);
    opt_parse.add_opt("outprefix", 'o', "output file prefix", true, outprefix);
    opt_parse.add_opt("min-fragsize", 'f', "min CpG in fragments to output "
                      "(default: " + toa(min_sites_per_frag) + ")", false,
                      min_sites_per_frag);

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
    /******************** END COMMAND LINE OPTIONS ********************/

    /******************** LOAD PHYLOGENETIC TREE ******************************/
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in)
      throw SMITHLABException("cannot read: " + tree_file);
    // ADS: is there a need to check if the tree is consistent with
    // the one given in the params file (assuming params file
    // specified)?
    PhyloTreePreorder t;
    tree_in >> t;

    vector<size_t> subtree_sizes;
    t.get_subtree_sizes(subtree_sizes);
    vector<double> branches;
    t.get_branch_lengths(branches);

    t.assign_missing_node_names();
    vector<string> node_names;
    t.get_node_names(node_names);

    const size_t n_nodes = subtree_sizes.size();
    const size_t n_leaves = count_leaves(subtree_sizes);

    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);

    if (VERBOSE)
      cerr << "[tree:]\n" << t.tostring() << endl;

    /******************** INITIALIZE PARAMETERS *******************************/
    param_set initial_params;
    if (!paramfile.empty()) {
      initial_params.read(paramfile, t);
      if (VERBOSE)
        cerr << "[given params={" << initial_params << "}]" << endl;
    }
    else {
      // ADS: must we assume branch lengths are provided here?
      vector<double> branches;
      t.get_branch_lengths(branches);
      for (size_t i = 0; i < branches.size(); ++i)
        branches[i] = 1.0 - 1.0/exp(branches[i]);
      const double pi0 = 0.5;   // MAGIC
      const double rate0 = 0.5; // MAGIC
      const double g0 = 0.9;    // MAGIC
      const double g1 = 0.9;    // MAGIC
      initial_params = param_set(pi0, rate0, g0, g1, branches);
      if (VERBOSE)
        cerr << "[starting params={" << initial_params << "}]" << endl;
    }

    /******************* READ THE METHYLATION DATA *****************************/
    if (VERBOSE)
      cerr << "[reading methylation data (mode="
           << (assume_complete_data ? "complete" : "missing") << ")]" << endl;

    vector<MSite> sites;
    vector<vector<double> > tree_probs;
    vector<string> meth_table_species;
    read_meth_table(meth_table_file, sites, meth_table_species, tree_probs);
    const size_t n_sites = tree_probs.size();

    if (assume_complete_data) {
      if (meth_table_species.size() != n_nodes) {
        cerr << "complete data specified but inconsistent tree sizes:" << endl
             << meth_table_file << endl
             << tree_file << endl;
        return EXIT_SUCCESS;
      }
    }
    else {
      // make sure meth data and tree info is in sync
      if (meth_table_species.size() != n_leaves ||
          !has_same_species_order(t, meth_table_species))
        throw SMITHLABException("inconsistent species counts, names or order");
      add_internal_node_probs(subtree_sizes, tree_probs);

      if (VERBOSE) {
        cerr << "number of leaf species: " << meth_table_species.size() << endl;
        vector<size_t> species_in_order;
        subtree_sizes_to_leaves_preorder(subtree_sizes, species_in_order);
        for (size_t i = 0; i < species_in_order.size(); ++i)
          cerr << meth_table_species[i] << "\t" << species_in_order[i] << endl;
        cerr << "[total_sites=" << sites.size() << "]" << endl;
      }
    }

    if (VERBOSE)
      cerr << "[separating deserts]" << endl;
    vector<pair<size_t, size_t> > blocks;
    separate_regions(desert_size, sites, blocks);
    if (VERBOSE)
      cerr << "number of blocks: " << blocks.size() << endl;

    if (paramfile.empty()) {
      // heuristic starting points for G and pi
      estimate_g0_g1(tree_probs, initial_params.g0, initial_params.g1);
      initial_params.pi0 = estimate_pi0(tree_probs);
    }

    // true counts
    pair<double, double> root_start_counts;
    pair_state root_counts;
    vector<pair_state> start_counts;
    vector<triple_state> triad_counts;

    vector<vector<bool> > tree_states(n_sites, vector<bool>(n_nodes, false));
    // for complete data, the sampling will have probability 1 or 0
    sample_initial_states(tree_probs, tree_states);

    param_set params(initial_params);

    if (assume_complete_data) {
      count_triads(subtree_sizes, parent_ids, tree_states, blocks,
                   root_start_counts, root_counts, start_counts, triad_counts);
      if (paramfile.empty())
        maximization_step(VERBOSE, opt_max_iterations, subtree_sizes,
                          root_start_counts, root_counts, start_counts,
                          triad_counts, params);
      if (VERBOSE) {
        const double llk =
          log_likelihood(subtree_sizes, params, root_start_counts,
                         root_counts, start_counts, triad_counts);
        // ADS: for complete data and params provided, this is all
        // that will be output
        cerr << "log_likelihood=" << llk << endl;
      }
    }
    else {

      const epiphy_mcmc sampler(mh_max_iterations, 0);

      if (paramfile.empty())
        expectation_maximization(VERBOSE, em_max_iterations, opt_max_iterations,
                                 mh_max_iterations, sampler, subtree_sizes, parent_ids,
                                 tree_probs, blocks, params,
                                 root_start_counts, root_counts,
                                 start_counts, triad_counts, tree_states);

      if (VERBOSE)
        cerr << "[computing posterior probabilities]" << endl;
      vector<vector<double> > posteriors;
      estimate_posterior(VERBOSE, mh_max_iterations, sampler,
                         subtree_sizes, parent_ids, tree_probs,
                         blocks, params, tree_states, posteriors);

      /************************** Output states ********************************/
      if (VERBOSE)
        cerr << "[building domains of contiguous state]" << endl;
      vector<string> states;
      vector<GenomicRegion> domains;
      tree_prob_to_states(subtree_sizes, tree_probs, states);
      build_domain(min_sites_per_frag, desert_size, sites, states, domains);
      if (VERBOSE)
        cerr << "total domains: " << domains.size() << endl;

      // output the domains as contiguous intervals
      const string interval_outfile(outprefix + "." + intervals_file_suffix);
      std::ofstream out(interval_outfile.c_str());
      if (!out)
        throw SMITHLABException("bad output file: " + interval_outfile);
      copy(domains.begin(), domains.end(),
           ostream_iterator<GenomicRegion>(out, "\n"));

      // output the posteriors on states
      const string post_outfile(outprefix + "." + posteriors_file_suffix);
      write_treeprob_states(subtree_sizes, sites, posteriors,
                            node_names, post_outfile);
    }

    if (VERBOSE) {
      cerr << "[sufficient statistics]" << endl;
      cerr << "root_start_counts:\n"
           << root_start_counts.first << '\t'
           << root_start_counts.second << endl
           << "root_counts:\n" << root_counts << endl
           << "start_counts:\n";
      copy(start_counts.begin(), start_counts.end(),
           ostream_iterator<pair_state>(cerr, "\n"));
      cerr << "triad_counts:\n";
      copy(triad_counts.begin(), triad_counts.end(),
           ostream_iterator<triple_state>(cerr, "\n"));
    }

    const string param_outfile(outprefix + "." + param_file_suffix);
    params.write(t, param_outfile);
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

  (4) Determine whether any states should not have their values
  printed in the final output, due to their distance from any actual
  leaf data making them meaningless.

  separate site filters for training and posterior estimation.

  (3) The part of this program that actually identifies contiguous
  domains should be removed and put in a separate program. The reason
  is that this is the program that we expect should take a long time
  to run, and afterwards we should be able to use all learned
  parameters and posteriors for doing any quick things -- without
  assuming this program gets rerun.

  (5) Determine whether any states should not be involved in parameter
  estimation, as they are so far from any observed data that they just
  waste time.

  (6) Decide whether we need to have the ternary root, or if there is
  enough information to estimate both branches out of it. This would
  mean we have more information than is assumed in the papers that
  show non-identifiability, but I wouldn't be too surprised.

  (7) I'm still concerned about convergence of the mcmc, and we can
  implement multiple chains without taking too much extra space. We
  can encode many boolean values inside a machine word, and have them
  updated concurrently.

*/
