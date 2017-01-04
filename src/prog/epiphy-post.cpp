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
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <random>
#include <algorithm>  //std::max, min
#include <cmath>      //std::abs
#include <limits>     //std::numeric_limits
#include <random>

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

/* from methpipe */
#include "MethpipeSite.hpp"

/* headers for epigenomic evolution */
#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"
#include "param_set.hpp"
#include "epiphy_utils.hpp"
#include "epiphy_mcmc.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::pair;
using std::numeric_limits;
using std::min;
using std::max;
using std::to_string;
using std::ostream_iterator;


static void
separate_regions(const size_t desert_size, vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &blocks) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      blocks.push_back(std::make_pair(i, i));
    else blocks.back().second = i;
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
/////////////// POSTERIOR ESTIMATION     ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



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


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

template <class T> static void
sample_initial_states(const size_t rng_seed,
                      const vector<vector<double> > &tree_probs,
                      vector<vector<T> > &sampled_states) {

  std::random_device rd; // random devide to generate seed
  std::mt19937_64 gen(rng_seed == numeric_limits<size_t>::max() ? rd() : rng_seed);

  std::uniform_real_distribution<> dis(0, 1);

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

    size_t rng_seed = numeric_limits<size_t>::max();

    string outfile;

    size_t desert_size = 1000;
    size_t mh_max_iterations = 500;         // MAGIC

    // run mode flags
    bool VERBOSE = false;
    bool assume_complete_data = false;

    /********************* COMMAND LINE OPTIONS ***********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "estimate posteriors on methylation states",
                           "<param-file> <meth-table>");
    opt_parse.add_opt("mcmc-iter", 'h', "max mcmc iterations (default: " +
                      to_string(mh_max_iterations) + ")",
                      false, mh_max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info "
                      "(default: " + string(VERBOSE ? "true" : "false") + ")",
                      false, VERBOSE);
    opt_parse.add_opt("seed", 's', "rng seed (default: none)",
                      false, rng_seed);

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
    const string paramfile = leftover_args.front();
    const string meth_table_file = leftover_args.back();
    /******************** END COMMAND LINE OPTIONS ********************/

    /******************** INITIALIZE PARAMETERS *******************************/
    PhyloTreePreorder t;
    param_set params;
    params.read(paramfile, t);
    if (VERBOSE)
      cerr << "[given params={" << params << "}]" << endl;

    /******************** LOAD PHYLOGENETIC TREE ******************************/
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

    /******************* READ THE METHYLATION DATA *****************************/
    if (VERBOSE)
      cerr << "[reading methylation data (mode="
           << (assume_complete_data ? "complete" : "missing") << ")]" << endl;

    vector<MSite> sites;
    vector<vector<double> > tree_probs;
    vector<string> meth_table_species;
    read_meth_table(meth_table_file, sites, meth_table_species, tree_probs);
    const size_t n_sites = tree_probs.size();

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

    if (VERBOSE)
      cerr << "[separating deserts]" << endl;
    vector<pair<size_t, size_t> > blocks;
    separate_regions(desert_size, sites, blocks);
    if (VERBOSE)
      cerr << "number of blocks: " << blocks.size() << endl;

    vector<vector<bool> > tree_states(n_sites, vector<bool>(n_nodes, false));
    // for complete data, the sampling will have probability 1 or 0
    sample_initial_states(rng_seed, tree_probs, tree_states);

    const epiphy_mcmc sampler(mh_max_iterations, 0);

    if (VERBOSE)
      cerr << "[computing posterior probabilities]" << endl;
    vector<vector<double> > posteriors;
    estimate_posterior(VERBOSE, mh_max_iterations, sampler,
                       subtree_sizes, parent_ids, tree_probs,
                       blocks, params, tree_states, posteriors);

    /************************** Output states ********************************/

    // output the posteriors on states
    write_treeprob_states(subtree_sizes, sites, posteriors,
                          node_names, outfile);
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