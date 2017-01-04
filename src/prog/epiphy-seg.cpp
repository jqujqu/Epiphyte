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

#include <random>
using std::uniform_real_distribution;
//std::random_device rd; //seed generator
//std::mt19937_64 gen(rd()); //generator initialized with seed from rd
std::mt19937_64 gen(0); //generator initialized with seed 0

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


static void
separate_regions(const size_t desert_size, vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &blocks) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      blocks.push_back(std::make_pair(i, i));
    else blocks.back().second = i;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////// FORMATTING AND WRITING OUTPUT /////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


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
//////////////////////          MAIN             ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    string outfile;

    size_t desert_size = 1000;
    size_t min_sites_per_frag = 5;  // ADS: should this be a
                                    // post-processing filter step in
                                    // a separate program?
    // size_t minCpG = 10;                     // MAGIC

    // run mode flags
    bool VERBOSE = false;

    /********************* COMMAND LINE OPTIONS ***********************/
    OptionParser opt_parse(strip_path(argv[0]), "Estimate phylogeny shape "
                           "and methylation state transition rates for "
                           "methylome evolution",
                           "<newick> <hypoprob-tab>");
    opt_parse.add_opt("verbose", 'v', "print more run info "
                      "(default: " + string(VERBOSE ? "true" : "false") + ")",
                      false, VERBOSE);
    opt_parse.add_opt("out", 'o', "output file name (default: stdout)",
                      true, outfile);
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
      cerr << "[reading methylation data (states or posteriors)]" << endl;

    vector<MSite> sites;
    vector<vector<double> > tree_probs;
    vector<string> meth_table_species;
    read_meth_table(meth_table_file, sites, meth_table_species, tree_probs);
    const size_t n_sites = tree_probs.size();

    if (meth_table_species.size() != n_nodes) {
      cerr << "complete data specified but inconsistent tree sizes:" << endl
           << meth_table_file << endl
           << paramfile << endl;
      return EXIT_SUCCESS;
    }

    if (VERBOSE)
      cerr << "[separating deserts]" << endl;
    vector<pair<size_t, size_t> > blocks;
    separate_regions(desert_size, sites, blocks);
    if (VERBOSE)
      cerr << "number of blocks: " << blocks.size() << endl;

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
    std::ofstream out(outfile.c_str());
    if (!out)
      throw SMITHLABException("bad output file: " + outfile);
    copy(domains.begin(), domains.end(),
         ostream_iterator<GenomicRegion>(out, "\n"));

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
