/*****************************************************************************
 *  epiphy-test:
 *
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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <random>
#include <iomanip>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "MethpipeSite.hpp"

#include "PhyloTreePreorder.hpp"
#include "param_set.hpp"
#include "sufficient_statistics_helpers.hpp"
#include "optimize_params.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::pair;

using std::max;
using std::min;


static bool
parse_line(const string &line,
           vector<MSite> &sites, vector<vector<bool> > &states) {

  std::istringstream iss(line);

  sites.push_back(MSite());
  if (!(iss >> sites.back().chrom >> sites.back().pos))
    return false;

  states.push_back(vector<bool>(std::istream_iterator<bool>(iss),
                                std::istream_iterator<bool>()));

  return true;
}

static void
read_meth_table(const string &table_file,
                vector<MSite> &sites,
                vector<vector<bool> > &states) {

  std::ifstream table_in(table_file.c_str());
  if (!table_in)
    throw std::runtime_error("bad table file: " + table_file);

  string line;
  while (getline(table_in, line))
    if (line[0] != '#')
      if (!parse_line(line, sites, states))
        throw std::runtime_error("bad table file line: " + line);
}

static size_t
distance(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom ? max(a.pos, b.pos) - min(a.pos, b.pos) :
    std::numeric_limits<size_t>::max();
}


static void
separate_regions(const size_t desert_size,
                 vector<MSite> &sites, vector<size_t> &reset_points) {
  reset_points.push_back(0);
  for (size_t i = 1; i < sites.size(); ++i)
    if (distance(sites[i - 1], sites[i]) > desert_size)
      reset_points.push_back(i);
  reset_points.push_back(sites.size());
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
count_triads(const vector<size_t> &subtree_sizes,
             const vector<size_t> &parent_ids,
             const vector<vector<bool> > &tree_state_table,
             const vector<size_t> &reset_points,
             vector<triple_state> &triad_counts,
             vector<pair_state> &start_counts,
             pair_state &root_counts,
             pair<double, double> &root_start_counts) {

  triad_counts = vector<triple_state>(subtree_sizes.size());
  start_counts = vector<pair_state>(subtree_sizes.size());
  root_counts = pair_state();

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const size_t start = reset_points[i];
    const size_t end = reset_points[i + 1];

    root_start_counts.first += !tree_state_table[start][0];
    root_start_counts.second += tree_state_table[start][0];

    for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
      const size_t parent = tree_state_table[start][parent_ids[node_id]];
      const size_t curr = tree_state_table[start][node_id];
      start_counts[node_id](parent, curr) += 1.0;
    }

    for (size_t pos = start + 1; pos < end; ++pos) {

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


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile("/dev/stdout");
    size_t desert_size = 1000;
    bool independent_sites = false;

    bool test_pi0 = false;
    bool test_G = false;
    bool test_lambda = false;
    size_t test_branch = 0;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "test harness for "
                           "some functions used in epiphy",
                           "<meth-table-file> <params-file>");
    opt_parse.add_opt("indep", 'I', "sites are independent",
                      false, independent_sites);
    opt_parse.add_opt("desert", 'd',
                      "desert size (default " + std::to_string(1000) + ")",
                      false, desert_size);
    opt_parse.add_opt("pi0", 'p', "optimize initial meth distribution pi0",
                      false, test_pi0);
    opt_parse.add_opt("horiz", 'G', "optimize horizontal transitions G",
                      false, test_G);
    opt_parse.add_opt("lambda", 'l', "optimize vertical rate lambda",
                      false, test_lambda);
    opt_parse.add_opt("branch", 'b', "optimize this branch",
                      false, test_branch);
    opt_parse.add_opt("output", 'o', "name of output file "
                      "(default: stdout)", true, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
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
    const size_t flag_sum = test_pi0 + test_G + test_lambda;
    if (flag_sum > 1 && test_branch > 0) {
      cerr << "pi0, G, lambda and branch are mutually exclusive" << endl;
      return EXIT_SUCCESS;
    }
    if (flag_sum + test_branch == 0) {
      cerr << "select one parameter to test" << endl;
      return EXIT_SUCCESS;
    }
    const string meth_table_file(leftover_args.front());
    const string param_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/


    if (VERBOSE)
      cerr << "[reading parameters]" << endl;
    PhyloTreePreorder t;
    param_set ps;
    ps.read(param_file, t);
    if (VERBOSE)
      cerr << "==> parameters={" << ps << "}" << endl;

    vector<size_t> subtree_sizes;
    t.get_subtree_sizes(subtree_sizes);

    vector<string> leaf_names;
    t.get_leaf_names(leaf_names);

    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);

    const size_t n_nodes = t.get_size();

    if (VERBOSE)
      cerr << "[reading states]" << endl;

    vector<MSite> sites;
    vector<vector<bool> > states;
    read_meth_table(meth_table_file, sites, states);

    if (VERBOSE)
      cerr << "==> n_sites=" << states.size() << endl;

    if (VERBOSE)
      cerr << "[separating by deserts]" << endl;
    vector<size_t> reset_points;
    separate_regions(desert_size, sites, reset_points);
    if (VERBOSE)
      cerr << "==> n_resets=" << reset_points.size() << endl;

    vector<triple_state> triad_counts;
    vector<pair_state> start_counts;
    pair_state root_counts;
    pair<double, double> root_start_counts;

    count_triads(subtree_sizes, parent_ids,
                 states, reset_points,
                 triad_counts,
                 start_counts,
                 root_counts,
                 root_start_counts);

    if (VERBOSE) {
      cerr << "root_counts:\n" << root_counts << endl
           << endl << "start_counts:" << endl;
      for (size_t i = 0; i < start_counts.size(); ++i)
        cerr << start_counts[i] << endl;
      cerr << endl << "triad_counts:" << endl;
      for (size_t i = 0; i < triad_counts.size(); ++i)
        cerr << triad_counts[i] << endl;
    }

    param_set optimized_ps(ps);

    if (test_pi0) {
      if (VERBOSE)
        cerr << "optimizing pi0" << endl;
      max_likelihood_pi0(VERBOSE, root_start_counts, optimized_ps);
      if (VERBOSE)
        cerr << optimized_ps << endl;
    }

    if (test_G) {
      if (VERBOSE)
        cerr << "optimizing G" << endl;
      max_likelihood_horiz(VERBOSE, subtree_sizes,
                           root_counts, triad_counts, optimized_ps);
      if (VERBOSE)
        cerr << optimized_ps << endl;
    }

    if (test_lambda) {
      if (VERBOSE)
        cerr << "optimizing lambda" << endl;
      max_likelihood_rate(VERBOSE, subtree_sizes,
                          start_counts, triad_counts, optimized_ps);
      if (VERBOSE)
        cerr << optimized_ps << endl;
    }

    if (test_branch > 0) {
      if (VERBOSE)
        cerr << "optimizing branch: " << test_branch << endl;
      max_likelihood_branch(VERBOSE, subtree_sizes, test_branch,
                            start_counts, triad_counts, optimized_ps);
    }

    optimized_ps.write(t, outfile);
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
