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
using std::to_string;
using std::vector;
using std::endl;
using std::cerr;
using std::pair;

using std::max;
using std::min;

#include <random>
using std::uniform_real_distribution;

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

template <class T> void
read_meth_table(const string &table_file,
                vector<MSite> &sites,
                vector<string> &species_names,
                vector<vector<T> > &states) {

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
separate_regions(const size_t desert_size, vector<MSite> &sites,
                 vector<pair<size_t, size_t> > &blocks) {
  for (size_t i = 0; i < sites.size(); ++i)
    if (i == 0 || distance(sites[i - 1], sites[i]) > desert_size)
      blocks.push_back(std::make_pair(i, i));
    else blocks.back().second = i;
}


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile("/dev/stdout");
    size_t desert_size = 1000;
    bool independent_sites = false;

    bool test_pi0 = false;
    bool test_G = false;
    bool test_root_G = false;
    bool test_lambda = false;
    size_t test_branch = 0;

    /* model of horizontal rates: 
       1 same across all; 
       2 root vs non-root
       3 node-specific */
    size_t horiz_mode = 2;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "test harness for "
                           "some functions used in epiphy",
                           "<meth-table-file> <params-file>");
    opt_parse.add_opt("indep", 'I', "sites are independent",
                      false, independent_sites);
     opt_parse.add_opt("mode", 'm', "mode of horizontal transition rates\n"
                      "1 -- all nodes share the same parameters\n"
                      "2 -- all non-root nodes share the same parameters\n"
                      "3 -- parameters are node-specific\n"
                      "(default: " + to_string(horiz_mode) + ")",
                      false, horiz_mode);
    opt_parse.add_opt("desert", 'd',
                      "desert size (default " + std::to_string(1000) + ")",
                      false, desert_size);
    opt_parse.add_opt("pi0", 'p', "optimize initial meth distribution pi0",
                      false, test_pi0);
    opt_parse.add_opt("horiz", 'G', "optimize horizontal transitions G",
                      false, test_G);
    opt_parse.add_opt("f0f1", 'f', "optimize root horizontal transitions f0 f1",
                      false, test_root_G);
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
    if (!(horiz_mode == 1 || horiz_mode ==2 || horiz_mode == 3)) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string meth_table_file(leftover_args.front());
    const string param_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (independent_sites)
      desert_size = 0;

    // initialize the random number generator
    std::random_device rd;
    std::mt19937_64 gen(rd());

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

    if (VERBOSE)
      cerr << "[reading states]" << endl;

    vector<MSite> sites;
    vector<vector<bool> > states;
    vector<string> meth_table_species;
    read_meth_table(meth_table_file, sites, meth_table_species, states);

    if (VERBOSE)
      cerr << "==> n_sites=" << states.size() << endl;

    if (VERBOSE)
      cerr << "[separating by deserts]" << endl;
    vector<pair<size_t, size_t> > blocks;
    separate_regions(desert_size, sites, blocks);
    if (VERBOSE)
      cerr << "==> n_resets=" << blocks.size() << endl;

    vector<triple_state> triad_counts;
    vector<pair_state> start_counts;
    pair_state root_counts;
    pair<double, double> root_start_counts;

    count_triads(subtree_sizes, parent_ids,
                 states, blocks,
                 root_start_counts,
                 root_counts,
                 start_counts,
                 triad_counts);

    param_set optimized_ps(ps);

    if (test_pi0) {
      // sample arbitrary value
      optimized_ps.pi0 = std::uniform_real_distribution<>(0.0, 1.0)(gen);
      if (VERBOSE)
        cerr << "optimizing pi0" << endl
             << "initial params: " << optimized_ps << endl;
      max_likelihood_pi0(VERBOSE, root_start_counts, optimized_ps);
      if (VERBOSE)
        cerr << optimized_ps << endl;
    }

    // if (test_root_G) {
    //   // sample arbitrary value (but making it > 0.5)
    //   optimized_ps.g0[0] = std::uniform_real_distribution<>(0.5, 1.0)(gen);
    //   optimized_ps.g1[0] = std::uniform_real_distribution<>(0.5, 1.0)(gen);
    //   if (VERBOSE)
    //     cerr << "optimizing root_G" << endl
    //          << "initial params: " << optimized_ps << endl;
    //   max_likelihood_f0_f1(VERBOSE, root_counts, optimized_ps);
    //   if (VERBOSE)
    //     cerr << optimized_ps << endl;
    // }

    if (test_G) {
      // sample arbitrary value (but making it > 0.5)
      for (size_t i = 0; i < ps.T.size(); ++i) {
        optimized_ps.g0[i] = std::uniform_real_distribution<>(0.5, 1.0)(gen);
        optimized_ps.g1[i] = std::uniform_real_distribution<>(0.5, 1.0)(gen);
      }
      if (VERBOSE)
        cerr << "optimizing G" << endl
             << "initial params: " << optimized_ps << endl;
      max_likelihood_horiz(VERBOSE, horiz_mode, subtree_sizes, root_counts,
                           triad_counts, optimized_ps);
      if (VERBOSE)
        cerr << optimized_ps << endl;
    }

    if (test_lambda) {
      // sample arbitrary value
      // ADS: is this the right way to get a simulated value???
      optimized_ps.rate0 = std::uniform_real_distribution<>(0.0, 1.0)(gen);
      if (VERBOSE)
        cerr << "optimizing lambda" << endl
             << "initial params: " << optimized_ps << endl;
      max_likelihood_rate(VERBOSE, subtree_sizes,
                          start_counts, triad_counts, optimized_ps);
      if (VERBOSE)
        cerr << optimized_ps << endl;
    }

    if (test_branch > 0) {
      // ADS: we need a way to sample a random branch
      if (VERBOSE)
        cerr << "optimizing branch: " << test_branch << endl;
      max_likelihood_branch(VERBOSE, subtree_sizes, test_branch,
                            start_counts, triad_counts, optimized_ps);
    }

    optimized_ps.write(t, outfile);


    if (VERBOSE) {
      cerr << "[sufficient statistics]" << endl;
      cerr << "root_start_counts:\n"
           << root_start_counts.first << '\t'
           << root_start_counts.second << endl
           << "root_counts:\n" << root_counts << endl
           << "start_counts:\n";
      copy(start_counts.begin(), start_counts.end(),
           std::ostream_iterator<pair_state>(cerr, "\n"));
      if (!independent_sites) {
        cerr << "triad_counts:\n";
        copy(triad_counts.begin(), triad_counts.end(),
             std::ostream_iterator<triple_state>(cerr, "\n"));
      }
      const double llk =
        log_likelihood(subtree_sizes, optimized_ps, root_start_counts,
                       root_counts, start_counts, triad_counts);
      cerr << "log_likelihood=" << llk << endl;
    }

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
