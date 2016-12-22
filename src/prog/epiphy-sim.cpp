/*****************************************************************************
 *  epiphy-sim: a program to simulate methylation states for species in
 *  a phylogenetic tree according to inheritance and auto-correlation
 *  of methylation states, and WGBS data sets for extant species.
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

using std::string;
using std::vector;
using std::endl;
using std::cerr;

using std::max;
using std::min;


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


class single_edge_sampler {
public:
  single_edge_sampler(const pair_state &P) { // need to allow seed
    dis = std::uniform_real_distribution<>(0, 1);
    std::random_device rd;
    gen = std::mt19937_64(rd());
    u_gvn_parent_u = P(0, 0)/max(1.0, P(0, 0) + P(0, 1));
    u_gvn_parent_m = P(1, 0)/max(1.0, P(1, 0) + P(1, 1));
  }
  bool operator()(const bool parent_is_u) const {
    return dis(gen) < (parent_is_u ? u_gvn_parent_u : u_gvn_parent_m);
  }
  string tostring() const {
    std::ostringstream oss;
    oss << std::setprecision(3);
    oss << '['
        << u_gvn_parent_u << ',' << 1.0 - u_gvn_parent_u << "]["
        << u_gvn_parent_m << ',' << 1.0 - u_gvn_parent_m << ']';
    return oss.str();
  }
private:
  mutable std::uniform_real_distribution<> dis;
  mutable std::mt19937_64 gen;
  double u_gvn_parent_u;
  double u_gvn_parent_m;
};

std::ostream&
operator<<(std::ostream &out, const single_edge_sampler &ses) {
  return out << ses.tostring();
}


class two_edge_sampler {
public:
  two_edge_sampler(const triple_state &GP) { // need to allow seed
    dis = std::uniform_real_distribution<>(0, 1);
    std::random_device rd;
    gen = std::mt19937_64(rd());

    u_gvn_left_u_par_u = GP(0, 0, 0)/max(1.0, GP(0, 0, 0) + GP(0, 0, 1));
    u_gvn_left_u_par_m = GP(0, 1, 0)/max(1.0, GP(0, 1, 0) + GP(0, 1, 1));
    u_gvn_left_m_par_u = GP(1, 0, 0)/max(1.0, GP(1, 0, 0) + GP(1, 0, 1));
    u_gvn_left_m_par_m = GP(1, 1, 0)/max(1.0, GP(1, 1, 0) + GP(1, 1, 1));
  }
  bool operator()(const bool left_is_u, const bool par_is_u) const {
    return dis(gen) < (left_is_u ?
                       (par_is_u ? u_gvn_left_u_par_u : u_gvn_left_u_par_m) :
                       (par_is_u ? u_gvn_left_m_par_u : u_gvn_left_m_par_m));
  }
  string tostring() const {
    std::ostringstream oss;
    oss << std::setprecision(3);
    oss << '['
        << u_gvn_left_u_par_u << ',' << 1.0 - u_gvn_left_u_par_u << "]["
        << u_gvn_left_u_par_m << ',' << 1.0 - u_gvn_left_u_par_m << "]\n["
        << u_gvn_left_m_par_u << ',' << 1.0 - u_gvn_left_m_par_u << "]["
        << u_gvn_left_m_par_m << ',' << 1.0 - u_gvn_left_m_par_m << ']';
    return oss.str();
  }
private:
  mutable std::uniform_real_distribution<> dis;
  mutable std::mt19937_64 gen;

  double u_gvn_left_u_par_u;
  double u_gvn_left_u_par_m;
  double u_gvn_left_m_par_u;
  double u_gvn_left_m_par_m;
};

std::ostream&
operator<<(std::ostream &out, const two_edge_sampler &tes) {
  return out << tes.tostring();
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
simulate_site_start(const vector<size_t> &subtree_sizes,
                    const vector<single_edge_sampler> &Psamp,
                    const size_t node_id, vector<bool> &states) {
  if (!is_leaf(subtree_sizes[node_id])) {
    const bool current_state = states[node_id];
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      states[child_id] = Psamp[child_id](current_state);
      simulate_site_start(subtree_sizes, Psamp, child_id, states);
      count += subtree_sizes[child_id];
    }
  }
}

static void
simulate_site_start(const vector<size_t> &subtree_sizes, const param_set &ps,
                    const vector<single_edge_sampler> &Psamp, vector<bool> &states) {

  states.resize(subtree_sizes.size());

  // sample root at start position
  std::random_device rd;
  std::mt19937_64 gen(rd());
  states[0] = std::uniform_real_distribution<>(0, 1)(gen) < ps.pi0;

  simulate_site_start(subtree_sizes, Psamp, 0, states);
}


static void
simulate_site(const vector<size_t> &subtree_sizes,
              const vector<two_edge_sampler> &GPsamp,
              const vector<bool> &prev_states, const size_t node_id,
              vector<bool> &states) {

  if (!is_leaf(subtree_sizes[node_id])) {
    const bool current_state = states[node_id];
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      states[child_id] = GPsamp[child_id](prev_states[child_id], current_state);
      simulate_site(subtree_sizes, GPsamp, prev_states, child_id, states);
      count += subtree_sizes[child_id];
    }
  }
}

static void
simulate_site(const vector<size_t> &subtree_sizes,
              const single_edge_sampler &Gsamp,
              const vector<two_edge_sampler> &GPsamp,
              const vector<bool> &prev_states, vector<bool> &states) {

  states.resize(subtree_sizes.size());
  states[0] = Gsamp(prev_states[0]); // sample root using state to left

  simulate_site(subtree_sizes, GPsamp, prev_states, 0, states);
}


int main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;
    size_t desert_size = 1000;
    bool independent_sites = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate methylomes "
                           "related by a phylogenetic tree",
                           "<sites-file> <params-file>");
    opt_parse.add_opt("indep", 'I', "sites are independent",
                      false, independent_sites);
    opt_parse.add_opt("desert", 'd',
                      "desert size (default " + std::to_string(1000) + ")",
                      false, desert_size);
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
    const string cpgs_file(leftover_args.front());
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

    const size_t n_nodes = t.get_size();

    if (VERBOSE)
      cerr << "[computing transition matrices]" << endl;

    pair_state G(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);
    vector<pair_state> P;
    vector<triple_state> GP;
    get_transition_matrices(ps, P, GP);

    single_edge_sampler Gsamp(G);
    vector<single_edge_sampler> Psamp;
    vector<two_edge_sampler> GPsamp;
    for (size_t i = 0; i < n_nodes; ++i) {
      Psamp.push_back(single_edge_sampler(P[i]));
      GPsamp.push_back(two_edge_sampler(GP[i]));
    }

    if (VERBOSE) {
      cerr << "==> horizontal transitions:" << endl;
      cerr << Gsamp << endl;
      cerr << "==> vertical transitions:" << endl;
      for (size_t i = 0; i < n_nodes; ++i)
        cerr << Psamp[i] << endl;
      cerr << "==> combined transitions:" << endl;
      for (size_t i = 0; i < n_nodes; ++i)
        cerr << GPsamp[i] << endl;
    }


    if (VERBOSE)
      cerr << "[loading sites]" << endl;
    std::ifstream sites_in(cpgs_file.c_str());
    if (!sites_in)
      throw std::runtime_error("bad input file: " + cpgs_file);
    vector<MSite> sites;
    MSite site;
    while (sites_in >> site) sites.push_back(site);
    const size_t n_sites = sites.size();
    if (VERBOSE)
      cerr << "==> n_sites=" << n_sites << endl;

    if (VERBOSE)
      cerr << "[separating by deserts]" << endl;
    vector<size_t> reset_points;
    separate_regions(desert_size, sites, reset_points);
    if (VERBOSE)
      cerr << "==> n_resets=" << reset_points.size() << endl;

    if (VERBOSE)
      cerr << "[simulating]" << endl;

    vector<vector<bool> > states(n_sites);
    for (size_t i = 0; i < reset_points.size() - 1; ++i) {
      const size_t start = reset_points[i];
      const size_t end = reset_points[i + 1]; // !!!!
      assert(start < end);
      simulate_site_start(subtree_sizes, ps, Psamp, states[start]);
      for (size_t pos = start + 1; pos < end; ++pos)
        simulate_site(subtree_sizes, Gsamp, GPsamp, states[pos - 1], states[pos]);
    }

    vector<size_t> leaves_preorder;
    subtree_sizes_to_leaves_preorder(subtree_sizes, leaves_preorder);

    std::ofstream out(outfile.c_str());
    for (size_t i = 0; i < states.size(); ++i) {
      out << sites[i].chrom << ':' << sites[i].pos;
      for (size_t j = 0; j < states[i].size(); ++j)
        out << '\t' << states[i][j];
      out << endl;
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
