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
#include <random>

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"


/* from methpipe */
#include "MethpipeFiles.hpp"

/* headers for epigenomic evolution */
#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"
#include "param_set.hpp"
#include "epiphy_utils.hpp"
#include "sufficient_statistics_helpers.hpp"

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

std::random_device rd; //seed generator
//std::mt19937_64 gen(rd()); //generator initialized with seed from rd
std::mt19937_64 gen(0); //generator initialized with seed from rd

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////////////    FOR LOADING DATA     //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// methtab format: COL1=chr:start, COLn=HypoProb (-1 if no coverage)
struct Site {
  string chrom;
  size_t pos;

  Site() {}
  Site(const string &chr, size_t position) :
    chrom(chr), pos(position) {}
  Site(const string &line) {
    const size_t chrom_end = line.find_first_of(':');
    chrom = line.substr(0, chrom_end);
    const size_t pos_end = line.find_first_of(':', chrom_end + 1);
    pos = stoi(line.substr(chrom_end + 1, pos_end - chrom_end));
  }

  void
  assign_from_methtab_first_column(const string &line) {
    const size_t chrom_end = line.find_first_of(':');
    chrom = line.substr(0, chrom_end);
    const size_t pos_end = line.find_first_of(':', chrom_end + 1);
    pos = stoi(line.substr(chrom_end + 1, pos_end - chrom_end));
  }
  static size_t
  distance(const Site &s, const Site &t) {
    size_t d = std::numeric_limits<size_t>::max();
    if (s.chrom == t.chrom)
      d = (t.pos > s.pos) ? t.pos - s.pos : s.pos - t.pos;
    return d;
  }
};


static bool
parse_table_line(istringstream &iss, vector<vector<double> > &meth) {
  double val;
  vector<double> meth_fields;
  size_t observed = 0;

  while (iss >> val) {
    // hypomethylation probability should be -1.0 or [0,1] -1.0
    // indicates missing data
    if (!valid_meth_value(val))
      throw SMITHLABException("bad probability value: [" + iss.str() + "]");
    if (!missing_meth_value(val)) {
      ++observed;
      val = away_from_extremes(val, PROBABILITY_GUARD);
    }
    meth_fields.push_back(val);
  }

  if (meth_fields.empty())
    throw SMITHLABException("bad line format: " + iss.str());

  if (observed > 0)
    meth.push_back(meth_fields);

  return (observed > 0);
}


template <class T> void
parse_meth_table_line(const string &line, vector<Site> &sites,
                      vector<vector<T> > &meth) {

  istringstream iss(line);

  // take care of the site location (chrom and position)
  string first_column;
  iss >> first_column;
  Site site(first_column);

  // now get the methylation information, either as pairs of read
  // counts or posteriors, depending on the file type (run mode)
  if (parse_table_line(iss, meth))
    sites.push_back(site);
}


template <class T> void
load_meth_table(const string &filename, vector<Site> &sites,
                vector<vector<T> > &meth,
                vector<string> &species_names) {

  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot read: " + filename);

  // get the species from the header
  string line;
  getline(in, line);
  istringstream iss(line);
  string token;
  while (iss >> token)
    species_names.push_back(token);

  const size_t n_species = species_names.size();
  while (getline(in, line)) {
    parse_meth_table_line(line, sites, meth);
    if (meth.back().size() != n_species)
      throw SMITHLABException("inconsistent line length: " + line);
  }
}


static void
load_full_states_as_prob(const string filename, const size_t n_nodes,
                         vector<Site> &sites,
                         vector<vector<double> > &tree_prob_table) {
  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot read:" + filename);

  string line;
  while (getline(in, line)) {
    if (line[0] != '#') { //skip header lines
      istringstream iss(line);
      string first_column;
      iss >> first_column;
      sites.push_back(Site(first_column));

      string state;
      iss >> state;
      vector<double> probs;
      if (state.length() != n_nodes)
        throw SMITHLABException("inconsistent state length: " + line);

      for (size_t i = 0; i < n_nodes; ++i)
        probs.push_back(away_from_extremes(state[i] == '0' ? 1.0 : 0.0,
                                           PROBABILITY_GUARD));

      tree_prob_table.push_back(probs);
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////// ABOVE CODE IS FOR LOADING DATA  //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////// BELOW CODE IS FOR PROCESSING TREE ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
get_parent_id(const vector<size_t> &subtree_sizes, vector<size_t> &parent_id) {
  parent_id = vector<size_t>(subtree_sizes.size(), 0);
  for (size_t i = 0; i < subtree_sizes.size(); ++i)
    for (size_t count = 1; count < subtree_sizes[i];) {
      const size_t child_id = i + count;
      parent_id[child_id] = i;
      count += subtree_sizes[child_id];
    }
}

static bool
has_same_species_order(const PhyloTreePreorder &the_tree,
                       const vector<string> &meth_table_species) {
  vector<string> leaf_names;
  the_tree.get_leaf_names(leaf_names);
  return leaf_names == meth_table_species;
}

static void
fill_leaf_prob(const bool VERBOSE,
               const vector<Site> &sites, const size_t desert_size,
               vector<vector<double> > &hypo_prob_table, double &pi0_est) {

  const size_t n_sites = hypo_prob_table.size();
  const size_t n_leaf = hypo_prob_table[0].size();

  pi0_est = 0.0;
  for (size_t i = 0; i < n_leaf; ++i) {
    double mean_hypo_prob = 0.0;
    size_t n_obs = 0;
    for (size_t j = 0; j < n_sites; ++j)
      if (hypo_prob_table[j][i] >= 0.0) {
        mean_hypo_prob += hypo_prob_table[j][i];
        ++n_obs;
      }
    mean_hypo_prob = mean_hypo_prob/n_obs;
    pi0_est += mean_hypo_prob;

    if (VERBOSE)
      cerr << "[mean_hypo_prob: " << mean_hypo_prob << "]" << endl;

    // ADS: not clear on what is the purpose of the loop below
    size_t prev = 0;
    size_t next = 0;
    size_t count = 0;
    for (size_t j = 0; j < n_sites; ++j) {
      if (hypo_prob_table[j][i] >= 0) {
        prev = j;
        next = j;
      } else {
        if (j >= next) {
          next = j + 1;
          while (next < n_sites && hypo_prob_table[next][i] < 0.0)
            ++next;
        }
        const size_t d1 = Site::distance(sites[j], sites[prev]);
        const size_t d2 = (next < n_sites) ?
          Site::distance(sites[j], sites[next]) : desert_size;

        if (prev < j && j < next && d1 < desert_size && d2 < desert_size) {
          const double w1 = static_cast<double>(d2)/(d1 + d2);
          hypo_prob_table[j][i] = (hypo_prob_table[prev][i]*w1 +
                                   hypo_prob_table[next][i]*(1.0 - w1));
        } else if (prev < j && d1 < desert_size) {
          const double w1 = d1/desert_size;
          hypo_prob_table[j][i] = (mean_hypo_prob*w1 +
                                   hypo_prob_table[prev][i]*(1.0 - w1));
        } else if (d2 < desert_size) {
          const double w1 = d2/desert_size;
          hypo_prob_table[j][i] = (mean_hypo_prob*w1 +
                                   hypo_prob_table[next][i]*(1.0 - w1));
        } else {
          hypo_prob_table[j][i] = mean_hypo_prob;
        }
        ++count;
      }
    }
    if (VERBOSE)
      cerr << "[filled " << count << " "
           << "missing sites in leaf: " << i << "]" << endl;
  }
  pi0_est = pi0_est/n_leaf;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////// PHYLOTREEPREORDER HELPERS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
subtree_sizes_to_leaves_preorder(const vector<size_t> &subtree_sizes,
                                 vector<size_t> &leaves_preorder) {
  for (size_t i = 0; i < subtree_sizes.size(); ++i)
    if (subtree_sizes[i] == 1)
      leaves_preorder.push_back(i);
}


bool
is_binary(const vector<size_t> &subtree_sizes) {
  // ADS: this function seems not to be used
  return (subtree_sizes[0] == 1 + subtree_sizes[1] +
          subtree_sizes[subtree_sizes[1] + 1]);
}

static void
get_node_degrees(const vector<size_t> &subtree_sizes,
                 const size_t tree_start,
                 vector<size_t> &degrees) {

  assert(degrees.size() == subtree_sizes.size());

  if (subtree_sizes[tree_start] == 1) {
    degrees[tree_start] = 1;
  } else {
    if (tree_start > 0)
      degrees[tree_start] = 1;

    for (size_t count = 1; count < subtree_sizes[tree_start];) {
      const size_t next_child = tree_start + count;
      get_node_degrees(subtree_sizes, next_child, degrees);
      degrees[tree_start] += 1;
      count += subtree_sizes[next_child];
    }
  }
}


static void
get_degrees(const vector<size_t> &subtree_sizes, vector<size_t> &degrees) {
  degrees = vector<size_t>(subtree_sizes.size(), 0);
  get_node_degrees(subtree_sizes, 0, degrees);
}


bool
is_semi_binary(const vector<size_t> &degrees) {

  if (degrees[0] < 2 || degrees[0] > 3)
    return false;

  for (size_t i = 1; i < degrees.size(); ++i )
    if (degrees[i] != 1 && degrees[i] != 3)
      return false;

  return true;
}

static bool
is_root(const size_t node_id) {return node_id == 0;}

static bool
is_leaf(const size_t subtree_size) {return subtree_size == 1;}

static size_t
count_leaves(const vector<size_t> &subtree_sizes) {
  size_t n_leaf = 0;
  for (size_t i = 0; i < subtree_sizes.size(); ++i)
    n_leaf += is_leaf(subtree_sizes[i]);
  return n_leaf;
}

static void
get_children(const size_t node_id, const vector<size_t> &subtree_sizes,
             vector<size_t> &children) {
  for (size_t c = 1; c < subtree_sizes[node_id]; c += subtree_sizes[node_id + c])
    children.push_back(node_id + c);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/// CODE BELOW IS FOR SEPARATING SITES THAT ARE FAR FROM EACH OTHER ////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

static void
separate_regions(const bool VERBOSE, const size_t desert_size,
                 const size_t min_sites_per_block,
                 vector<vector<double> > &meth,
                 vector<vector<double> > &meth_full_leaf,
                 vector<Site> &sites,
                 vector<size_t> &reset_points,
                 double &g0_est, double &g1_est, double &pi0_est) {

  const size_t n_leaves = meth[0].size();
  const size_t n_sites = sites.size();

  vector<vector<double> > meth_aug = meth;
  fill_leaf_prob(VERBOSE, sites, desert_size, meth_aug, pi0_est);

  //scan desert
  //uncovered site has prob value set to -1.0
  vector<bool> is_desert(n_sites, false);
  for (size_t i = 0; i < n_leaves; ++i) {
    if (VERBOSE)
      cerr << "[processing sample: " << i;
    size_t j = 0;
    while (j < n_sites && meth[j][i] == -1.0) ++j;

    if (j == n_sites)
      throw SMITHLABException("No valid observation.");

    if (Site::distance(sites[0], sites[j]) > desert_size )
      for (size_t w = 0; w < j; ++w)
        is_desert[w] = true;

    Site prev_obs = sites[j];
    for (size_t k = j+1; k < n_sites; ++k) {

      while (k < n_sites && meth[k][i] == -1.0) ++k;

      if (k < n_sites) {
        if (Site::distance(prev_obs, sites[k]) > desert_size)
          for (size_t w = j + 1; w < k; ++w)
            is_desert[w] = true;
        j = k;
        prev_obs = sites[j];
      } else {
        for (size_t w = j + 1; w < n_sites; ++w)
          is_desert[w] = true;
      }
    }
    if (VERBOSE)
      cerr << " done]" << endl;
  }

  size_t good = 0;
  for (size_t i = 0; i < n_sites; ++i)
    if (!is_desert[i]) ++good;

  if (VERBOSE)
    cerr << "[n_good_sites: " << good << "]" << endl;

  vector<Site> sites_copy;
  Site last_site;
  vector<vector<double> > meth_copy;
  meth_full_leaf.clear();

  size_t last_block_size = 0;
  for (size_t i = 0; i < n_sites; ++i) {
    if (!is_desert[i]) {
      // Add new site
      sites_copy.push_back(sites[i]);
      meth_full_leaf.push_back(meth_aug[i]);
      meth_copy.push_back(meth[i]);

      ++last_block_size;
      // Maintain blocks
      if (sites_copy.size() == 1) {
        reset_points.push_back(0);
      } else if (Site::distance(last_site, sites[i]) > desert_size) {
        if (last_block_size < min_sites_per_block) { //remove small blocks
          const size_t start = reset_points.back();
          sites_copy.erase(sites_copy.begin() + start, sites_copy.end());
          meth_copy.erase(meth_copy.begin() + start, meth_copy.end());
          meth_full_leaf.erase(meth_full_leaf.begin() + start, meth_full_leaf.end());
          reset_points.pop_back();
          last_block_size =
            reset_points.empty() ? 0 : meth_copy.size() - reset_points.back();
        } else { //end block, and start new block
          reset_points.push_back(sites_copy.size() - 1);
          last_block_size = 1;
        }
      }
      if (sites_copy.size() > 0)
        last_site = sites_copy.back();
    }
  }

  if (last_block_size < min_sites_per_block) { //remove small blocks
    const size_t start = reset_points.back();
    sites_copy.erase(sites_copy.begin() + start, sites_copy.end() );
    meth_copy.erase(meth_copy.begin() + start, meth_copy.end());
    meth_full_leaf.erase(meth_full_leaf.begin() + start, meth_full_leaf.end());
    reset_points.pop_back();
    last_block_size =
      reset_points.empty() ? 0 : meth_copy.size() - reset_points.back();
  } else {
    reset_points.push_back(sites_copy.size());
  }

  meth.swap(meth_copy);
  sites.swap(sites_copy);

  // ADS: not sure why this code is in this function
  double g00 = 0.0, g01 = 0.0, g10 = 0.0, g11 = 0.0;
  for (size_t i = 0; i < n_leaves; ++i)
    for (size_t j = 0; j < meth.size() - 1; ++j)
      if (meth[j][i] >= 0 && meth[j + 1][i] >= 0) {
        g00 += meth[j][i]*meth[j + 1][i];
        g01 += meth[j][i]*(1.0 - meth[j + 1][i]);
        g10 += (1.0 - meth[j][i])*meth[j + 1][i];
        g11 += (1.0 - meth[j][i])*(1.0 - meth[j + 1][i]);
      }

  g0_est = (g00 + 1.0)/(g01 + g00 + 1.0);
  g1_est = (g11 + 1.0)/(g11 + g10 + 1.0);
}


//for complete data, assume each block have enough sites, only set reset_points
static void
separate_regions(const size_t desert_size,
                 vector<Site> &sites, vector<size_t> &reset_points) {
  reset_points.push_back(0);
  for (size_t i = 1; i < sites.size(); ++i)
    if (Site::distance(sites[i-1], sites[i]) > desert_size)
      reset_points.push_back(i);
  reset_points.push_back(sites.size());
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////       Units grouped       ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


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


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////      DATA AUGMENTATION        /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// uninitialized nodes have value -1.0
// assuming all leaves are initialized already
static void
init_internal_prob(const vector<size_t> &subtree_sizes, const size_t node_id,
                   vector<double> &node_probs) {
  // Recursion
  size_t n_desc_leaf = 0;
  double sum_desc_leaf_prob = 0;
  for (size_t count = 1; count < subtree_sizes[node_id]; ) {
    const size_t child_id = node_id + count;
    if (node_probs[child_id] < 0) // ADS: is this for missing value? [missing_meth_value(const double x)]
      init_internal_prob(subtree_sizes, child_id, node_probs);
    count += subtree_sizes[child_id];
  }

  double direction = 0.0;
  for (size_t i = 1; i < subtree_sizes[node_id]; ++i) {
    if (subtree_sizes[i + node_id]==1) {
      ++n_desc_leaf;
      sum_desc_leaf_prob += node_probs[i + node_id];
      direction += (node_probs[i + node_id] > 0.5) ? 1.0 : -1.0;
    }
  }

  if (node_id != 0 && abs(direction) < n_desc_leaf) {
    double sum_extra_leaf_hypo_prob = 0.0;
    size_t n_extra_leaves = 0;
    // children disagree, find majority state of leaf species outside this subtree
    for (size_t i = 0; i < subtree_sizes.size(); ++i)
      if (is_leaf(subtree_sizes[i]) &&
          (i < node_id || i >= node_id + subtree_sizes[node_id])) {
        sum_extra_leaf_hypo_prob += node_probs[i];
        ++n_extra_leaves;
      }

    if (n_extra_leaves > 0)
      sum_desc_leaf_prob += sum_extra_leaf_hypo_prob/n_extra_leaves;

    ++n_desc_leaf;
  }

  node_probs[node_id] = away_from_extremes(sum_desc_leaf_prob/n_desc_leaf,
                                            PROBABILITY_GUARD);
}

static void
copy_leaf_to_tree_prob(const vector<size_t> &subtree_sizes,
                       const vector<vector<double> > &meth_prob_table,
                       vector<vector<double> > &tree_prob_table) {
  vector<size_t> leaves_preorder;
  subtree_sizes_to_leaves_preorder(subtree_sizes, leaves_preorder);

  const size_t n_leaves = meth_prob_table.front().size();
  for (size_t i = 0; i < meth_prob_table.size(); ++i)
    for (size_t j = 0; j < n_leaves; ++j)
      if (valid_probability(meth_prob_table[i][j]))
        tree_prob_table[i][leaves_preorder[j]] =
          away_from_extremes(meth_prob_table[i][j],
                             PROBABILITY_GUARD);
}

static void
leaf_to_tree_prob(const vector<size_t> &subtree_sizes,
                  const vector<vector<double> > &meth_prob_table,
                  vector<vector<double> > &tree_prob_table) {
  // copy leaves first
  copy_leaf_to_tree_prob(subtree_sizes, meth_prob_table, tree_prob_table);
  // initialize internal probs
  for (size_t i = 0; i < tree_prob_table.size(); ++i)
    init_internal_prob(subtree_sizes, 0, tree_prob_table[i]);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////         APPROX POSTERIOR        /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// The next several functions are for estimating the marginal
// posterior for each site and each node. In each case, we want to
// compute the probability for a given state based on all possible
// combinations of states in the Markov blanket. The functions are
// organized according to whether the current site is the "start", the
// "middle" or the "end". These situations are for the first site, a
// regular site (with neighbors on both sides), or the last
// site. Within each of these functions, there are three major
// conditions corresponding to the current site being at the root, a
// leaf, or an internal node of the tree.

static double
val_or_flip(const double val, int flip) {
  return flip ? 1.0 - val : val;
}

static void
posterior_start(const bool update_leaves,
                const vector<size_t> &subtree_sizes, const pair_state &G,
                const vector<pair_state> &P, const vector<triple_state > &GP,
                const size_t pos, const size_t node_id, const size_t parent_id,
                const vector<double> &obs_meth_curr,
                vector<double> &meth_curr, vector<double> &meth_next,
                double &diff) {
  if (is_leaf(subtree_sizes[node_id])) { //leaf
    if (update_leaves && missing_meth_value(obs_meth_curr[node_id])) {
      const double p_par_next_orig = meth_next[parent_id];
      const double p_par_orig = meth_curr[parent_id];
      const double p_next_orig = meth_next[node_id];

      double p_cur = 0.0;
      for (size_t next = 0; next < 2; ++next) {
        const double p_next = val_or_flip(p_next_orig, next);
        for (size_t par_next = 0; par_next < 2; ++par_next) {
          const double p_par_next = val_or_flip(p_par_next_orig, par_next);
          for (size_t par = 0; par < 2; ++par) {
            const double p_par = val_or_flip(p_par_orig, par);
            const double marginal = p_next*p_par_next*p_par;
            const double p_cur_u =
              P[node_id](par, 0)*GP[node_id](0, par_next, next);
            const double p_cur_m =
              P[node_id](par, 1)*GP[node_id](1, par_next, next);
            const double p_cur_mb = p_cur_u/(p_cur_u + p_cur_m);
            p_cur += p_cur_mb*marginal;
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    }
  } else { // current node is internal
    vector<size_t> children;
    get_children(node_id, subtree_sizes, children);

    for (size_t i = 0; i < children.size(); ++i)
      posterior_start(update_leaves, subtree_sizes, G, P, GP, pos,
                      children[i], node_id, obs_meth_curr,
                      meth_curr, meth_next, diff);

    if (node_id == 0) { //root node  //root has 3 children!
      const double p_next_orig = meth_next[node_id];
      const double p_c0_orig = meth_curr[children[0]];
      const double p_c1_orig = meth_curr[children[1]];
      const double p_c2_orig = meth_curr[children[2]];

      double p_cur = 0.0;
      for (size_t c0 = 0; c0 < 2; ++c0) {
        const double p_c0 = val_or_flip(p_c0_orig, c0);
        for (size_t c1 = 0; c1 < 2; ++c1) {
          const double p_c1 = val_or_flip(p_c1_orig, c1);
          for (size_t c2 = 0; c2 < 2; ++c2) {
            const double p_c2 = val_or_flip(p_c2_orig, c2);
            for (size_t next = 0; next < 2; ++next) {
              const double p_next = val_or_flip(p_next_orig, next);
              const double marginal = p_c0*p_c1*p_c2*p_next;
              const double p_cur_u = (P[children[0]](0, c0)*
                                      P[children[1]](0, c1)*
                                      P[children[2]](0, c2)*G(0, next));
              const double p_cur_m = (P[children[0]](1, c0)*
                                      P[children[1]](1, c1)*
                                      P[children[2]](1, c2)*G(1, next));
              p_cur += p_cur_u/(p_cur_m + p_cur_u)*marginal;
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    } else { // internal node // two children
      const double p_c0_orig = meth_curr[children[0]];
      const double p_c1_orig = meth_curr[children[1]];
      const double p_next_orig = meth_next[node_id];
      const double p_par_orig = meth_curr[parent_id];
      const double p_par_next_orig = meth_next[parent_id];

      double p_cur = 0.0;
      for (size_t c0 = 0; c0 < 2; ++c0) {
        const double p_c0 = val_or_flip(p_c0_orig, c0);
        for (size_t c1 = 0; c1 < 2; ++c1) {
          const double p_c1 = val_or_flip(p_c1_orig, c1);
          for (size_t next = 0; next < 2; ++next) {
            const double p_next = val_or_flip(p_next_orig, next);
            for (size_t par_cur = 0; par_cur < 2; ++par_cur) {
              const double p_par = val_or_flip(p_par_orig, par_cur);
              for (size_t par_next = 0; par_next < 2; ++par_next) {
                const double p_par_next = val_or_flip(p_par_next_orig, par_next);
                const double marginal = p_c0*p_c1*p_next*p_par*p_par_next;
                const double p_cur_u = (P[node_id](par_cur, 0)*
                                        P[children[0]](0, c0)*
                                        P[children[1]](0, c1)*
                                        GP[node_id](0, par_next, next));
                const double p_cur_m = (P[node_id](par_cur, 1)*
                                        P[children[0]](1, c0)*
                                        P[children[1]](1, c1)*
                                        GP[node_id](1, par_next, next));
                p_cur += p_cur_u/(p_cur_u + p_cur_m)*marginal;
              }
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    }
  }
}


static void
posterior_middle(const bool update_leaves,
                 const vector<size_t> &subtree_sizes, const pair_state &G,
                 const vector<pair_state> &P, const vector<triple_state> &GP,
                 const size_t pos, const size_t node_id, const size_t parent_id,
                 const vector<double> &obs_meth_curr,
                 vector<double> &meth_prev, vector<double> &meth_curr,
                 vector<double> &meth_next, double &diff) {

  double p_cur = 0.0;
  if (is_leaf(subtree_sizes[node_id])) { // leaf
    if (update_leaves && missing_meth_value(obs_meth_curr[node_id])) {
      const double p_prev_orig = meth_prev[node_id];
      const double p_next_orig = meth_next[node_id];
      const double p_par_orig = meth_curr[parent_id];
      const double p_par_next_orig = meth_next[parent_id];

      for (size_t prev = 0; prev < 2; ++prev) {
        const double p_prev = val_or_flip(p_prev_orig, prev);
        for (size_t next = 0; next < 2; ++next) {
          const double p_next = val_or_flip(p_next_orig, next);
          for (size_t par_next = 0; par_next <2; ++par_next) {
            const double p_par_next = val_or_flip(p_par_next_orig, par_next);
            for (size_t par = 0; par < 2; ++par) {
              const double p_par = val_or_flip(p_par_orig, par);
              const double marginal = p_prev*p_next*p_par_next*p_par;
              const double p_cur_u = (GP[node_id](prev, par, 0)*
                                      GP[node_id](0, par_next, next));
              const double p_cur_m = (GP[node_id](prev, par, 1)*
                                      GP[node_id](1, par_next, next));
              const double p_cur_mb = p_cur_u/(p_cur_u + p_cur_m);
              p_cur += p_cur_mb*marginal;
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    }
  } else {
    // get children and recurse
    vector<size_t> children;
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      posterior_middle(update_leaves, subtree_sizes, G, P, GP, pos, child_id,
                       node_id, obs_meth_curr,
                       meth_prev, meth_curr, meth_next, diff);
      children.push_back(child_id);
      count += subtree_sizes[child_id];
    }

    if (is_root(node_id)) { //root node  //root has 3 children!
      const double p_prev_orig = meth_prev[node_id];
      const double p_prev_c0_orig = meth_prev[children[0]];
      const double p_prev_c1_orig = meth_prev[children[1]];
      const double p_prev_c2_orig = meth_prev[children[2]];
      const double p_c0_orig = meth_curr[children[0]];
      const double p_c1_orig = meth_curr[children[1]];
      const double p_c2_orig = meth_curr[children[2]];
      const double p_next_orig = meth_next[node_id];

      for (size_t prev = 0; prev < 2; ++prev) {
        const double p_prev = val_or_flip(p_prev_orig, prev);
        for (size_t prev_c0 = 0; prev_c0 < 2; ++prev_c0) {
          const double p_prev_c0 = val_or_flip(p_prev_c0_orig, prev_c0);
          for (size_t prev_c1 = 0; prev_c1 < 2; ++prev_c1) {
            const double p_prev_c1 = val_or_flip(p_prev_c1_orig, prev_c1);
            for (size_t prev_c2 = 0; prev_c2 < 2; ++prev_c2) {
              const double p_prev_c2 = val_or_flip(p_prev_c2_orig, prev_c2);
              for (size_t c0 = 0; c0 < 2; ++c0) {
                const double p_c0 = val_or_flip(p_c0_orig, c0);
                for (size_t c1 = 0; c1 < 2; ++c1) {
                  const double p_c1 = val_or_flip(p_c1_orig, c1);
                  for (size_t c2 = 0; c2 < 2; ++c2) {
                    const double p_c2 = val_or_flip(p_c2_orig, c2);
                    for (size_t next = 0; next < 2; ++next) {
                      const double p_next = val_or_flip(p_next_orig, next);
                      const double marginal = (p_prev*p_prev_c0*p_prev_c1*
                                               p_prev_c2*p_c0*p_c1*p_c2*p_next);
                      const double p_cur_u = (GP[children[0]](prev, 0, c0)*
                                              GP[children[1]](prev, 0, c1)*
                                              GP[children[2]](prev, 0, c2)*
                                              G(0, next)*G(prev, 0));
                      const double p_cur_m = (GP[children[0]](prev, 1, c0)*
                                              GP[children[1]](prev, 1, c1)*
                                              GP[children[2]](prev, 1, c2)*
                                              G(1, next)*G(prev, 1));
                      p_cur += p_cur_u/(p_cur_m + p_cur_u)*marginal;
                    }
                  }
                }
              }
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    } else { // internal node // two children
      const double p_prev_c0_orig = meth_prev[children[0]];
      const double p_prev_c1_orig = meth_prev[children[1]];
      const double p_prev_orig = meth_prev[node_id];
      const double p_c0_orig = meth_curr[children[0]];
      const double p_c1_orig = meth_curr[children[1]];
      const double p_next_orig = meth_next[node_id];
      const double p_par_orig = meth_curr[parent_id];
      const double p_par_next_orig = meth_next[parent_id];

      for (size_t prev_c0 = 0; prev_c0 < 2; ++prev_c0) {
        const double p_prev_c0 = val_or_flip(p_prev_c0_orig, prev_c0);
        for (size_t prev_c1 = 0; prev_c1 < 2; ++prev_c1) {
          const double p_prev_c1 = val_or_flip(p_prev_c1_orig, prev_c1);
          for (size_t prev = 0; prev < 2; ++prev) {
            const double p_prev = val_or_flip(p_prev_orig, prev);
            for (size_t c0 = 0; c0 < 2; ++c0) {
              const double p_c0 = val_or_flip(p_c0_orig, c0);
              for (size_t c1 = 0; c1 < 2; ++c1) {
                const double p_c1 = val_or_flip(p_c1_orig, c1);
                for (size_t next = 0; next < 2; ++next) {
                  const double p_next = val_or_flip(p_next_orig, next);
                  for (size_t par = 0; par < 2; ++par) {
                    const double p_par = val_or_flip(p_par_orig, par);
                    for (size_t par_next = 0; par_next < 2; ++par_next) {
                      const double p_par_next =
                        val_or_flip(p_par_next_orig, par_next);
                      const double marginal = (p_prev_c0*p_prev_c1*p_prev*p_c0*
                                               p_c1*p_next*p_par*p_par_next);
                      const double p_cur_u = (GP[node_id](prev, par, 0)*
                                              GP[children[0]](prev_c0, 0, c0)*
                                              GP[children[1]](prev_c1, 0, c1)*
                                              GP[node_id](0, par_next, next));
                      const double p_cur_m = (GP[node_id](prev, par, 1)*
                                              GP[children[0]](prev_c0, 1, c0)*
                                              GP[children[1]](prev_c1, 1, c1)*
                                              GP[node_id](1, par_next, next));
                      p_cur += p_cur_u/(p_cur_u + p_cur_m)*marginal;
                    }
                  }
                }
              }
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    }
  }
}


static void
posterior_end(const bool update_leaves,
              const vector<size_t> &subtree_sizes, const pair_state &G,
              const vector<pair_state> &P, const vector<triple_state> &GP,
              const size_t pos, const size_t node_id, const size_t parent_id,
              const vector<double> &obs_meth_curr,
              vector<double> &meth_prev, vector<double> &meth_curr,
              double &diff) {

  if (is_leaf(subtree_sizes[node_id])) { //leaf
    if (update_leaves && missing_meth_value(obs_meth_curr[node_id])) {
      const double p_par_orig = meth_curr[parent_id];
      const double p_prev_orig = meth_prev[node_id];

      double p_cur = 0.0;
      for (size_t prev = 0; prev < 2; ++prev) {
        const double p_prev = val_or_flip(p_prev_orig, prev);
        for (size_t par = 0; par < 2; ++par) {
          const double p_par = val_or_flip(p_par_orig, par);
          const double marginal = p_prev*p_par;
          const double p_cur_u = GP[node_id](prev, par, 0);
          const double p_cur_m = GP[node_id](prev, par, 1);
          const double p_cur_mb = p_cur_u/(p_cur_u + p_cur_m);
          p_cur += p_cur_mb*marginal;
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    }
  } else {
    vector<size_t> children;
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      posterior_end(update_leaves, subtree_sizes, G, P, GP, pos,
                    child_id, node_id, obs_meth_curr,
                    meth_prev, meth_curr, diff);
      children.push_back(child_id);
      count += subtree_sizes[child_id];
    }
    if (is_root(node_id)) { //root node  //root has 3 children!
      const double p_prev_orig = meth_prev[node_id];
      const double p_c0_orig = meth_curr[children[0]];
      const double p_c1_orig = meth_curr[children[1]];
      const double p_c2_orig = meth_curr[children[2]];
      const double p_prev_c0_orig = meth_prev[children[0]];
      const double p_prev_c1_orig = meth_prev[children[1]];
      const double p_prev_c2_orig = meth_prev[children[2]];

      double p_cur = 0.0;
      for (size_t prev = 0; prev < 2; ++prev) {
        const double p_prev = val_or_flip(p_prev_orig, prev);
        for (size_t c0 = 0; c0 < 2; ++c0) {
          const double p_c0 = val_or_flip(p_c0_orig, c0);
          for (size_t c1 = 0; c1 < 2; ++c1) {
            const double p_c1 = val_or_flip(p_c1_orig, c1);
            for (size_t c2 = 0; c2 < 2; ++c2) {
              const double p_c2 = val_or_flip(p_c2_orig, c2);
              for (size_t prev_c0 = 0; prev_c0 < 2; ++prev_c0) {
                const double p_prev_c0 = val_or_flip(p_prev_c0_orig, prev_c0);
                for (size_t prev_c1 = 0; prev_c1 < 2; ++prev_c1) {
                  const double p_prev_c1 = val_or_flip(p_prev_c1_orig, prev_c1);
                  for (size_t prev_c2 = 0; prev_c2 < 2; ++prev_c2) {
                    const double p_prev_c2 = val_or_flip(p_prev_c2_orig, prev_c2);
                    const double marginal = (p_prev*p_c0*p_c1*p_c2*
                                             p_prev_c0*p_prev_c1*p_prev_c2);
                    const double p_cur_u = (GP[children[0]](prev_c0, 0, c0)*
                                            GP[children[1]](prev_c1, 0, c1)*
                                            GP[children[2]](prev_c2, 0, c2)*
                                            G(prev, 0));
                    const double p_cur_m = (GP[children[0]](prev_c0, 1, c0)*
                                            GP[children[1]](prev_c1, 1, c1)*
                                            GP[children[2]](prev_c2, 1, c2)*
                                            G(prev, 1));
                    p_cur += p_cur_u/(p_cur_m + p_cur_u)*marginal;
                  }
                }
              }
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    } else { // internal node: has two children
      const double p_c0_orig = meth_curr[children[0]];
      const double p_c1_orig = meth_curr[children[1]];
      const double p_prev_c0_orig = meth_prev[children[0]];
      const double p_prev_c1_orig = meth_prev[children[1]];
      const double p_prev_orig = meth_prev[node_id];
      const double p_par_orig = meth_curr[parent_id];

      double p_cur = 0.0;
      for (size_t c0 = 0; c0 < 2; ++c0) {
        const double p_c0 = val_or_flip(p_c0_orig, c0);
        for (size_t c1 = 0; c1 < 2; ++c1) {
          const double p_c1 = val_or_flip(p_c1_orig, c1);
          for (size_t prev = 0; prev < 2; ++prev) {
            const double p_prev = val_or_flip(p_prev_orig, prev);
            for (size_t par = 0; par < 2; ++par) {
              const double p_par = val_or_flip(p_par_orig, par);
              for (size_t prev_c0 = 0; prev_c0 < 2; ++prev_c0) {
                const double p_prev_c0 = val_or_flip(p_prev_c0_orig, prev_c0);
                for (size_t prev_c1 = 0; prev_c1 < 2; ++prev_c1) {
                  const double p_prev_c1 = val_or_flip(p_prev_c1_orig, prev_c1);
                  const double marginal = (p_c0*p_c1*p_prev*p_par*
                                           p_prev_c0*p_prev_c1);
                  const double p_cur_u = (GP[node_id](prev, par, 0)*
                                          GP[children[0]](prev_c0, 0, c0)*
                                          GP[children[1]](prev_c1, 0, c1));
                  const double p_cur_m = (GP[node_id](prev, par, 1)*
                                          GP[children[0]](prev_c0, 1, c0)*
                                          GP[children[1]](prev_c1, 1, c1));
                  p_cur += p_cur_u/(p_cur_u + p_cur_m)*marginal;
                }
              }
            }
          }
        }
      }
      diff += abs(meth_curr[node_id] - p_cur);
      meth_curr[node_id] = p_cur;
    }
  }
}


static void
approx_posterior(const vector<size_t> &subtree_sizes,
                 const param_set &ps,
                 const vector<size_t> &reset_points,
                 const double tolerance,
                 const size_t MAXITER,
                 const vector<vector<double> > &obs_hypo_prob_table,
                 vector<vector<double> > &hypo_prob_table) {

  // whether to update posterior at leaves with missing observation
  static const bool update_leaves = false;

  // collect probabilities and derivatives by branch
  pair_state G(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);
  vector<pair_state> P;
  vector<triple_state> GP;
  get_transition_matrices(ps, P, GP);

  const size_t n_nodes = subtree_sizes.size();

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const size_t start = reset_points[i];
    const size_t end = reset_points[i + 1] - 1;
    if (start < end) {
      const double tol = (update_leaves ? // ADS: is it possible for this to be true?
                          tolerance*n_nodes*(end - start) :
                          tolerance*n_nodes*(end - start)/2);

      double diff = std::numeric_limits<double>::max();
      for (size_t iter = 0; iter < MAXITER && diff > tol; ++iter) {
        diff = 0.0; // ADS: why set this value here??
        posterior_start(update_leaves, subtree_sizes, G, P, GP, start, 0, 0,
                        obs_hypo_prob_table[start],
                        hypo_prob_table[start], hypo_prob_table[start + 1], diff);
        for (size_t j = start + 1; j < end; ++j)
          posterior_middle(update_leaves, subtree_sizes, G, P, GP, j, 0, 0,
                           obs_hypo_prob_table[j],
                           hypo_prob_table[j - 1], hypo_prob_table[j],
                           hypo_prob_table[j + 1], diff);
        posterior_end(update_leaves, subtree_sizes, G, P, GP, end, 0, 0,
                      obs_hypo_prob_table[end],
                      hypo_prob_table[end - 1], hypo_prob_table[end], diff);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////   Markov blanket probabilities    ///////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static double
MB_middle_site(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
               const pair_state &logG, const vector<triple_state> &logGP,
               const vector<size_t> &meth_prev,
               const vector<size_t> &meth_curr,
               const vector<size_t> &meth_next, const size_t node_id) {
  const size_t parent = parent_ids[node_id];
  double lp0 = 0.0, lp1 = 0.0;

  //parents
  if (is_root(node_id)) {
    const size_t prev_state = meth_prev[node_id];
    lp0 += logG(prev_state, 0);
    lp1 += logG(prev_state, 1);
  } else { // (node_id > 0)
    const size_t parent_state = meth_curr[parent];
    const size_t prev_state = meth_prev[node_id];
    lp0 += logGP[node_id](prev_state, parent_state, 0);
    lp1 += logGP[node_id](prev_state, parent_state, 1);
  }

  // children and children's parents in graph
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const size_t prev_state = meth_prev[child_id];
    const size_t current_state = meth_curr[child_id];
    lp0 += logGP[child_id](prev_state, 0, current_state);
    lp1 += logGP[child_id](prev_state, 1, current_state);
    count += subtree_sizes[child_id];
  }

  if (is_root(node_id)) {
    const size_t current_state = meth_next[node_id];
    lp0 += logG(0, current_state);
    lp1 += logG(1, current_state);
  } else { // (node_id > 0)
    const size_t parent_state = meth_next[parent];
    const size_t current_state = meth_next[node_id];
    lp0 += logGP[node_id](0, parent_state, current_state);
    lp1 += logGP[node_id](1, parent_state, current_state);
  }
  return lp0 - lp1;
}


// ADS: why is pi0 only used for start nodes?
static double
MB_end_site(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
            const pair_state &logG, const vector<triple_state> &logGP,
            const vector<size_t> &meth_prev, const vector<size_t> &meth_curr,
            const size_t node_id) {
  double lp0 = 0.0, lp1 = 0.0;
  // parents
  if (is_root(node_id)) { // is_root(node_id)
    const size_t prev_state = meth_prev[node_id];
    lp0 += logG(prev_state, 0);
    lp1 += logG(prev_state, 1);
  } else { // (node_id > 0)
    const size_t parent_state = meth_curr[parent_ids[node_id]];
    const size_t prev_state = meth_prev[node_id];
    lp0 += logGP[node_id](prev_state, parent_state, 0);
    lp1 += logGP[node_id](prev_state, parent_state, 1);
  }

  // process children
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const size_t prev_child_state = meth_prev[child_id];
    const size_t curr_child_state = meth_curr[child_id];
    lp0 += logGP[child_id](prev_child_state, 0, curr_child_state);
    lp1 += logGP[child_id](prev_child_state, 1, curr_child_state);
    count += subtree_sizes[child_id];
  }
  // ADS: what happens for leaf nodes?

  return lp0 - lp1;
}


static double
MB_start_site(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
              const double pi0, const pair_state &logG,
              const vector<pair_state> &logP, const vector<triple_state> &logGP,
              const vector<size_t> &meth_curr, const vector<size_t> &meth_next,
              const size_t node_id) {

  double lp0 = 0.0, lp1 = 0.0;
  const size_t parent = parent_ids[node_id];
  // transition from parents (or root distribution)
  if (is_root(node_id)) { // is_root(node_id)
    lp0 += log(pi0);
    lp1 += log(1.0 - pi0);
  } else { // (node_id > 0)
    const size_t parent_state = meth_curr[parent];
    lp0 += logP[node_id](parent_state, 0);
    lp1 += logP[node_id](parent_state, 1);
  }

  // iterate over children and process transitions to them
  for (size_t count = 1; count < subtree_sizes[node_id];) {
    const size_t child_id = node_id + count;
    const size_t current_state = meth_curr[child_id];
    lp0 += logP[child_id](0, current_state);
    lp1 += logP[child_id](1, current_state);
    count += subtree_sizes[child_id];
  }

  // transitions to next node horizontally
  if (is_root(node_id)) {
    const size_t next_state = meth_next[node_id];
    lp0 += logG(0, next_state);
    lp1 += logG(1, next_state);
  } else { // non-root node
    const size_t parent_state = meth_next[parent];
    const size_t current_state = meth_next[node_id];
    lp0 += logGP[node_id](0, parent_state, current_state);
    lp1 += logGP[node_id](1, parent_state, current_state);
  }
  return lp0 - lp1;
}


template <class T> static void
accept_or_reject_proposal(const double log_ratio, T &state) {
  const T tmp_state = state;
  const double ratio = (state == 0) ? exp(-log_ratio) : exp(log_ratio);
  std::uniform_real_distribution<> dis(0, 1);  // interval [0, 1)
  if (dis(gen) < ratio)
    state = 1 - tmp_state;
}


static void
MH_update_site_middle(const vector<size_t> &subtree_sizes,
                      const vector<size_t> &parent_ids,
                      const pair_state &logG,
                      const vector<triple_state> &logGP,
                      const vector<double> &probs_table,
                      const vector<size_t> &states_prev,
                      vector<size_t> &states_curr,
                      const vector<size_t> &states_next,
                      const size_t node_id) {

  /* if leaf observed, sample from the observed probability*/
  if (is_leaf(subtree_sizes[node_id]) && probs_table[node_id] >= 0.0) {
    std::uniform_real_distribution<> dis(0, 1);
    states_curr[node_id] = dis(gen) < probs_table[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      MH_update_site_middle(subtree_sizes, parent_ids, logG, logGP,
                            probs_table, states_prev, states_curr,
                            states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio = MB_middle_site(subtree_sizes, parent_ids, logG, logGP,
                                            states_prev, states_curr, states_next,
                                            node_id);

    accept_or_reject_proposal(log_ratio, states_curr[node_id]);
  }
}


static void
MH_update_site_end(const vector<size_t> &subtree_sizes,
                   const vector<size_t> &parent_ids,
                   const pair_state &logG,
                   const vector<triple_state> &logGP,
                   const vector<double> &probs_table,
                   const vector<size_t> &states_prev,
                   vector<size_t> &states_curr,
                   const size_t node_id) {

  /* if leaf observed, sample from the observed probability*/
  if (is_leaf(subtree_sizes[node_id]) && probs_table[node_id] >= 0) {
    std::uniform_real_distribution<> dis(0, 1);
    states_curr[node_id] = dis(gen) < probs_table[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id];) {
      const size_t child_id = node_id + count;
      MH_update_site_end(subtree_sizes, parent_ids, logG, logGP,
                         probs_table, states_prev, states_curr, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio = MB_end_site(subtree_sizes, parent_ids, logG, logGP,
                                         states_prev, states_curr, node_id);

    accept_or_reject_proposal(log_ratio, states_curr[node_id]);
  }
}


static void
MH_update_site_start(const vector<size_t> &subtree_sizes,
                     const vector<size_t> &parent_ids,
                     const double pi0,
                     const pair_state &logG,
                     const vector<pair_state> &logP,
                     const vector<triple_state> &logGP,
                     const vector<double> &probs_table,
                     vector<size_t> &states_curr,
                     const vector<size_t> &states_next,
                     const size_t node_id) {

  /* if leaf observed, sample from the observed probability*/
  if (is_leaf(subtree_sizes[node_id]) && probs_table[node_id] >= 0) {
    std::uniform_real_distribution<> dis(0, 1);
    states_curr[node_id] = dis(gen) < probs_table[node_id] ? 0 : 1;
  } else {
    for (size_t count = 1; count < subtree_sizes[node_id]; ) {
      const size_t child_id = node_id + count;
      MH_update_site_start(subtree_sizes, parent_ids, pi0, logG, logP, logGP,
                           probs_table, states_curr, states_next, child_id);
      count += subtree_sizes[child_id];
    }
    const double log_ratio =
      MB_start_site(subtree_sizes, parent_ids, pi0, logG, logP, logGP,
                    states_curr, states_next, node_id);

    accept_or_reject_proposal(log_ratio, states_curr[node_id]);
  }
}


static void
MH_update(const vector<size_t> &subtree_sizes, const vector<size_t> &parent_ids,
          const param_set &ps, const vector<size_t> &reset_points,
          const vector<vector<double> > &probs_table,
          vector<vector<size_t> > &states_table) {

  pair_state logG(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);
  logG.make_logs();

  vector<pair_state> logP;
  vector<triple_state> logGP;
  get_transition_matrices(ps, logP, logGP);
  for (size_t i = 0; i < logP.size(); ++i) {
    logP[i].make_logs();
    logGP[i].make_logs();
  }

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const size_t start = reset_points[i];
    const size_t end = reset_points[i + 1] - 1; // !!!!
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
             const vector<vector<size_t> > &tree_state_table,
             const vector<size_t> &reset_points,
             vector<triple_state> &triad_counts,
             vector<pair_state> &start_counts,
             pair_state &root_counts) {

  triad_counts = vector<triple_state>(subtree_sizes.size());
  start_counts = vector<pair_state>(subtree_sizes.size());
  root_counts = pair_state();

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {

    const size_t start = reset_points[i];
    const size_t end = reset_points[i + 1];

    for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
      const size_t anc = tree_state_table[start][parent_ids[node_id]];
      const size_t cur = tree_state_table[start][node_id];
      start_counts[node_id](anc, cur) += 1.0;
    }

    for (size_t pos = start + 1; pos < end; ++pos) {

      root_counts(tree_state_table[pos - 1][0], tree_state_table[pos][0]) += 1.0;

      for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
        const size_t anc = tree_state_table[pos][parent_ids[node_id]];
        const size_t prev = tree_state_table[pos - 1][node_id];
        const size_t cur = tree_state_table[pos][node_id];
        triad_counts[node_id](prev, anc, cur) += 1.0;
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////         OPTIMIZATION          /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static double
log_likelihood(const vector<size_t> &subtree_sizes, const param_set &ps,
               const vector<triple_state> &triad_weights, // treesize x 2 x 2 x 2
               const vector<pair_state > &start_weights,// treesize x 2 x 2
               const pair_state &root_weights) {

  vector<pair_state> P;
  vector<triple_state> GP;
  get_transition_matrices(ps, P, GP);

  double llk =
    (root_weights(0, 0)*log(ps.g0) + root_weights(0, 1)*log(1.0 - ps.g0) +
     root_weights(1, 0)*log(1.0 - ps.g1) + root_weights(1, 1)*log(ps.g1));

  for (size_t node = 1; node < subtree_sizes.size(); ++node)
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        llk += start_weights[node](j, k)*log(P[node](j, k));
        for (size_t i = 0; i < 2; ++i)
          llk += triad_weights[node](i, j, k)*log(GP[node](i, j, k));
      }

  return llk;
}

void
objective_branch(const param_set &ps,
                 const vector<triple_state> &triad_weights,
                 const vector<pair_state> &start_weights,
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
        F += triad_weights[node_id](i, j, k)*log(GP[node_id](i, j, k));
        deriv += triad_weights[node_id](i, j, k)*
          (GP_dT[node_id](i, j, k)/GP[node_id](i, j, k));
      }

  double rate0 = ps.rate0;
  pair_state P_dT(-rate0, rate0, 1.0 - rate0, rate0 - 1.0);

  for (size_t j = 0; j < 2; ++j)
    for (size_t k = 0; k < 2; ++k) {
      F += start_weights[node_id](j, k)*log(P[node_id](j, k));
      deriv += start_weights[node_id](j, k)*(P_dT(j, k)/P[node_id](j, k));
    }
}


template <class T> static bool
btwn_01_inclusive(const T x, const double &epsilon) {
  return x > (0.0 + epsilon) && x < (1.0 - epsilon);
}


static void
update_branch(const bool VERBOSE, const vector<size_t> &subtree_sizes,
              const param_set &orig_ps, size_t node_id,
              const vector<triple_state> &triad_weights,
              const vector<pair_state> &start_weights, param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(orig_ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  param_set prev_ps(orig_ps);
  double prev_F, new_F;
  double prev_deriv, new_deriv;

  objective_branch(prev_ps, triad_weights, start_weights, P, GP, GP_dT,
                   node_id, prev_F, prev_deriv);

  double frac = 1.0;
  bool CONVERGE = false;
  double improvement = 0.0;
  while (!CONVERGE) {
    while (frac > param_set::tolerance &&
           !btwn_01_inclusive(frac*sign(prev_deriv) +
                              prev_ps.T[node_id], param_set::tolerance)) {
      frac = frac/2;
    }

    ps.T[node_id] = prev_ps.T[node_id] + frac*sign(prev_deriv);

    // update trans mats and derivs with new-params
    get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

    objective_branch(ps, triad_weights, start_weights, P, GP,
                     GP_dT, node_id, new_F, new_deriv);

    if (new_F > prev_F) {
      improvement += new_F - prev_F;
      if (VERBOSE)
        cerr << "[update_branch: delta=" << improvement
             << ", branch[" << node_id << "]=" << ps.T[node_id] << ']' << endl;
      prev_F = new_F;
      prev_deriv = new_deriv;
      prev_ps.T[node_id] = ps.T[node_id];
    }
    frac = frac/2;
    if (frac < param_set::tolerance)
      CONVERGE = true;
  }
  ps.T[node_id] = prev_ps.T[node_id];
}


static void
objective_rate(const vector<size_t> &subtree_sizes,
               const param_set &ps, const vector<triple_state> &triad_weights,
               const vector<pair_state> &start_weights,
               const vector<pair_state> &P, const vector<triple_state> &GP,
               const vector<triple_state> &GP_drate,
               double &F, double &deriv_rate) {

  const size_t n_nodes = subtree_sizes.size();

  F = 0.0;
  deriv_rate = 0.0;
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k) {
          F += triad_weights[node_id](i, j, k)*log(GP[node_id](i, j, k));
          deriv_rate += triad_weights[node_id](i, j, k)*
            (GP_drate[node_id](i, j, k)/GP[node_id](i, j, k));
        }

    const double T_val = ps.T[node_id];
    const pair_state P_drate(-T_val, T_val, -T_val, T_val);
    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        F += start_weights[node_id](j, k)*log(P[node_id](j, k));
        deriv_rate += start_weights[node_id](j, k)*
          P_drate(j, k)/P[node_id](j, k);
      }
  }
}


static void
update_rate(const bool VERBOSE, const vector<size_t> &subtree_sizes,
            const param_set &orig_ps, const vector<triple_state> &triad_weights,
            const vector<pair_state> &start_weights, param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(orig_ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

  param_set prev_ps = orig_ps;
  ps = orig_ps;

  double prev_F = 0.0;
  double new_F = 0.0;
  double prev_deriv, new_deriv;
  objective_rate(subtree_sizes, prev_ps, triad_weights,
                 start_weights, P, GP, GP_drate, prev_F, prev_deriv);

  double denom = abs(prev_deriv);
  double frac = 1.0;
  bool CONVERGE = false;
  double improvement = 0.0;
  while (!CONVERGE) {
    while (frac > param_set::tolerance &&
           !btwn_01_inclusive(frac*sign(prev_deriv) +
                              prev_ps.rate0, param_set::tolerance)) {
      frac = frac/2;
    }
    ps.rate0 = prev_ps.rate0 + frac*(prev_deriv/denom);
    get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);

    objective_rate(subtree_sizes, ps, triad_weights, start_weights,
                   P, GP, GP_drate, new_F, new_deriv);

    if (new_F > prev_F) {
      improvement += new_F - prev_F;
      cerr << "[update_rate: delta=" << improvement
           << ", rate=" << ps.rate0 << ']' << endl;
      prev_F = new_F;
      prev_deriv = new_deriv;
      prev_ps.rate0 = ps.rate0;
      denom = abs(new_deriv);
    }

    frac = frac/2;
    if (frac < param_set::tolerance)
      CONVERGE = true;
  }
  ps.rate0 = prev_ps.rate0;
}


static void
objective_G(const vector<size_t> &subtree_sizes, const param_set &ps,
            const vector<triple_state> &triad_weights,
            const pair_state &root_weights,
            const vector<triple_state> &GP,
            const vector<triple_state> &GP_dg0,
            const vector<triple_state> &GP_dg1,
            double &F, vector<double> &deriv_G) {

  const size_t n_nodes = subtree_sizes.size();

  F = 0.0;
  deriv_G = vector<double>(2, 0.0);
  for (size_t node_id = 1; node_id < n_nodes; ++node_id) {

    for (size_t i = 0; i < 2; ++i)
      for (size_t j = 0; j < 2; ++j)
        for (size_t k = 0; k < 2; ++k)
          F += triad_weights[node_id](i, j, k)*log(GP[node_id](i, j, k));


    for (size_t j = 0; j < 2; ++j)
      for (size_t k = 0; k < 2; ++k) {
        deriv_G[0] += triad_weights[node_id](0, j, k)*
          (GP_dg0[node_id](0, j, k)/GP[node_id](0, j, k));
        deriv_G[1] += triad_weights[node_id](1, j, k)*
          (GP_dg1[node_id](1, j, k)/GP[node_id](1, j, k));
      }
  }

  const pair_state G(ps.g0, 1.0 - ps.g0, 1.0 - ps.g1, ps.g1);

  for (size_t i = 0; i < 2; ++i)
    for (size_t j = 0; j < 2; ++j)
      F += root_weights(i, j)*log(G(i, j));

  deriv_G[0] += root_weights(0, 0)/G(0, 0) - root_weights(0, 1)/G(0, 1);
  deriv_G[1] += -1.0*root_weights(1, 0)/G(1, 0) + root_weights(1, 1)/G(1, 1);
}


static void
update_G(const bool VERBOSE, const vector<size_t> &subtree_sizes,
         const param_set &orig_ps, const vector<triple_state> &triad_weights,
         const pair_state &root_weights, param_set &ps) {

  vector<pair_state> P;
  vector<triple_state> GP, GP_drate, GP_dg0, GP_dg1, GP_dT;
  get_transition_matrices_deriv(orig_ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);
  param_set prev_ps(orig_ps);
  ps = orig_ps;

  double prev_F = 0.0, new_F = 0.0; // ADS: naming is bad (never use "new")
  vector<double> prev_deriv, new_deriv;

  objective_G(subtree_sizes, prev_ps, triad_weights, root_weights,
              GP, GP_dg0, GP_dg1, prev_F, prev_deriv);

  double denom = abs(prev_deriv[0]) + abs(prev_deriv[1]);
  double frac = 1.0;
  bool converged = false;
  double improvement = 0.0;
  while (!converged) {
    while (frac > param_set::tolerance &&
           !(btwn_01_inclusive(frac*(prev_deriv[0]/denom) + prev_ps.g0, param_set::tolerance) &&
             btwn_01_inclusive(frac*(prev_deriv[1]/denom) + prev_ps.g1, param_set::tolerance))) {
      frac = frac/2;
    }
    ps.g0 = prev_ps.g0 + frac*(prev_deriv[0]/denom); //g0
    ps.g1 = prev_ps.g1 + frac*(prev_deriv[1]/denom); //g1

    get_transition_matrices_deriv(ps, P, GP, GP_drate, GP_dg0, GP_dg1, GP_dT);
    objective_G(subtree_sizes, ps, triad_weights, root_weights,
                GP, GP_dg0, GP_dg1, new_F, new_deriv);

    if (new_F > prev_F) {
      improvement += new_F - prev_F;
      if (VERBOSE)
        cerr << "[update_G: delta=" << improvement
             << ", G=(" << ps.g0 << ", " << ps.g1 << ")]" << endl;
      prev_F = new_F;
      prev_deriv = new_deriv;
      prev_ps.g0 = ps.g0;
      prev_ps.g1 = ps.g1;
      denom = abs(new_deriv[0]) + abs(new_deriv[1]);
    }
    frac = frac/2;
    if (frac < param_set::tolerance)
      converged = true;
  }
  ps.g0 = prev_ps.g0;
  ps.g1 = prev_ps.g1;
}


static void
update_pi0(const bool VERBOSE,
           const vector<pair_state> &start_weights, param_set &ps) {
  const double num = start_weights[1](0, 0) + start_weights[1](0, 1);
  const double denom = num + start_weights[1](1, 0) + start_weights[1](1, 1);
  ps.pi0 = num/denom;
  if (VERBOSE)
    cerr << "[update_pi0: pi0=" << ps.pi0 << ']' << endl;
}


static void
optimize_params(const bool VERBOSE, const vector<size_t> &subtree_sizes,
                const vector<size_t> &reset_points,
                const pair_state &root_weights,
                const vector<pair_state> &start_weights,
                const vector<triple_state> &triad_weights, param_set &ps) {

  // ADS: not clear why to keep these two variables below
  param_set cur_ps(ps);

  // update branches
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++node_id) {
    update_branch(VERBOSE, subtree_sizes, cur_ps, node_id,
                  triad_weights, start_weights, ps);
    cur_ps.T[node_id] = ps.T[node_id];
  }

  // update rate0
  update_rate(VERBOSE, subtree_sizes, cur_ps, triad_weights, start_weights, ps);
  cur_ps.rate0 = ps.rate0;

  //update pi0
  update_pi0(VERBOSE, start_weights, ps);
  cur_ps.pi0 = ps.pi0;

  //update G
  update_G(VERBOSE, subtree_sizes, cur_ps, triad_weights, root_weights, ps);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////// OPTIMIZING PARAMETERS ABOVE //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////// EXPECTATION MAXIMIZATION ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
maximization_step(const bool VERBOSE, const size_t MAXITER,
                  const vector<size_t> &subtree_sizes,
                  const vector<size_t> &reset_points,
                  const pair_state &root_weights,
                  const vector<pair_state> &start_weights,
                  const vector<triple_state> &triad_weights,
                  param_set &ps, PhyloTreePreorder &t) {

  // one M-step: optimize parameters
  double diff = std::numeric_limits<double>::max();
  for (size_t iter = 0; iter < MAXITER && diff > param_set::tolerance; ++iter) {
    if (VERBOSE)
      cerr << "[inside maximization: iter=" << iter << "]\n" << ps << endl;

    param_set prev_ps(ps);
    optimize_params(VERBOSE, subtree_sizes, reset_points,
                    root_weights, start_weights, triad_weights, ps);
    diff = param_set::absolute_difference(ps, prev_ps);
  }

  double llk = log_likelihood(subtree_sizes, ps, triad_weights,
                              start_weights, root_weights);
  if (VERBOSE)
    cerr << "log-likelihood=" << llk << endl;

}


static void
to_discrete_table(const vector<vector<double> > &tree_prob_table,
                  vector<vector<size_t> > &tree_prob_table_discrete) {
  const bool randomize = true;
  for (size_t i = 0; i < tree_prob_table.size(); ++i)
    for (size_t j = 0; j < tree_prob_table[i].size(); ++j) {
      if (randomize) {
        std::uniform_real_distribution<> dis(0, 1);  // interval [0, 1)
        tree_prob_table_discrete[i][j] = (dis(gen) < tree_prob_table[i][j]) ? 0: 1;
      } else {
        // 0==hypo state; 1==methylated state (opposite of the hypo prob)
        tree_prob_table_discrete[i][j] = (tree_prob_table[i][j] > 0.5)? 0: 1;
      }
    }
}

static void
expectation_step(const bool VERBOSE, const size_t mh_max_iterations,
                 const size_t max_app_iter,
                 const vector<size_t> &subtree_sizes,
                 const vector<size_t> &parent_ids,
                 const vector<size_t> &reset_points,
                 const param_set &ps,
                 pair_state &root_weights,
                 vector<pair_state> &start_weights,
                 vector<triple_state> &triad_weights,
                 vector<vector<double> > &tree_prob_table,
                 vector<vector<size_t> > &sampled_states) {

  const size_t n_nodes = subtree_sizes.size();

  triad_weights = vector<triple_state>(n_nodes);
  start_weights = vector<pair_state>(n_nodes);
  vector<triple_state> triad_weights_prev(n_nodes);

  bool converged = false;
  size_t mh_iter = 0;
  for (mh_iter = 0; mh_iter < mh_max_iterations && !converged; ++mh_iter) {
    if (VERBOSE)
      cerr << "[inside expectation: M-H (iter=" << mh_iter << ")]" << endl;

    pair_state root_weights_samp;
    vector<pair_state> start_weights_samp(n_nodes);
    vector<triple_state> triad_weights_samp(n_nodes);
    count_triads(subtree_sizes, parent_ids, sampled_states, reset_points,
                 triad_weights_samp, start_weights_samp, root_weights_samp);

    for (size_t i = 0; i < triad_weights.size(); ++i) {
      triad_weights[i] = triad_weights[i] + triad_weights_samp[i];
      start_weights[i] = start_weights[i] + start_weights_samp[i];
    }
    root_weights = root_weights + root_weights_samp;

    vector<double> divergence;
    kl_divergence(triad_weights_prev, triad_weights, divergence);

    converged =
      *max_element(divergence.begin(), divergence.end()) < KL_CONV_TOL;
    if (!converged) { // take next sample (all sites)
      /* "tree_prob_table" determines how leaf states are updated */
      MH_update(subtree_sizes, parent_ids, ps, reset_points,
                tree_prob_table, sampled_states);
    }
    triad_weights_prev = triad_weights;
  }

  root_weights.div(mh_iter);
  for (size_t i = 0; i < triad_weights.size(); ++i) {
    start_weights[i].div(mh_iter);
    triad_weights[i].div(mh_iter);
  }
  if (VERBOSE) {
    cerr << "[MH iterations=" << mh_iter << ']' << endl
         << "triad_weights:" << endl;
    copy(triad_weights.begin(), triad_weights.end(),
         ostream_iterator<triple_state>(cerr, "\n"));
    // ADS: should we print more summary statistics here?
  }
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
        break;
      }
      count += subtree_sizes[child_id];
    }
    if (SAME) {
      state[node_id] = s;
    } else {
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



static void
write_treeprob_states(const vector<size_t> &subtree_sizes,
                      const vector<Site> &sites,
                      const vector<vector<double> > &tree_prob_table,
                      std::ofstream &out) {

  const size_t n_nodes = subtree_sizes.size();
  for (size_t i = 0; i < sites.size(); ++i) {
    string state;
    string pme;
    for (size_t j = 0; j < n_nodes; ++j ) {
      state += (tree_prob_table[i][j] <= 0.5 ? '1' : '0');
      if (is_leaf(subtree_sizes[j]))
        pme += (tree_prob_table[i][j] <= 0.5 ? '1' : '0');
    }

    out << sites[i].chrom << "\t" << sites[i].pos << "\t"
        << sites[i].pos + 1 << "\t" << state << "\t0\t+\t";
    out << pme;
    for (size_t j = 0; j < n_nodes; ++j)
      out << "\t" << tree_prob_table[i][j];
    out << "\n";
  }
}


static void
build_domain(const size_t minCpG, const size_t desert_size,
             const vector<Site> &sites, const vector<string> &states,
             vector<GenomicRegion> &domains) {

  // first round collapsing
  GenomicRegion cpg(sites[0].chrom, sites[0].pos, sites[0].pos + 1,
                    states[0], 1.0, '+');
  for (size_t i = 1; i < sites.size(); ++i) {
    const size_t d = Site::distance(sites[i], sites[i - 1]);
    if (d < desert_size && states[i] == states[i - 1]) {
      cpg.set_score(1.0 + cpg.get_score());
      cpg.set_end(sites[i].pos + 1);
    } else {
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


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////          MAIN             ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {
  try{
    size_t desert_size = 1000;
    size_t minCpG = 10;
    size_t min_sites_per_frag = 5; // ADS: probably should be a post-filter step
    size_t MAXITER = 10;           // iterations inside the M-step
    size_t EMMAXITER = 30;         //rounds of EM iterations

    string outfile;
    string paramfile, outparamfile;

    // run mode flags
    bool VERBOSE = false;
    bool SINGLE = false;
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
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("min-fragsize", 'f', "min CpG in fragments to output "
                      "(default: " + toa(min_sites_per_frag) + ")", false,
                      min_sites_per_frag);
    opt_parse.add_opt("single", 's', "also output states by sites "
                      "(when -o is used)", false, SINGLE);

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
    vector<string> nodenames;
    vector<double> branches;
    t.get_subtree_sizes(subtree_sizes);
    t.get_branch_lengths(branches);
    get_degrees(subtree_sizes, node_degrees);

    if (!is_semi_binary(node_degrees))
      throw SMITHLABException("invalid tree structure");
    const size_t n_leaves = count_leaves(subtree_sizes);

    vector<size_t> parent_ids;
    get_parent_id(subtree_sizes, parent_ids);

    if (VERBOSE)
      cerr << "loaded tree (leaves=" << n_leaves << ")" << endl;

    /******************** INITIALIZE PARAMETERS *******************************/
    double pi0 = 0.48; // MAGIC
    double rate0 = 0.4;// MAGIC
    double g0 = 0.9;   // MAGIC
    double g1 = 0.95;  // MAGIC

    param_set init_ps(pi0, rate0, g0, g1);

    if (!paramfile.empty()) {
      init_ps.read(paramfile, t);
      if (VERBOSE)
        cerr << "[tree:]\n" << t.tostring() << endl
             << "[given params:]\n" << init_ps << endl;
    } else {
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

    /******************** LOAD METHYLATION DATA *******************************/
    vector<Site> sites;
    vector<size_t> reset_points;
    vector<vector<double> > tree_prob_table;
    vector<vector<double> > tree_state_table;

    if (COMPLETE) {
      if (VERBOSE)
        cerr << "============ COMPLETE MODE =============" << endl
             << "[Loading full table]" << endl;

      load_full_states_as_prob(meth_table_file, n_nodes, sites, tree_prob_table);

      if (VERBOSE)
        cerr << "[Separating deserts]" << endl;
      separate_regions(desert_size, sites, reset_points);
      const size_t n_sites = tree_prob_table.size();

      // true weights
      pair_state root_weights;
      vector<pair_state> start_weights;
      vector<triple_state> triad_weights;

      vector<vector<size_t> > tree_state_table(n_sites,
                                               vector<size_t>(n_nodes, 0));
      to_discrete_table(tree_prob_table, tree_state_table);
      count_triads(subtree_sizes, parent_ids, tree_state_table, reset_points,
                   triad_weights, start_weights, root_weights);

      if (VERBOSE) {
        cerr << "-----triad_weights_true------" << endl;
        for (size_t i = 0; i < triad_weights.size(); ++i)
          cerr << triad_weights[i];
      }

      // One M-step: optimize parameters
      /*M-step: optimize parameters */
      maximization_step(VERBOSE, MAXITER, subtree_sizes, reset_points,
                        root_weights, start_weights, triad_weights, init_ps, t);



    } else {
      if (VERBOSE)
        cerr << endl << "========= LEAF (INCOMPLETE DATA) MODE =========" << endl
             << "[Loading LEAF table]" << endl;

      vector<string> meth_table_species;
      vector<vector<double> > meth_prob_table;
      load_meth_table(meth_table_file, sites, meth_prob_table, meth_table_species);

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

      // separate by deserts
      vector<vector<double> > meth_prob_full_table;
      vector<size_t> reset_points;
      double g0_est, g1_est, pi0_est;
      /* meth_prob_full_table has missing data filled by a heuristic */
      separate_regions(VERBOSE, desert_size, minCpG,
                       meth_prob_table, meth_prob_full_table,
                       sites, reset_points, g0_est, g1_est, pi0_est);

      const size_t n_sites = sites.size();

      if (VERBOSE)
        cerr << "[post-filter=" << n_sites << "; blocks="
             << reset_points.size() - 1 << "]" << endl;

      /* heuristic estimates for a few parameters*/
      if (paramfile.empty()) {
        init_ps.pi0 = pi0_est;
        init_ps.g0 = g0_est;
        init_ps.g1 = g1_est;
      }

      /* meth_prob_table and tree_prob_table have missing data assigned to -1 */
      tree_prob_table = vector<vector<double> >(sites.size(),
                                                vector<double>(n_nodes, -1.0));
      copy_leaf_to_tree_prob(subtree_sizes, meth_prob_table, tree_prob_table);

      /* tree_prob_table_full have missing data filled*/
      vector<vector<double> > tree_prob_table_full(sites.size(),
                                                   vector<double>(n_nodes, -1.0));
      leaf_to_tree_prob(subtree_sizes, meth_prob_full_table, tree_prob_table_full);

      pair_state root_weights;
      vector<pair_state> start_weights;
      vector<triple_state> triad_weights;

      bool em_converged = !paramfile.empty(); // params fixed => converged
      const size_t max_app_iter = 100;        //MAGIC
      const size_t mh_max_iterations = 500;   //MAGIC

      // Initialize "sample" states from the heuristic probabilities
      vector<vector<size_t> > tree_state_table(n_sites, vector<size_t>(n_nodes));
      to_discrete_table(tree_prob_table_full, tree_state_table);

      param_set curr_ps(init_ps);
      for (size_t iter = 0; iter < EMMAXITER && !em_converged; ++iter) {

        if (VERBOSE)
          cerr << endl << "====================[EM ITERATION=" << iter
               << "]=======================" << endl;

        const param_set prev_ps(curr_ps);
        // next line is optional, just for consistency with older version.
        // tree_prob_table = tree_prob_table_full;
        expectation_step(VERBOSE, mh_max_iterations, max_app_iter,
                         subtree_sizes, parent_ids, reset_points, curr_ps,
                         root_weights, start_weights, triad_weights,
                         tree_prob_table_full, tree_state_table);

        if (VERBOSE) {
          cerr << "[E-step iter=" << iter << "]\ntriad weights:\n" << endl;
          for (size_t i = 0; i < triad_weights.size(); ++i)
            cerr << triad_weights[i] << "\tnode=" << i << endl;
        }

        /****************** M-step: optimize parameters ************************/
        maximization_step(VERBOSE, MAXITER, subtree_sizes, reset_points,
                          root_weights, start_weights, triad_weights, curr_ps, t);
        const double diff = param_set::absolute_difference(prev_ps, curr_ps);

        em_converged = (diff < param_set::tolerance*curr_ps.T.size());

        if (VERBOSE) {
          cerr << "[M-step iter=" << iter << ", params:" << endl
               << curr_ps << endl;
          cerr << "[EM iter=" << iter << ", delta=" << diff
               << ", conv=" << em_converged << ']' << endl;
        }
      } // end E-M iterations

      /* final parameters */
      init_ps = curr_ps;

      /* update phylotree with branch lenghts */
      vector<double> branches(init_ps.T.size());
      for (size_t i = 0; i < init_ps.T.size(); ++i)
        branches[i] = -log(1.0 - init_ps.T[i]);
      t.set_branch_lengths(branches);

      /* get posterior one last time*/
      if (VERBOSE)
        cerr << "[last round of posterior estimation]" << endl;
      approx_posterior(subtree_sizes, init_ps, reset_points, POST_CONV_TOL,
                       max_app_iter, tree_prob_table, tree_prob_table_full);

      /************************** Output states ********************************/
      if (!outfile.empty()) {
        std::ofstream out(outfile.c_str());
        if (!out)
          throw SMITHLABException("bad output file: " + outfile);

        /*build domain*/
        vector<string> states;
        vector<GenomicRegion> domains;
        double cutoff = 0.5; // MAGIC
        tree_prob_to_states(subtree_sizes, parent_ids, tree_prob_table_full,
                            cutoff, states);
        build_domain(min_sites_per_frag, desert_size, sites, states, domains);
        if (VERBOSE)
          cerr << "Built total " << domains.size() << " domains" << endl;

        copy(domains.begin(), domains.end(),
             ostream_iterator<GenomicRegion>(out, "\n"));

        /* output individual sites */
        if (SINGLE) {
          const string outssfile(outfile + "_bysite");
          std::ofstream outss(outssfile.c_str());
          if (!outss)
            throw std::runtime_error("bad output file: " + outssfile);
          write_treeprob_states(subtree_sizes, sites, tree_prob_table_full, outss);
        }
      }
    }

    if (!outparamfile.empty()) {
      std::ofstream out(outparamfile.c_str());
      if (!out)
        throw std::runtime_error("bad output file: " + outparamfile);
      // ADS: check that the tree always has correct branch lengths at
      // this point, even if they are not used
      out << t.Newick_format() << endl << init_ps << endl;
    }
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
