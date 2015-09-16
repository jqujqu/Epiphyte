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
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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
#include <iomanip> // std::setw
#include <fstream>
#include <numeric> //std::accumulate
#include <tr1/unordered_map>
#include <tr1/random>
#include <algorithm>
#include <cmath>    //std::abs, floor
#include <limits>   //std::numeric_limits
#include <iterator>     // std::distance
#include <unistd.h>

/* from smithlab_cpp */
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

/* from methpipe */
#include "MethpipeFiles.hpp"
//#include "BetaBin.hpp"

/* headers for epigenomic evolution */
#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::setw;
using std::pair;
using std::make_pair;
using std::tr1::unordered_map;
using std::abs;
using std::numeric_limits;
using std::accumulate;
using std::inner_product;
using std::min;
using std::max;
using std::istringstream;
using std::istream_iterator;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////        UTILITIES         //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

double
kronecker_delta( const size_t a, const size_t b) {
  return (a==b) ? 1.0 : 0.0;
}

double
sign(double val) {
  return (val >= 0)? 1.0 : -1.0;
}

struct Logscale {
  double logval;
  double symbol;
  Logscale(): logval(0.0), symbol(0.0) {}
  Logscale(const double &l, const double &s) :
    logval(l), symbol(s) {}

  void assign_from_val(const double val) {
    if (val == 0) {
      symbol = 0;
    } else {
      symbol = (val>0) ? 1.0 : -1.0;
      logval = log(abs(val));
    }
  }
};

// val= log[abs(exp(p)*sign_p + exp(q)*sign_q)]
// sign_val = sign(exp(p)*sign_p + exp(q)*sign_q)
static void
log_sum_log_sign(const Logscale p, const Logscale q,
                 Logscale &result) {
  if (p.symbol == 0) {
    result.logval = q.logval;
    result.symbol = q.symbol;
  } else if (q.symbol == 0) {
    result.logval = p.logval;
    result.symbol = p.symbol;
  } else {
    const double larger = (p.logval > q.logval) ? p.logval : q.logval;
    const double smaller = (p.logval > q.logval) ? q.logval : p.logval;
    result.symbol = (p.logval > q.logval)? p.symbol : q.symbol;
    if (p.symbol*q.symbol > 0) {
      result.logval = larger + log(1.0 + exp(smaller - larger));
    } else {
      result.logval = larger + log(1.0 - exp(smaller - larger));
    }
  }
}

// val= log[abs(exp(p)*sign_p + exp(q)*sign_q)]
// sign_val = sign(exp(p)*sign_p + exp(q)*sign_q)
static void
log_sum_log_sign(const double p, const double sign_p,
                 const double q, const double sign_q,
                 double &val, double &sign_val) {
  if (p == 0 || sign_p == 0) {
    val= q;
    sign_val = sign(sign_q);
  } else if (q == 0 || sign_q == 0) {
    val= p;
    sign_val = sign(sign_p);
  } else {
    const double larger = (p > q) ? p : q;
    const double smaller = (p > q) ? q : p;
    sign_val = (p > q)? sign(sign_p) : sign(sign_q);
    if (sign_p*sign_q > 0) {
      val = larger + log(1.0 + exp(smaller - larger));
    } else {
      val = larger + log(1.0 - exp(smaller - larger));
    }
  }
}


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
  Site(const string &chr, const size_t &position) :
    chrom(chr), pos(position) {}

  void
  assign_from_methtab_first_column(const string &line) {
    const size_t chrom_end = line.find_first_of(':');
    chrom = line.substr(0, chrom_end);
    const size_t pos_end = line.find_first_of(':', chrom_end + 1);
    pos = atoi(line.substr(chrom_end + 1, pos_end - chrom_end).c_str());
  }
};

size_t
distance_between_sites(Site s, Site t) {
  size_t d = std::numeric_limits<size_t>::max();
  if (s.chrom == t.chrom)
    d = (t.pos > s.pos)? t.pos - s.pos : s.pos - t.pos;
  return d;
}

void
parse_table_line(istringstream &iss,
                 vector<vector<pair<size_t, size_t> > > &meth) {
  meth.push_back(vector<pair<size_t, size_t> >());
  size_t meth_count = 0, unmeth_count = 0;
  while (iss >> meth_count) {
    if (iss >> unmeth_count) {
      meth.back().push_back(make_pair(meth_count, unmeth_count));
    } else {
      throw SMITHLABException("bad line format: " + iss.str());
    }
  }
}



void
parse_table_line(istringstream &iss,
                 vector<vector<double> > &meth) {
  meth.push_back(vector<double>(istream_iterator<double>(iss),
                                istream_iterator<double>()));
  for (size_t i = 0; i < meth.back().size(); ++ i) {
    if ((meth.back()[i] < 0 && meth.back()[i]!=-1.0) || meth.back()[i] > 1.0 )
      throw SMITHLABException("bad line content: [" + iss.str() +
                              "], expecting hypomethylation probabilities");
  }
  if (meth.back().empty())
    throw SMITHLABException("bad line format: " + iss.str());
}



template <class T> void
parse_meth_table_line(const string &line, vector<Site> &sites,
                      vector<vector<T> > &meth) {

  istringstream iss(line);

  // take care of the site location (chrom and position)
  string first_column;
  iss >> first_column;
  Site site;
  site.assign_from_methtab_first_column(first_column);
  sites.push_back(site);

  // now get the methylation information, either as pairs of read
  // counts or posteriors, depending on the file type (run mode)
  parse_table_line(iss, meth);
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
read_params(const bool VERBOSE, const string &paramfile,
            double &root_unmeth_prob, double &lam,
            PhyloTreePreorder &t) {
  std::ifstream in(paramfile.c_str());
  if (!in)
    throw SMITHLABException("cannot read: " + paramfile);

  string line;
  getline(in, line);
  istringstream iss(line);
  iss >> t;
  if (VERBOSE)
    cerr << "Read in tree branch lengths:" << endl
         << t.tostring() << endl;

  getline(in, line);
  istringstream issp(line);
  issp >> root_unmeth_prob >> lam;
  if (VERBOSE)
    cerr << "Root_unmeth_prob=" << root_unmeth_prob
         << "\trate=" <<  lam << endl;
  if (root_unmeth_prob >= 1.0 || root_unmeth_prob <= 0.0 ||
      lam >= 1.0 || lam <= 0.0)
    throw SMITHLABException("wrong parameter range: " + paramfile);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////// ABOVE CODE IS FOR LOADING DATA  //////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static bool
has_same_species_order(const PhyloTreePreorder &the_tree,
                       const vector<string> &meth_table_species) {
  vector<string> leaf_names;
  the_tree.get_leaf_names(leaf_names);
  return leaf_names == meth_table_species;
}

static void
separate_regions(const size_t desert_size,
                 const size_t min_site_per_block,
                 vector<vector<double> > &meth, vector<Site> &sites,
                 vector<size_t> &reset_points) {
  const size_t nleaf = meth[0].size();
  const size_t totalsites = sites.size();

  //scan desert
  //uncovered site has prob value set to -1.0
  vector<bool> is_desert(totalsites, false);
  for (size_t i = 0; i < nleaf; ++i) {
    cerr << "sample" << i ;
    size_t j = 0;
    for (j = 0; j < totalsites && meth[j][i] == -1.0; ++j) {
      is_desert[j] = true;
    }

    if (j == totalsites)
      throw SMITHLABException("No valid observation.");

    Site prev_obs = sites[j];
    for (size_t k = j+1; k < totalsites; ++k) {
      while (k < totalsites && meth[k][i]== -1.0) {
        ++k;
      }
      if (k < totalsites) {
        if (distance_between_sites(prev_obs, sites[k]) > desert_size)  {
          for (size_t w = j+1; w < k; ++w) {
            is_desert[w] = true;
          }
        }
        j = k;
        prev_obs = sites[k];
      } else {
        for (size_t w = j+1; w < totalsites; ++w)
          is_desert[w] = true;
      }
    }
    cerr << endl;
  }
  cerr << "Scan complete" << endl;
  size_t good =0;
  for (size_t i = 0; i < totalsites; ++i) {
    if (!is_desert[i]) ++good;
  }
  cerr << good << " good sites" << endl;

  vector<Site> sites_copy;
  Site last_site;
  vector<vector<double> > meth_copy;
  const double empty_cutoff = -1.0*nleaf;

  size_t last_block_size = 0;
  for (size_t i = 0; i < totalsites; ++i) {
     if (accumulate(meth[i].begin(), meth[i].end(), 0.0) > empty_cutoff &&
        !is_desert[i]) {
      // Add new site
      sites_copy.push_back(sites[i]);
      meth_copy.push_back(meth[i]);
      ++ last_block_size;
      // Maintain blocks
      if (sites_copy.size()==1) {
        reset_points.push_back(0);
       } else if (distance_between_sites(last_site, sites[i]) > desert_size) {
        if (last_block_size < min_site_per_block) { //remover small block
          size_t start = reset_points.back();
          sites_copy.erase(sites_copy.begin()+start, sites_copy.end() );
          meth_copy.erase(meth_copy.begin()+start, meth_copy.end());
          reset_points.pop_back();
          last_block_size =
            reset_points.empty()? 0 : meth_copy.size()-reset_points.back();
        } else { //end block, and start new block
          reset_points.push_back(sites_copy.size()-1);
          last_block_size = 1;
        }
      }
      if (sites_copy.size()>0)
        last_site = sites_copy.back();
    }
  }
  reset_points.push_back(sites_copy.size());
  meth.swap(meth_copy);
  sites.swap(sites_copy);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////// PHYLOTREEPREORDER HELPERS ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
subtree_sizes_to_leaves_preorder(const vector<size_t> &subtree_sizes,
                                 vector<size_t> &leaves_preorder) {
  for (size_t i = 0; i < subtree_sizes.size(); ++i) {
    if (subtree_sizes[i] == 1) {
      leaves_preorder.push_back(i);
    }
  }
}


bool
is_binary(const vector<size_t> &subtree_sizes) {
  return (subtree_sizes[0] ==
          1 + subtree_sizes[1] + subtree_sizes[subtree_sizes[1]+1]);
}


static void
get_node_degrees(const vector<size_t> &subtree_sizes,
                 const size_t tree_start,
                 vector<size_t> &degrees) {
  if (subtree_sizes[tree_start] == 1) {
    degrees[tree_start] = 1;
  } else {
    size_t count = 1;
    if (tree_start > 0)
      degrees[tree_start] = 1;
    while (count < subtree_sizes[tree_start]) {
      size_t next_child = tree_start + count;
      get_node_degrees(subtree_sizes, next_child, degrees);
      count += subtree_sizes[next_child];
      degrees[tree_start] += 1;
    }
    assert (count == subtree_sizes[tree_start]);
  }
}


static void
get_degrees(const vector<size_t> &subtree_sizes,
            vector<size_t> &degrees) {
  degrees = vector<size_t>(subtree_sizes.size(), 0);
  get_node_degrees(subtree_sizes, 0, degrees);
}


bool
is_semi_binary(const vector<size_t> &degrees) {
  bool semi_binary = (degrees[0] ==2 || degrees[0] == 3);
  size_t i = 1;
  while (semi_binary && i < degrees.size()) {
    semi_binary = (degrees[i] ==1 || degrees[i] ==3);
    ++i;
  }
  return semi_binary;
}

size_t
leafsize(const vector<size_t> &subtree_sizes) {
  size_t n_leaf = 0;
  for (size_t i = 0; i < subtree_sizes.size(); ++i) {
    if (subtree_sizes[i] == 1)
      ++n_leaf;
  }
  return n_leaf;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////////       Likelihood          ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
double
evaluate_emission(const vector<size_t> &subtree_sizes,
                  const string states,
                  const vector<double> &hypo_prob) {
  const double TOL = 1e-10;
  double llk = 0;
  size_t leaf_count = 0;
  for (size_t i = 0; i < subtree_sizes.size(); ++i) {
    if (subtree_sizes[i] == 1 ) {
      double p = hypo_prob[leaf_count];
      if (p >= 0) {
        if (p < TOL) p = TOL;
        if (p > 1.0-TOL) p = 1.0-TOL;
        llk += (states[i] == '0') ? log(p) : log(1.0 - p);
      }
      ++ leaf_count;
    }
  }
  return llk;
}

static void
temporal_trans_prob_mat(const double exp_interval, const double rate0,
                        vector<vector<double> > &time_transition_matrix) {
  assert(rate0 > 0 && rate0 < 1 && exp_interval >1);
  const double rate1 = 1.0 - rate0;
  const double h = 1.0/exp_interval;
  time_transition_matrix = vector<vector<double> >(2, vector<double>(2, 0.0));
  time_transition_matrix[0][0] = (rate0*h + rate1)/(rate0 + rate1);
  time_transition_matrix[0][1] = 1.0 - time_transition_matrix[0][0];
  time_transition_matrix[1][1] = (rate0 + rate1*h)/(rate0 + rate1);
  time_transition_matrix[1][0] = 1.0 - time_transition_matrix[1][1];
}

static void
combined_trans_prob_mat(const double g0, const double g1,
                        const vector<vector<double> > &time_trans_mat,
                        vector<vector<vector<double> > > &combined_trans_mat) {
  // spatial transition probabilities
  vector<vector<double> >G(2, vector<double>(2, 0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  // normalization denominators
  vector<vector<double> > prev_anc_denom(2, vector<double>(2,0.0));
  for (size_t prev = 0; prev < 2; ++ prev) {
    for (size_t anc = 0; anc < 2; ++ anc) {
      prev_anc_denom[prev][anc] = std::inner_product(G[prev].begin(), G[prev].end(),
                                                     time_trans_mat[anc].begin(), 0.0);
    }
  }

  combined_trans_mat =
    vector<vector<vector<double> > >(2, vector<vector<double> >(2,vector<double>(2,0.0)));
  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      for (size_t cur = 0; cur < 2; ++ cur) {
        combined_trans_mat[prev][anc][cur] =
          G[prev][cur]*time_trans_mat[anc][cur]/prev_anc_denom[prev][anc];
      }
    }
  }
}




void
trans_deriv(const vector<vector<double> > &Q,
            const vector<vector<double> > &G,
            const double exp_branch,
            const size_t prev,
            const size_t anc,
            const size_t cur,
            const vector<vector<double> > &prev_anc_denom,
            const vector<vector<double> > &time_trans_mat,
            double &d_lam, double &d_g0, double &d_g1, double &d_T) {

  const double denom = prev_anc_denom[prev][anc];

  d_lam = exp_branch*G[prev][cur]/denom*
    (kronecker_delta(cur, 1) - time_trans_mat[anc][cur]*
     (G[prev][1]-G[prev][0])/denom);

  if (prev == 0) {
    d_g0 = time_trans_mat[anc][cur]/denom*
      (kronecker_delta(cur, 0)*2 - 1.0 -
       G[prev][cur]*(time_trans_mat[anc][0] -
                     time_trans_mat[anc][1])/denom);
  }
  if (prev == 1) {
    d_g1 = time_trans_mat[anc][cur]/denom*
      (kronecker_delta(cur, 1)*2 - 1.0 -
       G[prev][cur]*(time_trans_mat[anc][1] -
                     time_trans_mat[anc][0])/denom);
  }
  d_T = G[prev][cur]/denom*
    (Q[anc][cur] - time_trans_mat[anc][cur]*
     (G[prev][0]*Q[anc][0]+G[prev][1]*Q[anc][1])/denom);
}


static void
combined_trans_prob_mat_deriv(const double lam,
                              const double g0, const double g1,
                              const double exp_branch,
                              const vector<vector<double> > &time_trans_mat,
                              vector<vector<vector<double> > > &combined_trans_mat,
                              vector<vector<vector<double> > > &combined_trans_mat_drate,
                              vector<vector<vector<double> > > &combined_trans_mat_dg0,
                              vector<vector<vector<double> > > &combined_trans_mat_dg1,
                              vector<vector<vector<double> > > &combined_trans_mat_dT) {
  // deriv: (lam, g0, g1, T)
  vector<vector<double> > G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][1] = g1;
  G[1][0] = 1.0 - g1;
  vector<vector<double> > Q(2,vector<double>(2,0.0));
  Q[0][0] = -1.0*lam;
  Q[0][1] = lam;
  Q[1][0] = 1.0 - lam;
  Q[1][1] = lam - 1.0;

  vector<vector<double> > prev_anc_denom(2, vector<double>(2,0.0));
  for (size_t prev = 0; prev < 2; ++ prev) {
    for (size_t anc = 0; anc < 2; ++ anc) {
      prev_anc_denom[prev][anc] = std::inner_product(G[prev].begin(), G[prev].end(),
                                                     time_trans_mat[anc].begin(), 0.0);
    }
  }

  combined_trans_mat =
    vector<vector<vector<double> > >(2, vector<vector<double> >(2,vector<double>(2,0.0)));
  combined_trans_mat_drate = combined_trans_mat;
  combined_trans_mat_dg0 = combined_trans_mat;
  combined_trans_mat_dg1 = combined_trans_mat;
  combined_trans_mat_dT = combined_trans_mat;

  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      for (size_t cur = 0; cur < 2; ++ cur) {
        combined_trans_mat[prev][anc][cur] =
          G[prev][cur]*time_trans_mat[anc][cur]/prev_anc_denom[prev][anc];
      }
    }
  }

  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      for (size_t cur = 0; cur < 2; ++ cur) {
        trans_deriv(Q, G, exp_branch, prev, anc, cur,
                    prev_anc_denom, time_trans_mat,
                    combined_trans_mat_drate[prev][anc][cur],
                    combined_trans_mat_dg0[prev][anc][cur],
                    combined_trans_mat_dg1[prev][anc][cur],
                    combined_trans_mat_dT[prev][anc][cur]);
      }
    }
  }
}


static void
collect_transition_matrices(const double g0, const double g1,
                            const double rate0, const vector<double> &branches,
                            vector<vector<vector<double> > > &time_trans_mats,
                            vector<vector<vector<vector<double> > > > &combined_trans_mats) {
  assert(g0 > 0 && g0 <1 && g1 > 0 && g1 <1);
  vector<double> exp_branches;
  for (size_t i = 0; i < branches.size(); ++i) {
    exp_branches.push_back(exp(branches[i]));
  }
  const size_t n_nodes = branches.size();
  time_trans_mats = vector<vector<vector<double> > >(n_nodes, vector<vector<double> >(2,vector<double>(2,0.0)));
  combined_trans_mats = vector<vector<vector<vector<double> > > >(n_nodes, vector<vector<vector<double> > >(2, vector<vector<double> >(2, vector<double>(2,0.0))));

  for (size_t i = 1; i < branches.size(); ++i) {
    temporal_trans_prob_mat(exp_branches[i], rate0, time_trans_mats[i]);
    combined_trans_prob_mat(g0, g1, time_trans_mats[i],combined_trans_mats[i]);
  }
}


static void
collect_transition_matrices_deriv(const double g0, const double g1,
                                  const double rate0, const vector<double> &branches,
                                  vector<vector<vector<double> > > &time_trans_mats,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_drate,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_dg0,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_dg1,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_dT) {
  assert(g0 > 0 && g0 <1 && g1 > 0 && g1 <1);
  const size_t n_nodes = branches.size();
  vector<double> exp_branches;
  for (size_t i = 0; i < branches.size(); ++i) {
    exp_branches.push_back(exp(branches[i]));
  }

  vector<vector<vector<double> > > unit(2,vector<vector<double> >(2, vector<double>(2, 0.0)));
  time_trans_mats = vector<vector<vector<double> > > (n_nodes, vector<vector<double> >(2, vector<double>(2, 0.0)));
  combined_trans_mats = vector<vector<vector<vector<double> > > >(n_nodes, unit);
  combined_trans_mats_drate = vector<vector<vector<vector<double> > > >(n_nodes, unit);
  combined_trans_mats_dg0 = vector<vector<vector<vector<double> > > >(n_nodes, unit);
  combined_trans_mats_dg1 = vector<vector<vector<vector<double> > > >(n_nodes, unit);
  combined_trans_mats_dT = vector<vector<vector<vector<double> > > >(n_nodes, unit);

  for (size_t i = 1; i < branches.size(); ++i) {
    temporal_trans_prob_mat(exp_branches[i], rate0, time_trans_mats[i]);
    combined_trans_prob_mat_deriv(rate0, g0, g1, exp_branches[i], time_trans_mats[i],
                                  combined_trans_mats[i], combined_trans_mats_drate[i],
                                  combined_trans_mats_dg0[i], combined_trans_mats_dg1[i],
                                  combined_trans_mats_dT[i]);
  }
}


// hypo_prob: at leaf nodes only
double loglik_tree_start(const double pi0,
                         const vector<size_t> &subtree_sizes,
                         const vector<vector<vector<double> > > &time_trans_mats,
                         const string &states,
                         const vector<double> &hypo_prob) {
  const size_t n_nodes = subtree_sizes.size();

  double llk = (states[0] == '0')? log(pi0) : log(1.0-pi0);
  for (size_t i = 0; i < n_nodes; ++i) {
    size_t anc = (states[i] == '0')? 0 : 1;
    if (subtree_sizes[i] > 1) {
      size_t count = 1;
      while (count < subtree_sizes[i]) {
        size_t child = i + count;
        size_t cur = (states[child] == '0')? 0 : 1;
        llk += log(time_trans_mats[child][anc][cur]);
        count += subtree_sizes[child];
      }
    }
  }
  llk += evaluate_emission(subtree_sizes, states, hypo_prob);
  return llk;
}

// hypo_prob: at leaf nodes only
double
loglik_tree_start_deriv(const double pi0, const double rate0,
                        const vector<size_t> &subtree_sizes,
                        const vector<double> &exp_branches,
                        const vector<vector<vector<double> > > &time_trans_mats,
                        const string &states,
                        const vector<double> &hypo_prob,
                        double &d_pi, double &d_rate,
                        vector<double> &d_T) {

  const size_t n_nodes = subtree_sizes.size();
  double llk = (states[0] == '0')? log(pi0) : log(1.0-pi0);
  vector<vector<double> > Q (2, vector<double>(2,0.0));
  Q[0][0] = -1.0*rate0;
  Q[0][1] = rate0;
  Q[1][0] = 1.0 - rate0;
  Q[1][1] = rate0 - 1.0;


  d_rate = 0.0;
  d_T = vector<double>(n_nodes, 0.0);

  for (size_t i = 0; i < n_nodes; ++i) {
    size_t anc = (states[i] == '0') ? 0 : 1;
    if (subtree_sizes[i] > 1) {
      size_t count = 1;
      while (count < subtree_sizes[i]) {
        size_t child = i + count;
        size_t cur = (states[child] == '0') ? 0 : 1;
        const double transprob = time_trans_mats[child][anc][cur];
        llk += log(transprob);
        d_rate += exp_branches[child]*(kronecker_delta(cur, 1)*2 - 1.0)/transprob;
        d_T[child] = Q[anc][cur]/transprob;
        count += subtree_sizes[child];
      }
    }
  }
  llk += evaluate_emission(subtree_sizes, states, hypo_prob);

  d_rate = d_rate*exp(llk);
  d_pi = (states[0] == '0')? exp(llk)/pi0 : -1.0*exp(llk)/(1.0-pi0);
  for (size_t i = 1; i < n_nodes; ++i) {
    d_T[i] = d_T[i]*exp(llk);
  }

  return llk;
}

double loglik_tree_transition(const vector<vector<double> > &G,
                              const vector<size_t> &subtree_sizes,
                              const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                              const string prev_states,
                              const string cur_states) {
  const size_t n_nodes = subtree_sizes.size();

  size_t prev = (prev_states[0] == '0') ? 0 : 1;
  size_t cur = (cur_states[0] == '0') ? 0 : 1;
  double llk = log(G[prev][cur]);

  for (size_t i = 0; i < n_nodes; ++i) {
    size_t anc = cur_states[i];
    if (subtree_sizes[i] > 1) {
      size_t count = 1;
      while (count < subtree_sizes[i]) {
        size_t child = i + count;
        prev = (prev_states[child] == '0') ? 0 : 1;
        cur = (cur_states[child] == '0') ? 0 : 1;
        llk += log(combined_trans_mats[child][prev][anc][cur]);
        count += subtree_sizes[child];
      }
    }
  }
  return llk;
}


double
loglik_tree_transition_deriv(vector<vector<double> > &G,
                             const vector<size_t> &subtree_sizes,
                             const vector<vector<vector<double> > > &time_trans_mats,
                             const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                             const vector<vector<vector<vector<double> > > > &combined_trans_mats_drate,
                             const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg0,
                             const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg1,
                             const vector<vector<vector<vector<double> > > > &combined_trans_mats_dT,
                             const string prev_states,
                             const string cur_states,
                             double &trans_deriv_rate,
                             double &trans_deriv_g0,
                             double &trans_deriv_g1,
                             vector<double> &trans_deriv_T) {

  const size_t n_nodes = subtree_sizes.size();
  trans_deriv_T = vector<double>(n_nodes, 0.0);
  trans_deriv_rate = 0;
  trans_deriv_g0 = 0;
  trans_deriv_g1 = 0;

  size_t prev = (prev_states[0] == '0') ? 0 : 1;
  size_t cur = (cur_states[0] == '0') ? 0 : 1;
  double llk = log(G[prev][cur]);
  trans_deriv_g0 = kronecker_delta(prev, 0)*(kronecker_delta(cur, 0)*2 - 1.0)/G[prev][cur];
  trans_deriv_g1 = kronecker_delta(prev, 1)*(kronecker_delta(cur, 1)*2 - 1.0)/G[prev][cur];

  for (size_t i = 0; i < n_nodes; ++i) {
    size_t anc = (cur_states[i] == '0') ? 0 : 1;
    if (subtree_sizes[i] > 1) {
      size_t count = 1;
      while (count < subtree_sizes[i]) {
        size_t child = i + count;
        prev = (prev_states[child] == '0') ? 0 : 1;
        cur = (cur_states[child] == '0') ? 0 : 1;
        const double transprob = combined_trans_mats[child][prev][anc][cur];
        llk += log(transprob);
        trans_deriv_T[child] = combined_trans_mats_dT[child][prev][anc][cur]/transprob;
        trans_deriv_rate += combined_trans_mats_drate[child][prev][anc][cur]/transprob;
        if (prev == 0) {
          trans_deriv_g0 += combined_trans_mats_dg0[child][prev][anc][cur]/transprob;
        } else {
          trans_deriv_g1 += combined_trans_mats_dg1[child][prev][anc][cur]/transprob;
        }
        count += subtree_sizes[child];
      }
    }
  }

  trans_deriv_rate = trans_deriv_rate*exp(llk);
  trans_deriv_g0 = trans_deriv_g0*exp(llk);
  trans_deriv_g1 = trans_deriv_g1*exp(llk);
  for (size_t i = 1; i < n_nodes; ++i) {
    trans_deriv_T[i] = trans_deriv_T[i]*exp(llk);
  }

  return llk;
}


string
dec2binary(const size_t len, size_t n) {
  assert (n < pow(2, len));
  string result;
  do result.push_back( '0' + (n & 1) );
  while (n >>= 1);
  while (result.length()<len) {
    result.push_back('0');
  }
  reverse(result.begin(), result.end());
  return result;
}


size_t
binary2dec(const string states) {
  size_t id = 0;
  for (size_t i = 0; i < states.length(); ++i) {
    if (states[states.length()-i-1] == '1')
      id += pow(2, i);
  }
  return id;
}

static void
hme_trans_prob_mat(const double pi0, const double g0, const double g1,
                   const double rate0,
                   const vector<double> &branches,
                   const vector<size_t> &subtree_sizes,
                   const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                   vector<vector<double> > &hme_log_tpm) {
  size_t n_nodes = subtree_sizes.size();
  size_t n_hme = pow(2, n_nodes);
  hme_log_tpm = vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  vector<string> all_states;
  for (size_t i = 0; i < n_hme; ++i) {
    all_states.push_back(dec2binary(n_nodes, i));
  }

  vector<vector<double> >G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  for (size_t i = 0; i < n_hme; ++i) {
    for (size_t j = 0; j < n_hme; ++j) {
      hme_log_tpm[i][j] =
        loglik_tree_transition(G, subtree_sizes, combined_trans_mats,
                               all_states[i], all_states[j]);
    }
  }
}


static void
hme_trans_prob_mat_deriv(const double g0, const double g1,
                         const double rate0,
                         const vector<double> &branches,
                         const vector<size_t> &subtree_sizes,
                         const vector<vector<vector<double> > > &time_trans_mats,
                         const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                         const vector<vector<vector<vector<double> > > > &combined_trans_mats_drate,
                         const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg0,
                         const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg1,
                         const vector<vector<vector<vector<double> > > > &combined_trans_mats_dT,
                         vector<vector<double> > &hme_log_tpm,
                         vector<vector<double> > &hme_trans_drate,
                         vector<vector<double> > &hme_trans_dg0,
                         vector<vector<double> > &hme_trans_dg1,
                         vector<vector<vector<double> > > &hme_trans_dT) {
  size_t n_nodes = subtree_sizes.size();
  size_t n_hme = pow(2, n_nodes);
  hme_log_tpm = vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  hme_trans_dg0 = vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  hme_trans_dg1 = vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  hme_trans_drate = vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  hme_trans_dT = vector<vector<vector<double> > >(n_hme, vector<vector<double> >(n_hme, vector<double>(n_nodes, 0.0)));

  vector<string> all_states;
  for (size_t i = 0; i < n_hme; ++i) {
    all_states.push_back(dec2binary(n_nodes, i));
  }

  vector<vector<double> >G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  for (size_t i = 0; i < n_hme; ++i) {
    for (size_t j = 0; j < n_hme; ++j) {
      hme_log_tpm[i][j] =
        loglik_tree_transition_deriv(G, subtree_sizes,time_trans_mats,
                                     combined_trans_mats, combined_trans_mats_drate,
                                     combined_trans_mats_dg0, combined_trans_mats_dg1,
                                     combined_trans_mats_dT,
                                     all_states[i], all_states[j], hme_trans_drate[i][j],
                                     hme_trans_dg0[i][j], hme_trans_dg1[i][j],
                                     hme_trans_dT[i][j]);
    }
  }
}


// restricted transition states
// [0]: self; [1,...,n_nodes] change at each node;
// [n_node+1,2*n_nodes-1] change at root and each other node
static void
get_rstr_state_ids(const string states,
                   vector<size_t> &rstr_state_ids) {
  const size_t id = binary2dec(states);
  const size_t n_nodes = states.length();
  rstr_state_ids.push_back(id);
  for (size_t i = 0; i < n_nodes; ++i) {
    size_t delta =  pow(2, n_nodes-i-1);
    size_t sign = (states[i]=='0') ? 1 : -1;
    rstr_state_ids.push_back(id + sign*delta);
  }
  size_t root_delta =  pow(2, n_nodes-1);
  size_t root_sign = (states[0]=='0') ? 1 : -1;
  for (size_t i =1; i < n_nodes; ++i) {
    size_t delta =  pow(2, n_nodes-i-1);
    size_t sign = (states[i]=='0') ? 1 : -1;
    rstr_state_ids.push_back(id + sign*delta + root_delta*root_sign);
  }
}


// HME transitions correspond to at most one change in one branch,
// and zero or one change at the root
static void
rstr_hme_trans_prob_mat_deriv(const double pi0, const double g0,
                              const double g1, const double rate0,
                              const vector<double> &branches,
                              const vector<size_t> &subtree_sizes,
                              const vector<vector<vector<double> > > &time_trans_mats,
                              const vector<string> &all_states,
                              const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                              const vector<vector<vector<vector<double> > > > &combined_trans_mats_drate,
                              const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg0,
                              const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg1,
                              const vector<vector<vector<vector<double> > > > &combined_trans_mats_dT,
                              vector<vector<size_t> > &neighbors,
                              vector<vector<size_t> > &place_in_neighbor,
                              vector<vector<double> > &rstr_hme_log_tpm,
                              vector<vector<double> > &rstr_hme_tpm_drate,
                              vector<vector<double> > &rstr_hme_tpm_dg0,
                              vector<vector<double> > &rstr_hme_tpm_dg1,
                              vector<vector<vector<double> > > &rstr_hme_tpm_dT) {
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_hme = pow(2, n_nodes);
  const size_t n_neighbor = 2*n_nodes;

  rstr_hme_log_tpm = vector<vector<double> >(n_hme, vector<double>(n_neighbor, 0.0));
  rstr_hme_tpm_drate = vector<vector<double> >(n_hme, vector<double>(n_neighbor, 0.0));
  rstr_hme_tpm_dg0 = vector<vector<double> >(n_hme, vector<double>(n_neighbor, 0.0));
  rstr_hme_tpm_dg1 = vector<vector<double> >(n_hme, vector<double>(n_neighbor, 0.0));
  rstr_hme_tpm_dT = vector<vector<vector<double> > >(n_hme, vector<vector<double> >(n_neighbor, vector<double>(n_nodes, 0.0)));


  vector<vector<double> >G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  for (size_t i = 0; i < n_hme; ++i) {
    vector<size_t> rstr_state_ids;
    get_rstr_state_ids(all_states[i], rstr_state_ids);
    assert (rstr_state_ids.size()== n_neighbor);
    neighbors.push_back(rstr_state_ids);
    double row_sum = 0;

    vector<double> trans_deriv_rate(n_neighbor, 0.0);
    vector<double> trans_deriv_g0(n_neighbor, 0.0);
    vector<double> trans_deriv_g1(n_neighbor, 0.0);
    vector<vector<double> > trans_deriv_T(n_neighbor, vector<double>(n_nodes,0.0));

    // Compute transition probabilities
    for (size_t j = 0; j < n_neighbor; ++j) {
      const double log_prob =
        loglik_tree_transition_deriv(G, subtree_sizes, time_trans_mats,
                                     combined_trans_mats, combined_trans_mats_drate,
                                     combined_trans_mats_dg0, combined_trans_mats_dg1,
                                     combined_trans_mats_dT,
                                     all_states[i], all_states[rstr_state_ids[j]],
                                     trans_deriv_rate[j], trans_deriv_g0[j],
                                     trans_deriv_g1[j], trans_deriv_T[j]);
      rstr_hme_log_tpm[i][j] = log_prob;
      row_sum = log_sum_log(row_sum, log_prob);
    }

    //Compute restricted transition prob derivatives
    double sum_deriv_rate = 0.0;
    double sum_deriv_g0 = 0.0;
    double sum_deriv_g1 = 0.0;
    vector<double> sum_deriv_T(n_nodes, 0.0);
    for (size_t j = 0; j < n_neighbor; ++j) {
      sum_deriv_rate += trans_deriv_rate[j];
      sum_deriv_g0 += trans_deriv_g0[j];
      sum_deriv_g1 += trans_deriv_g1[j];
      for (size_t k = 1; k < n_nodes; ++k) {
        sum_deriv_T[k] += trans_deriv_T[j][k];
      }
    }
    for (size_t j = 0; j < n_neighbor; ++j ) {
      rstr_hme_tpm_drate[i][j] = trans_deriv_rate[j]/exp(row_sum) -
        exp(rstr_hme_log_tpm[i][j]-row_sum*2)*sum_deriv_rate;

      rstr_hme_tpm_dg0[i][j] = trans_deriv_g0[j]/exp(row_sum) -
        exp(rstr_hme_log_tpm[i][j]-row_sum*2)*sum_deriv_g0;

      rstr_hme_tpm_dg1[i][j] = trans_deriv_g1[j]/exp(row_sum) -
        exp(rstr_hme_log_tpm[i][j]-row_sum*2)*sum_deriv_g1;

      for (size_t k = 1; k < n_nodes; ++k) {
        rstr_hme_tpm_dT[i][j][k] = trans_deriv_T[j][k]/exp(row_sum) -
          exp(rstr_hme_log_tpm[i][j]-row_sum*2)*sum_deriv_T[k];
      }
    }

    // Normalize transition probabiliteis
    for (size_t j = 0; j < rstr_state_ids.size(); ++j) {
      rstr_hme_log_tpm[i][j] = rstr_hme_log_tpm[i][j]-row_sum;
    }
  }

  // my slots in my neighbor's lists
  place_in_neighbor =
    vector<vector<size_t> >(n_hme, vector<size_t>(2*n_nodes, 0.0));
  for (size_t i = 0; i < n_hme; ++i) {
    for (size_t j = 0; j < neighbors[0].size(); ++j) {
      const size_t neighbor = neighbors[i][j];
      const size_t myid = distance(neighbors[neighbor].begin(),
                                   find(neighbors[neighbor].begin(),
                                        neighbors[neighbor].end(), i));
      place_in_neighbor[i][j] = myid;
    }
  }
}


static void
rstr_hme_trans_prob_mat(const double pi0, const double g0,
                        const double g1, const double rate0,
                        const vector<double> &branches,
                        const vector<size_t> &subtree_sizes,
                        const vector<string> &all_states,
                        const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                        vector<vector<size_t> > &neighbors,
                        vector<vector<size_t> > &place_in_neighbor,
                        vector<vector<double> > &rstr_hme_log_tpm) {
  size_t n_nodes = subtree_sizes.size();
  size_t n_hme = pow(2, n_nodes);

  rstr_hme_log_tpm =
    vector<vector<double> >(n_hme, vector<double>(2*n_nodes, 0.0));

  vector<vector<double> >G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  for (size_t i = 0; i < n_hme; ++i) {
    vector<size_t> rstr_state_ids;
    get_rstr_state_ids(all_states[i], rstr_state_ids);
    assert (rstr_state_ids.size()== 2*n_nodes);
    neighbors.push_back(rstr_state_ids);
    double row_sum = 0;
    // Compute transition probabilities
    for (size_t j = 0; j < rstr_state_ids.size(); ++j) {
      const double log_prob =
        loglik_tree_transition(G, subtree_sizes, combined_trans_mats,
                               all_states[i], all_states[rstr_state_ids[j]]);
      rstr_hme_log_tpm[i][j] = log_prob;
      row_sum = log_sum_log(row_sum, log_prob);
    }
    // Normalize
    for (size_t j = 0; j < rstr_state_ids.size(); ++j) {
      rstr_hme_log_tpm[i][j] = rstr_hme_log_tpm[i][j]-row_sum;
      assert (rstr_hme_log_tpm[i][j] <0 );
    }
  }

  // my slots in my neighbor's lists
  // j.neighbor(i.myid[j]) = i
  place_in_neighbor =
    vector<vector<size_t> >(n_hme, vector<size_t>(2*n_nodes, 0.0));
  for (size_t i = 0; i < n_hme; ++i) {
    for (size_t j = 0; j < neighbors[0].size(); ++j) {
      const size_t neighbor = neighbors[i][j];
      const size_t myid = distance(neighbors[neighbor].begin(),
                                   find(neighbors[neighbor].begin(),
                                        neighbors[neighbor].end(), i));
      place_in_neighbor[i][j] = myid;
    }
  }
}






static void
forward_variables(const double pi0, const vector<size_t> &subtree_sizes,
                  const vector<string> &all_states,
                  const vector<vector<vector<double> > > &time_trans_mats,
                  const vector<vector<double> > &hme_log_tpm,
                  const vector<vector<double> > &meth_prob_table,
                  const size_t pos,
                  vector<double> &forward) {
  const size_t n_hme = all_states.size();
  forward = vector<double> (n_hme, 0.0);
  if (pos == 0) {
    for (size_t i = 0; i < n_hme; ++i) {
      forward[i] = loglik_tree_start(pi0, subtree_sizes, time_trans_mats,
                                     all_states[i], meth_prob_table[0]);
    }
  } else {
    vector<double> ff;
    forward_variables(pi0, subtree_sizes, all_states, time_trans_mats,
                      hme_log_tpm, meth_prob_table, pos-1, ff);
    for (size_t j = 0; j < n_hme; ++ j) {
      for (size_t i = 0; i < n_hme; ++ i) {
        forward[j] = log_sum_log(forward[j], ff[i] + hme_log_tpm[i][j]);
      }
      forward[j] += evaluate_emission(subtree_sizes, all_states[j],
                                      meth_prob_table[pos]);
    }
  }
  if (pos %1000 == 0) cerr << ".";
}



static void
forward_variables_deriv(const double pi0, const double rate0,
                        const vector<size_t> &subtree_sizes,
                        const vector<double> &exp_branches,
                        const vector<string> &all_states,
                        const vector<vector<vector<double> > > &time_trans_mats,
                        const vector<vector<double> > &hme_log_tpm,
                        const vector<vector<double> > &hme_trans_drate,
                        const vector<vector<double> > &hme_trans_dg0,
                        const vector<vector<double> > &hme_trans_dg1,
                        const vector<vector<vector<double> > > &hme_trans_dT,
                        const vector<vector<double> > &meth_prob_table,
                        const size_t pos,
                        vector<double> &forward,
                        vector<vector<Logscale> > &log_forward_deriv) {
  const size_t n_hme = all_states.size();
  const size_t n_nodes = subtree_sizes.size();
  forward = vector<double> (n_hme, 0.0);
  log_forward_deriv = vector<vector<Logscale> >(n_hme, vector<Logscale>(n_nodes+4, Logscale(0.0,0.0)));

  if (pos == 0) {
    for (size_t i = 0; i < n_hme; ++i) {
      double d_pi, d_rate;
      vector<double> d_T;
      forward[i] = loglik_tree_start_deriv(pi0, rate0, subtree_sizes, exp_branches,
                                           time_trans_mats, all_states[i],
                                           meth_prob_table[0], d_pi, d_rate, d_T);
      log_forward_deriv[i][0].logval = log(abs(d_pi));
      log_forward_deriv[i][0].symbol = sign(d_pi);
      log_forward_deriv[i][1].logval = log(abs(d_rate));
      log_forward_deriv[i][1].symbol = sign(d_rate);
      for (size_t j = 1; j < n_nodes; ++j) {
        log_forward_deriv[i][4+j].logval = log(abs(d_T[j]));
        log_forward_deriv[i][4+j].symbol = sign(d_T[j]);
      }
    }
  } else {
    vector<double> ff;
    vector<vector<Logscale> > ffdd;
    forward_variables_deriv(pi0, rate0, subtree_sizes, exp_branches,
                            all_states, time_trans_mats, hme_log_tpm, hme_trans_drate,
                            hme_trans_dg0, hme_trans_dg1,  hme_trans_dT,
                            meth_prob_table, pos-1, ff, ffdd);

    for (size_t j = 0; j < n_hme; ++ j) { //current site
      for (size_t i = 0; i < n_hme; ++ i) { //previous site
        double log_trans_prob = hme_log_tpm[i][j];
        forward[j] = log_sum_log(forward[j], ff[i] + log_trans_prob);

        vector<double> trans_derivs(n_nodes+4, 0.0);
        trans_derivs[1] = hme_trans_drate[i][j];
        trans_derivs[2] = hme_trans_dg0[i][j];
        trans_derivs[3] = hme_trans_dg1[i][j];
        for (size_t k = 1; k < n_nodes; ++k) {
          trans_derivs[4+k] = hme_trans_dT[i][j][k];
        }

        Logscale first, second, inc, result;
        for (size_t k = 0; k < 4+n_nodes; ++k) {
        first.symbol = ffdd[i][k].symbol;
        if (first.symbol != 0.0)
          first.logval = ffdd[i][k].logval + log_trans_prob;
        second.symbol = sign(trans_derivs[k]);
        if (second.symbol != 0.0)
          second.logval = ff[i] + log(abs(trans_derivs[k]));

        log_sum_log_sign(first, second, inc);
        log_sum_log_sign(log_forward_deriv[j][k], inc, result);
        log_forward_deriv[j][k].logval = result.logval;
        log_forward_deriv[j][k].symbol = result.symbol;
        }
      }
      double emi_log_prob = evaluate_emission(subtree_sizes, all_states[j], meth_prob_table[pos]);
      forward[j] += emi_log_prob;
      for (size_t k = 0; k < n_nodes+4; ++k) {
        if (k!=4)
          log_forward_deriv[j][k].logval += emi_log_prob;
      }
    }
  }
}


double
llk_by_forward_algorithm(const double pi0,
                         const vector<size_t> &subtree_sizes,
                         const vector<string> &all_states,
                         const vector<vector<vector<double> > > &time_trans_mats,
                         const vector<vector<double> > &hme_log_tpm,
                         const vector<vector<double> > &meth_prob_table) {
  const size_t n_hme = all_states.size();
  const size_t end = meth_prob_table.size()-1;
  double llk = 0;
  vector<double> forward;
  forward_variables(pi0, subtree_sizes, all_states, time_trans_mats,
                    hme_log_tpm, meth_prob_table, end, forward);
  for (size_t i = 0; i < n_hme; ++i ){
    llk = log_sum_log(llk, forward[i]);
  }
  return llk;
}


double
llk_by_forward_algorithm_deriv(const double pi0, const double rate0,
                               const vector<size_t> &subtree_sizes,
                               const vector<double> &exp_branches,
                               const vector<string> &all_states,
                               const vector<vector<vector<double> > > &time_trans_mats,
                               const vector<vector<double> > &hme_log_tpm,
                               const vector<vector<double> > &hme_trans_drate,
                               const vector<vector<double> > &hme_trans_dg0,
                               const vector<vector<double> > &hme_trans_dg1,
                               const vector<vector<vector<double> > > &hme_trans_dT,
                               const vector<vector<double> > &meth_prob_table,
                               vector<double> &deriv) {

  const size_t n_hme = all_states.size();
  const size_t n_nodes = subtree_sizes.size();
  const size_t end = meth_prob_table.size()-1;
  double llk = 0;
  vector<double> forward;
  vector<vector<Logscale> > forward_deriv;

  forward_variables_deriv(pi0, rate0, subtree_sizes, exp_branches, all_states,
                          time_trans_mats, hme_log_tpm, hme_trans_drate,
                          hme_trans_dg0, hme_trans_dg1, hme_trans_dT,
                          meth_prob_table, end, forward, forward_deriv);

  deriv = vector<double>(n_nodes+4, 0.0);
  vector<Logscale>deriv_copy(n_nodes+4, Logscale(0.0, 0.0));
  Logscale val;
  for (size_t i = 0; i < n_hme; ++i ){
    llk = log_sum_log(llk, forward[i]);
    for (size_t j = 0; j < deriv.size(); ++j){
      log_sum_log_sign(deriv_copy[j], forward_deriv[i][j],val);
      deriv_copy[j].logval = val.logval;
      deriv_copy[j].symbol = val.symbol;
    }
  }
  for (size_t i = 0; i < deriv.size(); ++i) {
    if (i!=4)
      deriv[i] = exp(deriv_copy[i].logval-llk)*deriv_copy[i].symbol;
  }
  return llk;
}


static void
forward_variables_rstr(const double pi0, const vector<size_t> &subtree_sizes,
                       const vector<string> &all_states,
                       const vector<vector<vector<double> > > &time_trans_mats,
                       const vector<vector<size_t> > &neighbors,
                       const vector<vector<size_t> > &place_in_neighbor,
                       const vector<vector<double> > &rstr_hme_log_tpm,
                       const vector<vector<double> > &meth_prob_table,
                       const size_t pos,
                       vector<double> &forward) {
  forward.clear();
  const size_t n_hme = all_states.size();

  if (pos == 0) {
    forward = vector<double>(n_hme, 0.0);
    for (size_t i = 0; i < n_hme; ++i) {
      forward[i] = loglik_tree_start(pi0, subtree_sizes, time_trans_mats,
                                     all_states[i], meth_prob_table[0]);
    }
  } else {
    vector<double> ff;
    forward_variables_rstr(pi0, subtree_sizes, all_states,
                           time_trans_mats, neighbors,
                           place_in_neighbor, rstr_hme_log_tpm,
                           meth_prob_table, pos-1, ff);
    forward = vector<double>(n_hme, 0.0);
    for (size_t j = 0; j < n_hme; ++ j) { //current site
      for (size_t i = 0; i < neighbors[0].size(); ++ i) { //previous site
        const size_t prev_state_id = neighbors[j][i];
        const size_t idx_in_neighbor = place_in_neighbor[j][i];
        double ele = ff[prev_state_id] +
          rstr_hme_log_tpm[prev_state_id][idx_in_neighbor] ;
        forward[j] = log_sum_log(forward[j], ele);
      }
      forward[j] += evaluate_emission(subtree_sizes, all_states[j],
                                      meth_prob_table[pos]);
    }
  }
  if ((pos+1)%1000 == 0) cerr << ".";
}



static void
forward_variables_rstr_deriv(const double pi0, const double rate0,
                             const vector<size_t> &subtree_sizes,
                             const vector<double> &exp_branches,
                             const vector<string> &all_states,
                             const vector<vector<vector<double> > > &time_trans_mats,
                             const vector<vector<size_t> > &neighbors,
                             const vector<vector<size_t> > &place_in_neighbor,
                             const vector<vector<double> > &rstr_hme_log_tpm,
                             const vector<vector<double> > &rstr_hme_trans_drate,
                             const vector<vector<double> > &rstr_hme_trans_dg0,
                             const vector<vector<double> > &rstr_hme_trans_dg1,
                             const vector<vector<vector<double> > > &rstr_hme_trans_dT,
                             const vector<vector<double> > &meth_prob_table,
                             const size_t pos,
                             vector<double> &forward,
                             vector<vector<Logscale> > &log_forward_deriv) {
  const size_t n_hme = all_states.size();
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_neighbor = neighbors[0].size();
  forward = vector<double> (n_hme, 0.0);
  log_forward_deriv = vector<vector<Logscale> >(n_hme, vector<Logscale>(n_nodes+4, Logscale(0.0, 0.0)));

  // begin
  if (pos == 0) {
    for (size_t i = 0; i < n_hme; ++i) {
      double d_pi, d_rate;
      vector<double> d_T;
      forward[i] = loglik_tree_start_deriv(pi0, rate0, subtree_sizes, exp_branches,
                                           time_trans_mats, all_states[i],
                                           meth_prob_table[0], d_pi, d_rate, d_T);
      log_forward_deriv[i][0].logval = log(abs(d_pi));
      log_forward_deriv[i][0].symbol = sign(d_pi);
      log_forward_deriv[i][1].logval = log(abs(d_rate));
      log_forward_deriv[i][1].symbol = sign(d_rate);

      for (size_t j = 1; j < n_nodes; ++j) {
        log_forward_deriv[i][4+j].logval = log(abs(d_T[j]));
        log_forward_deriv[i][4+j].symbol = sign(d_T[j]);
      }
    }
  } else {
    vector<double> ff;
    vector<vector<Logscale> > ffdd;
    forward_variables_rstr_deriv(pi0, rate0, subtree_sizes, exp_branches,
                                 all_states, time_trans_mats, neighbors,
                                 place_in_neighbor, rstr_hme_log_tpm,
                                 rstr_hme_trans_drate, rstr_hme_trans_dg0,
                                 rstr_hme_trans_dg1, rstr_hme_trans_dT,
                                 meth_prob_table, pos-1, ff, ffdd);

    for (size_t j = 0; j < n_hme; ++ j) { //current site
      for (size_t i = 0; i < n_neighbor; ++ i) { //previous site
        const size_t prev_state_id = neighbors[j][i];
        const size_t idx_in_neighbor = place_in_neighbor[j][i];
        const double log_trans_prob = rstr_hme_log_tpm[prev_state_id][idx_in_neighbor];
        forward[j] = log_sum_log(forward[j], ff[prev_state_id] +  log_trans_prob);

        vector<double> trans_derivs(n_nodes+4, 0.0);
        trans_derivs[1] = rstr_hme_trans_drate[prev_state_id][idx_in_neighbor];
        trans_derivs[2] = rstr_hme_trans_dg0[prev_state_id][idx_in_neighbor];
        trans_derivs[3] = rstr_hme_trans_dg1[prev_state_id][idx_in_neighbor];
        for (size_t k = 1; k < n_nodes; ++k) {
          trans_derivs[4+k] = rstr_hme_trans_dT[prev_state_id][idx_in_neighbor][k];
        }

        //update derivatives here
        Logscale first, second, inc, result;
        for (size_t k = 0; k < 4+n_nodes; ++k) {
          first.symbol = ffdd[prev_state_id][k].symbol;
          if (first.symbol != 0.0)
            first.logval = ffdd[prev_state_id][k].logval + log_trans_prob;
          second.symbol = sign(trans_derivs[k]);
          if (second.symbol != 0.0)
            second.logval = ff[prev_state_id] + log(abs(trans_derivs[k]));

          log_sum_log_sign(first, second, inc);
          log_sum_log_sign(log_forward_deriv[j][k], inc, result);
          log_forward_deriv[j][k].logval = result.logval;
          log_forward_deriv[j][k].symbol = result.symbol;
        }
      }
      double log_emission = evaluate_emission(subtree_sizes, all_states[j],
                                              meth_prob_table[pos]);
      forward[j] += log_emission;
      for (size_t n = 0; n < n_nodes+4; ++n) {
        if (n!=4) {
          log_forward_deriv[j][n].logval += log_emission;
        }
      }
    }
  }
  if ((pos+1)%100 == 0) cerr << ".";
}




double
llk_by_rstr_forward_algorithm(const double pi0,
                              const vector<size_t> &subtree_sizes,
                              const vector<string> &all_states,
                              const vector<vector<vector<double> > > &time_trans_mats,
                              const vector<vector<size_t> > &neighbors,
                              const vector<vector<size_t> > &place_in_neighbor,
                              const vector<vector<double> > &rstr_hme_log_tpm,
                              const vector<vector<double> > &meth_prob_table) {
  const size_t n_hme = all_states.size();
  const size_t end = meth_prob_table.size()-1;
  double llk = 0;
  vector<double> forward;
  forward_variables_rstr(pi0, subtree_sizes, all_states, time_trans_mats,
                         neighbors, place_in_neighbor, rstr_hme_log_tpm,
                         meth_prob_table, end, forward);
  for (size_t i = 0; i < n_hme; ++i ){
    llk = log_sum_log(llk, forward[i]);
  }
  return llk;
}



double
llk_by_rstr_forward_algorithm_deriv(const double pi0, const double rate0,
                                    const vector<size_t> &subtree_sizes,
                                    const vector<double> &exp_branches,
                                    const vector<string> &all_states,
                                    const vector<vector<vector<double> > > &time_trans_mats,
                                    const vector<vector<size_t> > &neighbors,
                                    const vector<vector<size_t> > &place_in_neighbor,
                                    const vector<vector<double> > &rstr_hme_log_tpm,
                                    const vector<vector<double> > &rstr_hme_trans_drate,
                                    const vector<vector<double> > &rstr_hme_trans_dg0,
                                    const vector<vector<double> > &rstr_hme_trans_dg1,
                                    const vector<vector<vector<double> > > &rstr_hme_trans_dT,
                                    const vector<vector<double> > &meth_prob_table,
                                    vector<double> &deriv) {

  const size_t n_hme = all_states.size();
  const size_t n_nodes = subtree_sizes.size();
  const size_t end = meth_prob_table.size()-1;
  double llk = 0;
  vector<double> forward;
  vector<vector<Logscale> > forward_deriv;

  forward_variables_rstr_deriv(pi0, rate0, subtree_sizes, exp_branches, all_states,
                               time_trans_mats, neighbors, place_in_neighbor,
                               rstr_hme_log_tpm, rstr_hme_trans_drate,
                               rstr_hme_trans_dg0, rstr_hme_trans_dg1,
                               rstr_hme_trans_dT, meth_prob_table, end,
                               forward, forward_deriv);

  vector<Logscale>deriv_copy(n_nodes+4, Logscale(0.0,0.0));
  for (size_t i = 0; i < n_hme; ++i ){
    llk = log_sum_log(llk, forward[i]);
    for (size_t j = 0; j < n_nodes+4; ++j){
      Logscale val;
      log_sum_log_sign(deriv_copy[j], forward_deriv[i][j], val);
      deriv_copy[j].logval = val.logval;
      deriv_copy[j].symbol = val.symbol;
    }
  }

  deriv = vector<double>(n_nodes+4, 0.0);
  for (size_t i = 0; i < deriv.size(); ++i) {
    if (i!=4)
      deriv[i] = exp(deriv_copy[i].logval-llk)*deriv_copy[i].symbol;
  }
  return llk;
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
    size_t minCpG = 50;

    // run mode flags
    bool VERBOSE = false;
    OptionParser opt_parse(strip_path(argv[0]), "test funcitonality of "
                           "phylo-methylome segmentation",
                           "<newick> <meth-tab>");
    opt_parse.add_opt("minCpG", 'm', "minimum observed #CpGs in a block"
                      "(default: 50)", false, minCpG);
    opt_parse.add_opt("verbose", 'v', "print more run info (default: false)",
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
    const string tree_file = leftover_args.front();
    const string meth_table_file = leftover_args.back();

    /******************** LOAD PHYLOGENETIC TREE ******************************/
    std::ifstream tree_in(tree_file.c_str());
    if (!tree_in)
      throw SMITHLABException("cannot read: " + tree_file);
    PhyloTreePreorder t;
    tree_in >> t;

    if (VERBOSE)
      cerr <<  t.tostring() << endl;
    vector<size_t> subtree_sizes;
    vector<size_t> node_degrees;
    vector<string> nodenames;
    vector<double> branches;
    t.get_subtree_sizes(subtree_sizes);
    t.get_branch_lengths(branches);
    get_degrees(subtree_sizes, node_degrees);
    if (!is_semi_binary(node_degrees))
      throw SMITHLABException("invalid tree structure");
    const size_t n_leaves =  leafsize(subtree_sizes);
    if (VERBOSE)
      cerr << "loaded tree (leaves=" << n_leaves << ")" << endl;

    /**************************************************************************/
    /******************** LOAD METHYLATION DATA *******************************/
    vector<Site> sites;
    vector<string> meth_table_species;
    vector<vector<double> > meth_prob_table;
    load_meth_table(meth_table_file, sites, meth_prob_table,
                    meth_table_species);
    if (VERBOSE)
      cerr << "loaded meth table (species="
           << meth_table_species.size() << ")" << endl;


    /**************************************************************************/
    /**********  ENSURE METH DATA AND TREE INFO IS IN SYNC ********************/
    if (meth_table_species.size() != n_leaves)
      throw SMITHLABException("inconsistent species counts");
    if (!has_same_species_order(t, meth_table_species))
      throw SMITHLABException("inconsistent species names or order");
    vector<size_t> species_in_order;
    subtree_sizes_to_leaves_preorder(subtree_sizes, species_in_order);
    for (size_t i = 0; i < species_in_order.size(); ++i) {
      cerr << meth_table_species[i] << "\t" << species_in_order[i] << endl;
    }

    if (VERBOSE)
      cerr << "Read in total " << sites.size() << " sites" << endl
           << "Separating sites" << endl;
    vector<size_t> reset_points;
    separate_regions(desert_size, minCpG, meth_prob_table,
                     sites, reset_points);
    if (VERBOSE)
      cerr << "After filtering: " << sites.size()
           << "in " << reset_points.size() << "blocks" << endl;

    const size_t n_nodes = subtree_sizes.size();
    const size_t n_hme = pow(2,n_nodes);
    vector<string> all_states;
    for (size_t i = 0; i < n_hme; ++i) {
      all_states.push_back(dec2binary(n_nodes, i));
    }

    //test transition probabilities.
    double pi0 = 0.1;
    double rate0 = 0.5;
    double g0 = 0.9;
    double g1 = 0.9;
    double llk = 0;

    vector<vector<vector<double> > > time_trans_mats;
    vector<vector<vector<vector<double> > > > combined_trans_mats;

    collect_transition_matrices(g0, g1, rate0, branches,
                                time_trans_mats, combined_trans_mats);

    vector<vector<vector<vector<double> > > > combined_trans_mats_drate;
    vector<vector<vector<vector<double> > > > combined_trans_mats_dg0;
    vector<vector<vector<vector<double> > > > combined_trans_mats_dg1;
    vector<vector<vector<vector<double> > > > combined_trans_mats_dT;

    cerr <<"[collect_transition_matrices_deriv]" << endl;
    collect_transition_matrices_deriv(g0, g1, rate0, branches, time_trans_mats, combined_trans_mats,
                                      combined_trans_mats_drate, combined_trans_mats_dg0,
                                      combined_trans_mats_dg1, combined_trans_mats_dT);

    vector<double> exp_branches;
    for (size_t i = 0; i < branches.size(); ++i) {
      exp_branches.push_back(exp(branches[i]));
    }

    bool FULL = false;
    if (FULL) { // Full transition
      vector<vector<double> > hme_log_tpm;
      bool DERIV = true;
      if (!DERIV) {
        cerr << "==============FULL TRANSITION========================" << endl;
        hme_trans_prob_mat(pi0, g0, g1, rate0, branches, subtree_sizes,
                           combined_trans_mats, hme_log_tpm);
        llk = llk_by_forward_algorithm(pi0, subtree_sizes, all_states,
                                       time_trans_mats, hme_log_tpm, meth_prob_table);
        cerr << "[Log-likelihood] " << llk << endl;
      } else { // with derivatives
        cerr << "==============FULL TRANSITION WITH DERIVATIVES========" << endl;
        vector<vector<double> > hme_trans_dg0;
        vector<vector<double> > hme_trans_dg1;
        vector<vector<double> > hme_trans_drate;
        vector<vector<vector<double> > > hme_trans_dT;
        hme_trans_prob_mat_deriv(g0, g1, rate0, branches, subtree_sizes, time_trans_mats,
                                 combined_trans_mats, combined_trans_mats_drate,
                                 combined_trans_mats_dg0, combined_trans_mats_dg1,
                                 combined_trans_mats_dT, hme_log_tpm, hme_trans_drate,
                                 hme_trans_dg0, hme_trans_dg1, hme_trans_dT);

        vector<double> deriv;
        llk = llk_by_forward_algorithm_deriv(pi0, rate0, subtree_sizes, exp_branches,
                                             all_states, time_trans_mats, hme_log_tpm,
                                             hme_trans_drate, hme_trans_dg0, hme_trans_dg1,
                                             hme_trans_dT, meth_prob_table, deriv);
        cerr << "[Log-likelihood] " << llk << endl;
        for (size_t i = 0; i < deriv.size(); ++i) {
          cerr << "d[" << i << "]" << deriv[i] << endl;
        }
      }
    }
    ///////////////////////////////////////////////
    bool RSTR =true;

    // Restricted neighbor transition
    if (RSTR) {
      vector<vector<double> > rstr_hme_log_tpm;
      vector<vector<size_t> > neighbors;
      vector<vector<size_t> > place_in_neighbor;

      cerr << "==============       RESTRICTED          ========"<< endl;
      bool DERIV = true;
      if (!DERIV){
        rstr_hme_trans_prob_mat(pi0, g0, g1, rate0, branches, subtree_sizes, all_states,
                                combined_trans_mats,
                                neighbors, place_in_neighbor, rstr_hme_log_tpm);

        llk = llk_by_rstr_forward_algorithm(pi0, subtree_sizes, all_states,
                                            time_trans_mats, neighbors,
                                            place_in_neighbor, rstr_hme_log_tpm,
                                            meth_prob_table);
        cerr << "[Log-likelihood] " << llk << endl;
      } else {
        vector<vector<double> > rstr_hme_tpm_drate;
        vector<vector<double> > rstr_hme_tpm_dg0;
        vector<vector<double> > rstr_hme_tpm_dg1;
        vector<vector<vector<double> > > rstr_hme_tpm_dT;
        cerr << "==============RESTRICTED WITH DERIVATIVES========" << endl;
        rstr_hme_trans_prob_mat_deriv(pi0, g0, g1, rate0, branches,
                                      subtree_sizes, time_trans_mats,
                                      all_states, combined_trans_mats, combined_trans_mats_drate,
                                      combined_trans_mats_dg0, combined_trans_mats_dg1,
                                      combined_trans_mats_dT,
                                      neighbors, place_in_neighbor,rstr_hme_log_tpm,
                                      rstr_hme_tpm_drate, rstr_hme_tpm_dg0,
                                      rstr_hme_tpm_dg1, rstr_hme_tpm_dT);
        vector<double> deriv;
        llk = llk_by_rstr_forward_algorithm_deriv(pi0, rate0, subtree_sizes, exp_branches,
                                                  all_states, time_trans_mats, neighbors,
                                                  place_in_neighbor, rstr_hme_log_tpm,
                                                  rstr_hme_tpm_drate, rstr_hme_tpm_dg0,
                                                  rstr_hme_tpm_dg1, rstr_hme_tpm_dT,
                                                  meth_prob_table, deriv);
        cerr << "[Log-likelihood] " << llk << endl;
        for (size_t i = 0; i < deriv.size(); ++i) {
          cerr << "d[" << i << "]" << deriv[i] << endl;
        }
      }
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
