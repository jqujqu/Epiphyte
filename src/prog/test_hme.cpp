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
#include <math.h>  //sqrt

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


static void
load_full_states(const string filename,
                 const size_t n_nodes,
                 vector<Site> &sites,
                 vector<string> &states) {

  std::ifstream in(filename.c_str());
  if (!in)
    throw SMITHLABException("cannot read:" + filename);
  string line;
  while (getline(in, line)) {
    if (line[0] == '#') { //skip header lines
      continue;
    } else {
      istringstream iss(line);
      string first_column;
      iss >> first_column;
      Site site;
      site.assign_from_methtab_first_column(first_column);
      sites.push_back(site);

      string state;
      iss >> state;
      if (state.length() != n_nodes)
        throw SMITHLABException("inconsistent state length: " + line);
      states.push_back(state);
    }
  }
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
            double &root_unmeth_prob, double &rate0,
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
  issp >> root_unmeth_prob >> rate0;
  if (VERBOSE)
    cerr << "Root_unmeth_prob=" << root_unmeth_prob
         << "\trate=" <<  rate0 << endl;
  if (root_unmeth_prob >= 1.0 || root_unmeth_prob <= 0.0 ||
      rate0 >= 1.0 || rate0 <= 0.0)
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
separate_regions(const bool VERBOSE,
                 const size_t desert_size,
                 const size_t min_site_per_block,
                 vector<vector<double> > &meth, vector<Site> &sites,
                 vector<size_t> &reset_points,
                 double &g0_est, double &g1_est) {
  const size_t nleaf = meth[0].size();
  const size_t totalsites = sites.size();

  //scan desert
  //uncovered site has prob value set to -1.0
  vector<bool> is_desert(totalsites, false);
  for (size_t i = 0; i < nleaf; ++i) {
    if (VERBOSE)
      cerr << "Processing sample" << i ;
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
    if (VERBOSE)
      cerr << "[Done]" << endl;
  }

  size_t good =0;
  for (size_t i = 0; i < totalsites; ++i) {
    if (!is_desert[i]) ++good;
  }
  if (VERBOSE)
    cerr << "Total" << good << " good sites" << endl;

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
        if (last_block_size < min_site_per_block) { //remove small block
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

  double g00= 0.0, g01=0.0, g10=0.0, g11=0.0;
  for (size_t i = 0; i < nleaf; ++i) {
    for (size_t j = 0; j < meth.size()-1;++j){
      if (meth[j][i] >= 0 && meth[j+1][i]>=0) {
        g00 += meth[j][i]*meth[j+1][i];
        g01 += meth[j][i]*(1.0-meth[j+1][i]);
        g10 += (1.0- meth[j][i])*meth[j+1][i];
        g11 += (1.0-meth[j][i])*(1.0-meth[j+1][i]);
      }
    }
  }
  g0_est = (g00+1)/(g01+g00+1);
  g1_est = (g11+1)/(g11+g10+1);
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


////////////////////////////////////////////////////////////////////////////////
//////////////////////      quantity units          ////////////////////////////
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
temporal_trans_prob_mat(const double T, const double rate0,
                        vector<vector<double> > &time_transition_matrix) {
  assert(rate0 > 0 && rate0 < 1 && T <1 && T>0);
  time_transition_matrix = vector<vector<double> >(2, vector<double>(2, 0.0));
  time_transition_matrix[0][0] = 1.0 - rate0*T;
  time_transition_matrix[0][1] = rate0*T;
  time_transition_matrix[1][0] = (1.0-rate0)*T;
  time_transition_matrix[1][1] = 1.0 - (1.0-rate0)*T;
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
      prev_anc_denom[prev][anc] = G[prev][0]*time_trans_mat[anc][0] +
        G[prev][1]*time_trans_mat[anc][1];
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



//T=1-exp(-branch)
static void
trans_deriv(const vector<vector<double> > &Q,
            const vector<vector<double> > &G,
            const double T,
            const size_t prev,
            const size_t anc,
            const size_t cur,
            const vector<vector<double> > &prev_anc_denom,
            const vector<vector<double> > &time_trans_mat,
            double &d_rate0, double &d_g0, double &d_g1, double &d_T) {

  const double denom = prev_anc_denom[prev][anc];

  double term = kronecker_delta(cur, 1)*2 - 1.0 - time_trans_mat[anc][cur]*(G[prev][1]-G[prev][0])/denom;
  d_rate0 = T*G[prev][cur]/denom*term;

  if (prev == 0) {
    d_g1 = 0;
    term = kronecker_delta(cur, 0)*2 - 1.0 - G[prev][cur]*(time_trans_mat[anc][0] - time_trans_mat[anc][1])/denom ;
    d_g0= time_trans_mat[anc][cur]/denom*term;
  } else {
    d_g0 = 0;
    term = kronecker_delta(cur, 1)*2 - 1.0 -
      G[prev][cur]*(time_trans_mat[anc][1] - time_trans_mat[anc][0])/denom;
    d_g1 = time_trans_mat[anc][cur]/denom*term;
  }
  term = Q[anc][cur] - time_trans_mat[anc][cur]*(G[prev][0]*Q[anc][0]+G[prev][1]*Q[anc][1])/denom;
  d_T = G[prev][cur]/denom*term;
}

// On a single branch
// T=1-exp(-branch)
static void
combined_trans_prob_mat_deriv(const double rate0,
                              const double g0, const double g1,
                              const double T,
                              const vector<vector<double> > &time_trans_mat,
                              vector<vector<vector<double> > > &combined_trans_mat, //prev x anc x cur
                              vector<vector<vector<double> > > &combined_trans_mat_drate,
                              vector<vector<vector<double> > > &combined_trans_mat_dg0,
                              vector<vector<vector<double> > > &combined_trans_mat_dg1,
                              vector<vector<vector<double> > > &combined_trans_mat_dT) {
  // deriv: (rate0, g0, g1, T)
  vector<vector<double> > G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][1] = g1;
  G[1][0] = 1.0 - g1;
  vector<vector<double> > Q(2,vector<double>(2,0.0));
  Q[0][0] = -1.0*rate0;
  Q[0][1] = rate0;
  Q[1][0] = 1.0 - rate0;
  Q[1][1] = rate0 - 1.0;

  vector<vector<double> > prev_anc_denom(2, vector<double>(2, 0.0));
  for (size_t prev = 0; prev < 2; ++ prev) {
    for (size_t anc = 0; anc < 2; ++ anc) {
      prev_anc_denom[prev][anc] =
        G[prev][0]*time_trans_mat[anc][0] + G[prev][1]*time_trans_mat[anc][1] ;
    }
  }

  combined_trans_mat =
    vector<vector<vector<double> > >(2, vector<vector<double> >(2,vector<double>(2, 0.0)));
  combined_trans_mat_drate = combined_trans_mat;
  combined_trans_mat_dg0 = combined_trans_mat;
  combined_trans_mat_dg1 = combined_trans_mat;
  combined_trans_mat_dT = combined_trans_mat;

  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      for (size_t cur = 0; cur < 2; ++cur) {
        combined_trans_mat[prev][anc][cur] =
          G[prev][cur]*time_trans_mat[anc][cur]/prev_anc_denom[prev][anc];
      }
    }
  }

  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      for (size_t cur = 0; cur < 2; ++ cur) {
        trans_deriv(Q, G, T, prev, anc, cur,
                    prev_anc_denom, time_trans_mat,
                    combined_trans_mat_drate[prev][anc][cur],
                    combined_trans_mat_dg0[prev][anc][cur],
                    combined_trans_mat_dg1[prev][anc][cur],
                    combined_trans_mat_dT[prev][anc][cur]);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//////////////////////       Units grouped       ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
static void
collect_transition_matrices(const double g0, const double g1,
                            const double rate0, const vector<double> &Ts,
                            vector<vector<vector<double> > > &time_trans_mats,
                            vector<vector<vector<vector<double> > > > &combined_trans_mats) {
  assert(g0 > 0 && g0 <1 && g1 > 0 && g1 <1);

  const size_t n_nodes = Ts.size();
  vector<vector<double> > mat2x2(2,vector<double>(2,0.0));
  time_trans_mats = vector<vector<vector<double> > >(n_nodes, mat2x2);
  vector<vector<vector<double> > > mat2x2x2(2, mat2x2);
  combined_trans_mats = vector<vector<vector<vector<double> > > >(n_nodes, mat2x2x2);

  for (size_t i = 1; i < n_nodes; ++i) {
    temporal_trans_prob_mat(Ts[i], rate0, time_trans_mats[i]);
    combined_trans_prob_mat(g0, g1, time_trans_mats[i], combined_trans_mats[i]);
  }
}


static void
collect_transition_matrices_deriv(const double rate0, const double g0, const double g1,
                                  const vector<double> &Ts,
                                  vector<vector<vector<double> > > &time_trans_mats,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_drate,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_dg0,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_dg1,
                                  vector<vector<vector<vector<double> > > > &combined_trans_mats_dT) {
  assert(g0 > 0 && g0 <1 && g1 > 0 && g1 <1);
  const size_t n_nodes = Ts.size();
  vector<vector<double> > mat2x2(2, vector<double>(2, 0.0));
  vector<vector<vector<double> > > mat2x2x2(2, mat2x2);

  time_trans_mats = vector<vector<vector<double> > > (n_nodes, mat2x2);
  vector<vector<vector<vector<double> > > > multi_mat2x2x2(n_nodes, mat2x2x2);
  combined_trans_mats =  multi_mat2x2x2; // node x prev x anc x cur
  combined_trans_mats_drate = multi_mat2x2x2;
  combined_trans_mats_dg0 = multi_mat2x2x2;
  combined_trans_mats_dg1 = multi_mat2x2x2;
  combined_trans_mats_dT = multi_mat2x2x2;

  for (size_t i = 1; i < n_nodes; ++i) {
    temporal_trans_prob_mat(Ts[i], rate0, time_trans_mats[i]);
    combined_trans_prob_mat_deriv(rate0, g0, g1, Ts[i], time_trans_mats[i],
                                  combined_trans_mats[i],
                                  combined_trans_mats_drate[i],
                                  combined_trans_mats_dg0[i],
                                  combined_trans_mats_dg1[i],
                                  combined_trans_mats_dT[i]);
  }
}

////////////////////////////////////////////////////////////////////////////////
//////////       All quantities for likelihood computation       ///////////////
////////////////////////////////////////////////////////////////////////////////
// probability and derivatives
static void
auxiliary(const vector<size_t> &subtree_sizes,
          const vector<string> &all_states,
          const vector<double> &params,
          vector<double> &start_log_probs,
          vector<vector<Logscale> >&start_log_prob_derivs,
          vector<vector<double> > &transition_log_probs,
          vector<vector<vector<Logscale> > >&transition_log_prob_derivs) {

  vector<vector<double> > Q(2,vector<double>(2,0.0));
  Q[0][0] = -1.0*params[1];
  Q[0][1] = params[1];
  Q[1][0] = 1.0 - params[1];
  Q[1][1] = params[1] - 1.0;

  vector<vector<double> > G(2,vector<double>(2,0.0));
  G[0][0] = params[2];
  G[0][1] = 1.0 - params[2];
  G[1][0] = 1.0 - params[3];
  G[1][1] = params[3];

  cerr << "In auxiliary" << endl;
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_hme = all_states.size();
  const vector<vector<double> > mat2x2(2, vector<double>(2,0.0));
  vector<vector<vector<double> > > time_trans_mats(n_nodes, mat2x2);

  // collect probabilities and derivatives by branch
  vector<vector<vector<vector<double> > > > combined_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats_drate;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg0;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg1;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dT;
  double rate0 = params[1];
  double g0 = params[2];
  double g1 = params[3];
  vector<double> Ts(params.begin()+4, params.end());
  collect_transition_matrices_deriv(rate0, g0, g1, Ts, time_trans_mats,
                                    combined_trans_mats,
                                    combined_trans_mats_drate,
                                    combined_trans_mats_dg0,
                                    combined_trans_mats_dg1,
                                    combined_trans_mats_dT);

  //single HME probs and derivs
  start_log_probs = vector<double>(n_hme, 0.0);
  const size_t n_params = params.size();
  start_log_prob_derivs = vector<vector<Logscale> >(n_hme, vector<Logscale>(n_params, Logscale(0.0, 0.0)));
  for (size_t i = 0; i < n_hme; ++i) {
    string states = all_states[i];
    start_log_probs[i] =
      (states[0] == '0') ? log(params[0]) : log(1.0-params[0]);
    // d_pi0 factor
    start_log_prob_derivs[i][0].logval =
      (states[0] == '0') ? - 1.0*log(params[0]) : -1.0*log(1.0-params[0]);
    start_log_prob_derivs[i][0].symbol = (states[0] == '0') ? 1.0 : -1.0;

    for (size_t j = 0; j < n_nodes; ++j) {
      size_t anc = (states[j] == '0')? 0 : 1;
      if (subtree_sizes[j] > 1) {
        size_t count = 1;
        Logscale inc(0.0, 0.0);
        while (count < subtree_sizes[j]) {
          size_t child = j + count;
          size_t cur = (states[child] == '0')? 0 : 1;

          // log-prob
          start_log_probs[i] += log(time_trans_mats[child][anc][cur]);

          //d_rate factor
          inc.logval = log(Ts[child]) - log(time_trans_mats[child][anc][cur]);
          inc.symbol = (cur == 1) ? 1.0 : -1.0;
          Logscale result(0.0, 0.0);
          log_sum_log_sign(start_log_prob_derivs[i][1], inc, result);
          start_log_prob_derivs[i][1].logval = result.logval;
          start_log_prob_derivs[i][1].symbol = result.symbol;

          // d_T factor
          start_log_prob_derivs[i][child+4].logval =
            log(abs(Q[anc][cur])) - log(time_trans_mats[child][anc][cur]);
          start_log_prob_derivs[i][child+4].symbol = sign(Q[anc][cur]);

          count += subtree_sizes[child];
        }
      }
    }

    // multiply deriv factors with prob to get derivs
    for (size_t j = 0; j < params.size(); ++j) {
      if (j < 2 || j > 4)
        start_log_prob_derivs[i][j].logval += start_log_probs[i];
    }

  }

  //transition and derivs
  transition_log_probs =
    vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  vector<Logscale> param_deriv_holder(n_params, Logscale(0.0, 0.0));
  transition_log_prob_derivs =
    vector<vector<vector<Logscale> > >(n_hme, vector<vector<Logscale> >(n_hme, param_deriv_holder));
  for (size_t i = 0; i < n_hme; ++i) {
    string prevstates = all_states[i];
    for (size_t j = 0; j < n_hme; ++j) {
      string curstates = all_states[j];
      size_t prev = (prevstates[0]=='0') ? 0 : 1;
      size_t cur = (curstates[0]=='0') ? 0 : 1;
      // transtion at root
      transition_log_probs[i][j] = log(G[prev][cur]);
      // transiton at each branch
      for (size_t k = 0; k < n_nodes; ++k) {
        size_t anc = (curstates[k] == '0')? 0 : 1;
        if (subtree_sizes[k] > 1) {
          size_t count = 1;
          while (count < subtree_sizes[k]) {
            size_t child = k + count;
            cur = (curstates[child] == '0')? 0 : 1;
            prev = (prevstates[child] == '0')? 0 : 1;
            double lp = log(combined_trans_mats[child][prev][anc][cur]);
            transition_log_probs[i][j] += lp;

            //d_pi0
            transition_log_prob_derivs[i][j][0].logval = 0.0;
            transition_log_prob_derivs[i][j][0].symbol = 0.0;

            //d_rate factor
            Logscale inc(log(abs(combined_trans_mats_drate[child][prev][anc][cur]))- lp,
                         sign(combined_trans_mats_drate[child][prev][anc][cur]));
            Logscale result;
            log_sum_log_sign(transition_log_prob_derivs[i][j][1], inc, result);
            transition_log_prob_derivs[i][j][1].logval = result.logval;
            transition_log_prob_derivs[i][j][1].symbol = result.symbol;

            //d_g0 factor
            if (prev==0) {
              inc.symbol = sign(combined_trans_mats_dg0[child][prev][anc][cur]);
              inc.logval = log(abs(combined_trans_mats_dg0[child][prev][anc][cur])) - lp;
              log_sum_log_sign(transition_log_prob_derivs[i][j][2],
                               inc, result);
              transition_log_prob_derivs[i][j][2].logval = result.logval;
              transition_log_prob_derivs[i][j][2].symbol = result.symbol;
            }

            //d_g1 factor
            if (prev==1) {
              inc.symbol = sign(combined_trans_mats_dg1[child][prev][anc][cur]);
              inc.logval =
                log(abs(combined_trans_mats_dg1[child][prev][anc][cur])) - lp;
              log_sum_log_sign(transition_log_prob_derivs[i][j][3],
                               inc, result);
              transition_log_prob_derivs[i][j][3].logval = result.logval;
              transition_log_prob_derivs[i][j][3].symbol = result.symbol;
            }

            //d_T factor
            inc.symbol = sign(combined_trans_mats_dT[child][prev][anc][cur]);
            inc.logval = log(abs(combined_trans_mats_dT[child][prev][anc][cur])) - lp;
            log_sum_log_sign(transition_log_prob_derivs[i][j][child+4],
                             inc, result);
            transition_log_prob_derivs[i][j][child+4].logval = result.logval;
            transition_log_prob_derivs[i][j][child+4].symbol = result.symbol;
            count += subtree_sizes[child];
          }
        }
      }

      // multiply factor to transition probability (rate, T, g0, g1)
      transition_log_prob_derivs[i][j][1].logval += transition_log_probs[i][j];
      if (transition_log_prob_derivs[i][j][2].symbol != 0)
        transition_log_prob_derivs[i][j][2].logval += transition_log_probs[i][j];
      if (transition_log_prob_derivs[i][j][3].symbol != 0)
        transition_log_prob_derivs[i][j][3].logval += transition_log_probs[i][j];
      for (size_t k = 1; k < n_nodes; ++k) {
        transition_log_prob_derivs[i][j][k+4].logval += transition_log_probs[i][j];
      }
    }
  }
}


//probability only
static void
auxiliary(const vector<size_t> &subtree_sizes,
          const vector<string> &all_states,
          const vector<double> &params,
          vector<double> &start_log_probs,
          vector<vector<double> > &transition_log_probs) {

  vector<vector<double> > Q(2,vector<double>(2,0.0));
  Q[0][0] = -1.0*params[1];
  Q[0][1] = params[1];
  Q[1][0] = 1.0 - params[1];
  Q[1][1] = params[1] - 1.0;

  vector<vector<double> > G(2,vector<double>(2,0.0));
  G[0][0] = params[2];
  G[0][1] = 1.0 - params[2];
  G[1][0] = 1.0 - params[3];
  G[1][1] = params[3];

  cerr << "In auxiliary prob" << endl;
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_hme = all_states.size();
  const vector<vector<double> > mat2x2(2, vector<double>(2,0.0));
  vector<vector<vector<double> > > time_trans_mats(n_nodes, mat2x2);

  // collect probabilities and derivatives by branch
  vector<vector<vector<vector<double> > > > combined_trans_mats;
  double rate0 = params[1];
  double g0 = params[2];
  double g1 = params[3];
  vector<double> Ts(params.begin()+4, params.end());
  collect_transition_matrices(rate0, g0, g1, Ts, time_trans_mats,
                              combined_trans_mats);

  //single HME probs
  start_log_probs = vector<double>(n_hme, 0.0);
  for (size_t i = 0; i < n_hme; ++i) {
    string states = all_states[i];
    start_log_probs[i] =
      (states[0] == '0') ? log(params[0]) : log(1.0-params[0]);
    for (size_t j = 0; j < n_nodes; ++j) {
      size_t anc = (states[j] == '0') ? 0 : 1;
      if (subtree_sizes[j] > 1) {
        size_t count = 1;
        Logscale inc(0.0, 0.0);
        while (count < subtree_sizes[j]) {
          size_t child = j + count;
          size_t cur = (states[child] == '0')? 0 : 1;
          start_log_probs[i] += log(time_trans_mats[child][anc][cur]);
          count += subtree_sizes[child];
        }
      }
    }
  }

  //transition
  transition_log_probs =
    vector<vector<double> >(n_hme, vector<double>(n_hme, 0.0));
  for (size_t i = 0; i < n_hme; ++i) {
    string prevstates = all_states[i];
    for (size_t j = 0; j < n_hme; ++j) {
      string curstates = all_states[j];
      size_t prev = (prevstates[0]=='0') ? 0 : 1;
      size_t cur = (curstates[0]=='0') ? 0 : 1;
      transition_log_probs[i][j] = log(G[prev][cur]);
      for (size_t k = 0; k < n_nodes; ++k) {
        size_t anc = (curstates[k] == '0')? 0 : 1;
        if (subtree_sizes[k] > 1) {
          size_t count = 1;
          while (count < subtree_sizes[k]) {
            size_t child = k + count;
            cur = (curstates[child] == '0')? 0 : 1;
            prev = (prevstates[child] == '0')? 0 : 1;
            double lp = log(combined_trans_mats[child][prev][anc][cur]);
            transition_log_probs[i][j] += lp;
            count += subtree_sizes[child];
          }
        }
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
///////////    transformation between states and states-id     /////////////////
////////////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////////////////
/////////////////         Neighbor restriction       ///////////////////////////
////////////////////////////////////////////////////////////////////////////////

// restricted transition states
// [0]: self; A&B differences can be explained by a subtree
static void
get_rstr_state_ids(const vector<size_t> &subtree_sizes,
                   const vector<string> all_states,
                   vector<vector<size_t> > &rstr_state_ids,
                   vector<vector<size_t> > &place_in_neighbor) {

  size_t n_hme = all_states.size();
  rstr_state_ids = vector<vector<size_t> >(n_hme, vector<size_t>(2*n_hme, 0));
  for (size_t k = 0; k < n_hme; ++k)
    rstr_state_ids[k].resize(0);

  for (size_t k = 0; k < n_hme; ++ k) {
    string states = all_states[k];
    const size_t n_nodes = subtree_sizes.size();
    rstr_state_ids[k].push_back(k); //self_transition
    for (size_t i = 0; i < n_nodes; ++i) {
      size_t start = i;
      for (size_t s = 0; s < 2; ++s) {
        string nb_states = states;
        string replacestr = (s == 0) ? "0": "1";
        for (size_t j = 0; j < subtree_sizes[i]; ++j) {
          nb_states.replace(start+j, 1, replacestr);
        }
        size_t neighbor_id = binary2dec(nb_states);
        rstr_state_ids[k].push_back(neighbor_id);
        rstr_state_ids[neighbor_id].push_back(k);
      }
    }
  }

  for (size_t k = 0; k < n_hme; ++k) {
    std::sort(rstr_state_ids[k].begin(), rstr_state_ids[k].end());
    vector<size_t>::iterator it;
    it = std::unique(rstr_state_ids[k].begin(), rstr_state_ids[k].end());
    rstr_state_ids[k].resize(std::distance(rstr_state_ids[k].begin(), it) );
    assert (rstr_state_ids[k].size() <= n_hme);
  }

  // my rank in each of my neighbor's neighbor list
  for (size_t i = 0; i < n_hme; ++i) {
    vector<size_t> places;
    for (size_t j = 0; j < rstr_state_ids[i].size(); ++j) {
      const size_t neighbor = rstr_state_ids[i][j];
      const size_t myid = distance(rstr_state_ids[neighbor].begin(),
                                   find(rstr_state_ids[neighbor].begin(),
                                        rstr_state_ids[neighbor].end(), i));
      places.push_back(myid);
    }
    place_in_neighbor.push_back(places);
  }
}


//neighboring transitions: self + change of state at one node and all nodes in its subtree.
static void
rstr_hme_trans_prob_mat(const vector<size_t> &subtree_sizes,
                        const vector<string> &all_states,
                        const vector<double> &start_log_probs,
                        const vector<vector<double> > &transition_log_probs,
                        vector<vector<size_t> > &neighbors,
                        vector<vector<size_t> > &place_in_neighbor,
                        vector<vector<double > >&rstr_trans_log_probs) {

  const size_t n_hme = all_states.size();
  get_rstr_state_ids(subtree_sizes, all_states, neighbors, place_in_neighbor);

  for (size_t i = 0; i < n_hme; ++i) { //prev state: i
    double log_row_sum = 0;
    size_t n_neighbor = neighbors[i].size();

    // resize transitions for each state
    rstr_trans_log_probs[i].resize(n_neighbor);

    // Compute transition probabilities
    for (size_t j = 0; j < n_neighbor; ++j) {
      size_t id = neighbors[i][j]; // cur state: neighbors[i][j]
      log_row_sum = log_sum_log(log_row_sum, transition_log_probs[i][id]);
    }

    // Normalize transition probabiliteis
    for (size_t j = 0; j < n_neighbor; ++j) {
      size_t id = neighbors[i][j];
      rstr_trans_log_probs[i][j] = transition_log_probs[i][id] - log_row_sum;
    }
  }

  for (size_t i = 0; i < n_hme; ++i) {
    cerr << all_states[i] << "\t" << neighbors[i].size() << "neighbors\t";
    for (size_t j = 0; j < neighbors[i].size(); ++j ) {
      cerr << exp(transition_log_probs[i][j]) << "\t";
    }
    cerr << endl;
  }
}


static void
rstr_hme_trans_prob_mat_deriv(const vector<size_t> &subtree_sizes,
                              const vector<string> &all_states,
                              const vector<double> &start_log_probs,
                              const vector<vector<Logscale> > &start_log_prob_derivs,
                              const vector<vector<double> > &trans_log_probs,
                              const vector<vector<vector<Logscale> > > &trans_log_prob_derivs,
                              vector<vector<size_t> > &neighbors,
                              vector<vector<size_t> > &place_in_neighbor,
                              vector<vector<double > > &rstr_trans_log_probs,
                              vector<vector<vector<Logscale> > > &rstr_trans_log_prob_derivs) {

  const size_t n_nodes = subtree_sizes.size();
  const size_t n_hme = all_states.size();
  get_rstr_state_ids(subtree_sizes, all_states, neighbors, place_in_neighbor);

  // initialize a big matrix
  rstr_trans_log_probs = vector<vector<double> > (n_hme, vector<double>(n_nodes+4, 0.0));
  rstr_trans_log_prob_derivs = vector<vector<vector<Logscale> > >(n_hme, vector<vector<Logscale> >(n_hme, vector<Logscale>(n_nodes+4, Logscale(0.0, 0.0))));

  for (size_t i = 0; i < n_hme; ++i) { //prev state: i
    double log_row_sum = 0;
    size_t n_neighbor = neighbors[i].size();
    // resize transitions for each state
    rstr_trans_log_probs[i].resize(n_neighbor);
    rstr_trans_log_prob_derivs[i].resize(n_neighbor);

    // Compute restricted transition probabilities
    for (size_t j = 0; j < n_neighbor; ++j) {
      size_t id = neighbors[i][j]; // cur state: neighbors[i][j]
      log_row_sum = log_sum_log(log_row_sum, trans_log_probs[i][id]);
    }

    // Normalize transition probabiliteis
    for (size_t j = 0; j < n_neighbor; ++j) {
      size_t id = neighbors[i][j];
      rstr_trans_log_probs[i][j] = trans_log_probs[i][id] - log_row_sum;
    }

    //Compute restricted transition prob derivatives
    vector<Logscale> rowsum_derivs_log(n_nodes+4, Logscale(0.0, 0.0));
    for (size_t j = 0; j < n_neighbor; ++j) {
      size_t id = neighbors[i][j];
      for (size_t k = 0; k < n_nodes+4; ++k) {
        Logscale result(0.0, 0.0);
        log_sum_log_sign(rowsum_derivs_log[k],
                         trans_log_prob_derivs[i][id][k], result);
        rowsum_derivs_log[k].logval = result.logval;
        rowsum_derivs_log[k].symbol = result.symbol;
      }
    }

    for (size_t j = 0; j < n_neighbor; ++j ) {
      size_t id = neighbors[i][j];
      for (size_t k = 1; k < n_nodes+4; ++k) {
        if (k !=4) {
          Logscale part1(trans_log_prob_derivs[i][id][k].logval - log_row_sum,
                         transition_log_prob_derivs[i][id][k].symbol);
          Logscale part2((trans_log_probs[i][id]-log_row_sum*2) +
                         rowsum_derivs_log[k].logval,
                         rowsum_derivs_log[k].symbol);
          log_sum_log_sign(part1, part2, rstr_trans_log_prob_derivs[i][j][k]);
        }
      }
    }
  }
}




static void
forward_variables(const size_t start, const size_t end,
                  const vector<size_t> &subtree_sizes,
                  const vector<string> &all_states,
                  const vector<double> &start_log_probs,
                  const vector<vector<double> > &transition_log_probs,
                  const vector<vector<double> > &meth_prob_table,
                  const size_t pos,
                  vector<double> &forward) {
  assert (pos >=start && pos < end);
  const size_t n_hme = all_states.size();
  forward = vector<double> (n_hme, 0.0);
  for (size_t i = 0; i < n_hme; ++i) {
    forward[i] = start_log_probs[i] +
      evaluate_emission(subtree_sizes, all_states[i], meth_prob_table[pos]);
  }

  for (size_t k = start; k < pos; ++k) {
    vector<double> ff;
    ff.swap(forward);
    forward = vector<double> (n_hme, 0.0);
    for (size_t j = 0; j < n_hme; ++j) {
      for (size_t i = 0; i < n_hme; ++i) {
        forward[j] = log_sum_log(forward[j],
                                 ff[i] + transition_log_probs[i][j]);
      }
      forward[j] += evaluate_emission(subtree_sizes, all_states[j],
                                      meth_prob_table[pos]);
    }
  }
  if (pos %1000 == 0) cerr << ".";
}

static void
forward_variables_deriv(const size_t start, const size_t end,
                        const vector<size_t> &subtree_sizes,
                        const vector<string> &all_states,
                        const vector<double> &start_log_probs,
                        const vector<vector<Logscale> >&start_log_prob_derivs,
                        const vector<vector<double> > &transition_log_probs,
                        const vector<vector<vector<Logscale> > >&transition_log_prob_derivs,
                        const vector<vector<double> > &meth_prob_table,
                        const size_t pos,
                        vector<double> &forward,
                        vector<vector<Logscale> > &log_forward_deriv) {

  assert (pos >= start && pos < end);

  const size_t n_hme = all_states.size();
  const size_t n_nodes = subtree_sizes.size();
  forward = vector<double> (n_hme, 0.0);
  log_forward_deriv =
    vector<vector<Logscale> >(n_hme,
                              vector<Logscale>(n_nodes+4, Logscale(0.0,0.0)));

  for (size_t i = 0; i < n_hme; ++i) {
    forward[i] =  start_log_probs[i] +
      evaluate_emission(subtree_sizes, all_states[i], meth_prob_table[pos]);
    for (size_t j = 0; j < 4+n_nodes; ++j) {
      log_forward_deriv[i][j].logval = start_log_prob_derivs[i][j].logval;
      log_forward_deriv[i][j].symbol = start_log_prob_derivs[i][j].symbol;
    }
  }

  for (size_t k =start; k < pos; ++k) {
    if (k %1000 == 0) cerr << ".";
    vector<double> ff;
    vector<vector<Logscale> > ffdd;
    ff.swap(forward);
    ffdd.swap(log_forward_deriv);
    forward = vector<double> (n_hme, 0.0);
    log_forward_deriv =
      vector<vector<Logscale> >(n_hme, vector<Logscale>(n_nodes+4,
                                                        Logscale(0.0,0.0)));

    for (size_t j = 0; j < n_hme; ++ j) { //current site
      for (size_t i = 0; i < n_hme; ++ i) { //previous site
        double log_trans_prob = transition_log_probs[i][j];
        forward[j] = log_sum_log(forward[j], ff[i] + log_trans_prob);

        Logscale first, second, inc, result;
        for (size_t k = 0; k < 4+n_nodes; ++k) {
        first.symbol = ffdd[i][k].symbol;
        if (first.symbol != 0.0)
          first.logval = ffdd[i][k].logval + log_trans_prob;
        second.symbol = transition_log_prob_derivs[i][j][k].symbol;
        if (second.symbol != 0.0)
          second.logval = ff[i] + transition_log_prob_derivs[i][j][k].logval;

        log_sum_log_sign(first, second, inc);
        log_sum_log_sign(log_forward_deriv[j][k], inc, result);
        log_forward_deriv[j][k].logval = result.logval;
        log_forward_deriv[j][k].symbol = result.symbol;
        }
      }
      double emi_log_prob = evaluate_emission(subtree_sizes, all_states[j],
                                              meth_prob_table[k]);
      forward[j] += emi_log_prob;
      for (size_t k = 0; k < n_nodes+4; ++k) {
        if (k!=4)
          log_forward_deriv[j][k].logval += emi_log_prob;
      }
    }
  }
}


double
llk_by_forward_algorithm(const vector<size_t> &reset_points,
                         const vector<size_t> &subtree_sizes,
                         const vector<string> &all_states,
                         const vector<double> &params,
                         const vector<vector<double> > &meth_prob_table) {

  vector<double> start_log_probs;
  vector<vector<double> > transition_log_probs;
  auxiliary(subtree_sizes, all_states, params,
            start_log_probs, transition_log_probs);

  const size_t n_hme = all_states.size();
  double llk = 0;
  for (size_t k = 0; k < reset_points.size()-1; ++k) {
    const size_t start = reset_points[k];
    const size_t end = reset_points[k+1];
    double block_llk = 0;
    vector<double> forward;

    forward_variables(start, end, subtree_sizes, all_states,
                      start_log_probs, transition_log_probs,
                      meth_prob_table, end-1, forward);

    for (size_t i = 0; i < n_hme; ++i ){
      block_llk = log_sum_log(block_llk, forward[i]);
    }
    llk += block_llk;
  }
  return llk;
}


double
llk_by_forward_algorithm_deriv(const vector<size_t> &reset_points,
                               const vector<size_t> &subtree_sizes,
                               const vector<string> &all_states,
                               const vector<double> &params,
                               const vector<vector<double> > &meth_prob_table,
                               vector<double> &deriv) {
  const size_t n_nodes = subtree_sizes.size();

  vector<double> start_log_probs;
  vector<vector<double> > transition_log_probs;
  vector<vector<Logscale> > start_log_prob_derivs;
  vector<vector<vector<Logscale> > > transition_log_prob_derivs;

  auxiliary(subtree_sizes, all_states, params,
            start_log_probs, start_log_prob_derivs,
            transition_log_probs, transition_log_prob_derivs);

  const size_t n_hme = all_states.size();
  double llk = 0;
  deriv = vector<double>(n_nodes+4, 0.0);
  for (size_t k = 0; k < reset_points.size()-1; ++k) {
    const size_t start = reset_points[k];
    const size_t end = reset_points[k+1];
    double block_llk = 0;
    vector<double> forward;
    vector<vector<Logscale> > forward_deriv;

  forward_variables_deriv(start, end, subtree_sizes, all_states,
                            start_log_probs, start_log_prob_derivs,
                            transition_log_probs, transition_log_prob_derivs,
                            meth_prob_table, end-1,
                            forward, forward_deriv) ;

    vector<Logscale>deriv_copy(n_nodes+4, Logscale(0.0, 0.0));
    Logscale val;
    for (size_t i = 0; i < n_hme; ++i ){
      block_llk = log_sum_log(block_llk, forward[i]);
      for (size_t j = 0; j < deriv.size(); ++j){
        log_sum_log_sign(deriv_copy[j], forward_deriv[i][j],val);
        deriv_copy[j].logval = val.logval;
        deriv_copy[j].symbol = val.symbol;
      }
    }

    llk += block_llk;

    for (size_t i = 0; i < deriv.size(); ++i) {
      if (i!=4)
        deriv[i] += exp(deriv_copy[i].logval-block_llk)*deriv_copy[i].symbol;
    }
  }
  return llk;
}


static void
forward_variables_rstr(const size_t start, const size_t end,
                       const vector<size_t> &subtree_sizes,
                       const vector<string> &all_states,
                       const vector<double> &start_log_probs,
                       const vector<vector<size_t> > &neighbors,
                       const vector<vector<size_t> > &place_in_neighbor,
                       const vector<vector<double > > &rstr_trans_log_probs,
                       const vector<vector<double> > &meth_prob_table,
                       const size_t pos,
                       vector<double> &forward) {

  assert(pos >=start && pos < end);
  forward.clear();

  // begin
  const size_t n_hme = all_states.size();
  forward = vector<double>(n_hme, 0.0);
  for (size_t i = 0; i < n_hme; ++i) {
    forward[i] = start_log_probs[i] +
      evaluate_emission(subtree_sizes, all_states[i], meth_prob_table[start]);
  }
  for (size_t k = start; k < pos; ++k) {
    if ((k+1)%1000 == 0) cerr << ".";
    vector<double> ff;
    ff.swap(forward);
    forward = vector<double>(n_hme, 0.0);

    for (size_t j = 0; j < n_hme; ++ j) { //current site
      for (size_t i = 0; i < neighbors[i].size(); ++ i) { //previous site
        const size_t prev_state_id = neighbors[j][i];
        const size_t idx_in_neighbor = place_in_neighbor[j][i];
        double ele = ff[prev_state_id] +
          rstr_trans_log_probs[prev_state_id][idx_in_neighbor] ;
        forward[j] = log_sum_log(forward[j], ele);
      }
      forward[j] += evaluate_emission(subtree_sizes, all_states[j],
                                      meth_prob_table[k]);
    }
  }
}



static void
forward_variables_rstr_deriv(const size_t start, const size_t end,
                             const vector<size_t> &subtree_sizes,
                             const vector<string> &all_states,
                             const vector<double> &start_log_probs,
                             const vector<vector<Logscale> > &start_log_prob_derivs,
                             const vector<vector<size_t> > &neighbors,
                             const vector<vector<size_t> > &place_in_neighbor,
                             const vector<vector<double > > &rstr_trans_log_probs,
                             const vector<vector<vector<Logscale> > > &rstr_trans_log_prob_derivs,
                             const vector<vector<double> > &meth_prob_table,
                             const size_t pos,
                             vector<double> &forward,
                             vector<vector<Logscale> > &log_forward_deriv) {

  assert(pos >=start && pos < end);
  const size_t n_hme = all_states.size();
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_params = 4+n_nodes;
  forward = vector<double> (n_hme, 0.0);
  log_forward_deriv = vector<vector<Logscale> >(n_hme, vector<Logscale>(n_params, Logscale(0.0, 0.0)));

  // begin
  for (size_t i = 0; i < n_hme; ++i) {
    double log_emission = evaluate_emission(subtree_sizes, all_states[i], meth_prob_table[start]);
    forward[i] = start_log_probs[i] + log_emission;

    for (size_t n = 0; n < n_params; ++n ) {
      log_forward_deriv[i][n].logval = start_log_prob_derivs[i][n].logval + log_emission;
      log_forward_deriv[i][n].symbol = start_log_prob_derivs[i][n].symbol;
    }
  }

  for (size_t  k = start; k < pos; ++k){
    if ((k+1)%100 == 0) cerr << ".";
    vector<double> ff;
    vector<vector<Logscale> > ffdd;
    ff.swap(forward);
    ffdd.swap(log_forward_deriv);

    forward = vector<double> (n_hme, 0.0);
    log_forward_deriv = vector<vector<Logscale> >(n_hme, vector<Logscale>(n_params, Logscale(0.0, 0.0)));

    for (size_t j = 0; j < n_hme; ++ j) { //current site
      size_t n_neighbor = neighbors[j].size();
      for (size_t i = 0; i < n_neighbor; ++ i) { //previous site
        const size_t prev_state_id = neighbors[j][i];
        const size_t idx_in_neighbor = place_in_neighbor[j][i];
        const double log_trans_prob =  rstr_trans_log_probs[prev_state_id][idx_in_neighbor];
        forward[j] = log_sum_log(forward[j], ff[prev_state_id] + log_trans_prob);

        //update derivatives here
        Logscale first, second, inc, result;
        for (size_t n = 0; n < n_params; ++n) {
          first.symbol = ffdd[prev_state_id][n].symbol;
          if (first.symbol != 0.0)
            first.logval = ffdd[prev_state_id][n].logval + log_trans_prob;
          second.symbol = rstr_trans_log_prob_derivs[prev_state_id][idx_in_neighbor][n].symbol;
          if (second.symbol != 0.0)
            second.logval = ff[prev_state_id] +
              rstr_trans_log_prob_derivs[prev_state_id][idx_in_neighbor][n].logval;

          log_sum_log_sign(first, second, inc);
          log_sum_log_sign(log_forward_deriv[j][n], inc, result);
          log_forward_deriv[j][n].logval = result.logval;
          log_forward_deriv[j][n].symbol = result.symbol;
        }
      }
      double log_emission = evaluate_emission(subtree_sizes, all_states[j],
                                              meth_prob_table[k]);
      forward[j] += log_emission;
      for (size_t n = 0; n < n_params; ++n) {
        if (n!=4) {
          log_forward_deriv[j][n].logval += log_emission;
        }
      }
    }
  }
}


double
llk_by_rstr_forward_algorithm(const vector<size_t> &reset_points,
                              const vector<size_t> &subtree_sizes,
                              const vector<string> &all_states,
                              const vector<double> &params,
                              const vector<vector<double> > &meth_prob_table) {

  vector<double> start_log_probs;
  vector<vector<double> > transition_log_probs;
  auxiliary(subtree_sizes, all_states, params,
            start_log_probs, transition_log_probs);

  vector<vector<size_t> > neighbors;
  vector<vector<size_t> > place_in_neighbor;
  vector<vector<double > > rstr_trans_log_probs;
  rstr_hme_trans_prob_mat(subtree_sizes, all_states,
                          start_log_probs, transition_log_probs,
                          neighbors, place_in_neighbor, rstr_trans_log_probs);

  const size_t n_hme = all_states.size();

  double llk = 0;
  for (size_t k = 0; k < reset_points.size()-1; ++k) {
    const size_t start = reset_points[k];
    const size_t end = reset_points[k+1];
    double block_llk = 0;
    vector<double> forward;
    forward_variables_rstr(start, end, subtree_sizes, all_states,
                           start_log_probs, neighbors, place_in_neighbor,
                           rstr_trans_log_probs, meth_prob_table,
                           end -1, forward);

    for (size_t i = 0; i < n_hme; ++i ){
      block_llk = log_sum_log(block_llk, forward[i]);
    }
    llk += block_llk;
  }
  return llk;
}


double
llk_by_rstr_forward_algorithm_deriv(const vector<size_t> &reset_points,
                                    const vector<size_t> &subtree_sizes,
                                    const vector<string> &all_states,
                                    const vector<double> &params,
                                    const vector<vector<double> > &meth_prob_table,
                                    vector<double> &deriv) {
  const size_t n_params = params.size();

  vector<double> start_log_probs;
  vector<vector<Logscale> > start_log_prob_derivs;
  vector<vector<double> > transition_log_probs;
  vector<vector<vector<Logscale> > > transition_log_prob_derivs;
  auxiliary(subtree_sizes, all_states, params,
            start_log_probs, start_log_prob_derivs,
            transition_log_probs, transition_log_prob_derivs);

  vector<vector<size_t> > neighbors;
  vector<vector<size_t> > place_in_neighbor;
  vector<vector<double > > rstr_trans_log_probs;
  vector<vector<vector<Logscale> > > rstr_trans_log_prob_derivs;
  rstr_hme_trans_prob_mat_deriv(subtree_sizes, all_states, start_log_probs,
                                start_log_prob_derivs, transition_log_probs,
                                transition_log_prob_derivs, neighbors,
                                place_in_neighbor, rstr_trans_log_probs,
                                rstr_trans_log_prob_derivs);

  const size_t n_hme = all_states.size();
  deriv = vector<double>(n_params, 0.0);
  double llk = 0;

  for (size_t k = 0; k < reset_points.size()-1; ++k) {
    const size_t start = reset_points[k];
    const size_t end = reset_points[k+1];

    double block_llk = 0;
    vector<double> forward;
    vector<vector<Logscale> > forward_deriv;
    forward_variables_rstr_deriv(start, end, subtree_sizes, all_states,
                                 start_log_probs, start_log_prob_derivs,
                                 neighbors, place_in_neighbor,
                                 rstr_trans_log_probs,
                                 rstr_trans_log_prob_derivs,
                                 meth_prob_table, end-1, forward,
                                 forward_deriv);

    vector<Logscale>deriv_copy(n_params, Logscale(0.0,0.0));
    for (size_t i = 0; i < n_hme; ++i ){
      block_llk = log_sum_log(block_llk, forward[i]);
      for (size_t j = 0; j < n_params; ++j){
        Logscale val;
        log_sum_log_sign(deriv_copy[j], forward_deriv[i][j], val);
        deriv_copy[j].logval = val.logval;
        deriv_copy[j].symbol = val.symbol;
      }
    }

    for (size_t i = 0; i < n_params; ++i) {
      if (i!=4) {
        deriv[i] += exp(deriv_copy[i].logval-block_llk)*deriv_copy[i].symbol;
      }
    }

    llk += block_llk;
  }
  return llk;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////         OPTIMIZATION          /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool
optimize_iteration(const bool FULL, const double TOL,
                   const vector<vector<double> > &meth_prob_table,
                   const vector<size_t> &reset_points,
                   const vector<size_t> &subtree_sizes,
                   const vector<string> &all_states,
                   const vector<double> &params, //pi0, rate0, g0, g1, 1.0-exp(-branches)
                   const double &llk,
                   const vector<double> &deriv,
                   bool &CONVERGED,
                   vector<double> &newparams,
                   double &newllk,
                   vector<double> &newderiv) {
  const size_t n_nodes = subtree_sizes.size();
  bool SUCCESS = false;

  double norm = deriv[0]*deriv[0] + deriv[1]*deriv[1];
  norm += deriv[2]*deriv[2] + deriv[3]*deriv[3];
  for (size_t i = 5; i < deriv.size(); ++i) {
    norm += deriv[i]*deriv[i];
  }
  double r = 1.0/sqrt(norm); //maximization (-1.0 if minimization)

  // Find proper starting factor
  const size_t n_params = n_nodes + 4;
  for (size_t i = 0; i < n_params; ++i) {
    if (i!=4) {
      double candidate_param = params[i] + deriv[i]*r;
      while (candidate_param < TOL || candidate_param > 1.0-TOL) {
        r = r/2;
        candidate_param = params[i] + deriv[i]*r;
     }
    }
  }

  newparams = params;
  // Test param point
  for (size_t i = 0; i < n_params; ++i) {
    if (i!=4)
      newparams[i] = params[i] + r*deriv[i];
  }

  // Reduce factor untill improvement or convergence
  while (!SUCCESS && !CONVERGED ) {
    r = r/2;
    cerr << r << endl;

    for (size_t i = 0; i < deriv.size(); ++i)
      cerr << "d[" << i << "]=" << deriv[i] << "\t";
    cerr << endl;


    if (r*sqrt(norm) <= TOL) {
      CONVERGED = true;
      return SUCCESS;
    }

    for (size_t i = 0; i < n_params; ++i) {
      if (i!=4)
        newparams[i] = params[i] + r*deriv[i];
    }

    vector<double> Ts;
    Ts.assign(newparams.begin()+4, newparams.end());
    if (FULL) {
      newllk = llk_by_forward_algorithm_deriv(reset_points, subtree_sizes,
                                              all_states, params,
                                              meth_prob_table, newderiv);
    } else {
      newllk = llk_by_rstr_forward_algorithm_deriv(reset_points, subtree_sizes,
                                                   all_states, newparams,
                                                   meth_prob_table, newderiv);
    }

    if (newllk > llk)
      SUCCESS = true;
  }

  return SUCCESS;
}



static void
optimize(const bool FULL,
         const double tolerance,
         const vector<vector<double> > &meth_prob_table,
         const vector<size_t> &reset_points,
         const vector<size_t> &subtree_sizes,
         const vector<string> &all_states,
         const vector<double> &start_params) {

  vector<double> deriv;
  double llk;
  if (FULL) {
    llk = llk_by_forward_algorithm_deriv(reset_points, subtree_sizes,
                                         all_states, start_params,
                                         meth_prob_table, deriv);
  } else {
    llk = llk_by_rstr_forward_algorithm_deriv(reset_points, subtree_sizes,
                                              all_states, start_params,
                                              meth_prob_table, deriv);
  }
  bool CONVERGED = false;
  bool SUCCESS = true;
  size_t iteration = 0;
  cerr << "Iteration " << iteration << "\t llk=" << llk << "\t";
  vector<double> params = start_params;
  for (size_t i = 0; i < params.size(); ++i) {
    if (i < 4) cerr << params[i] << "\t";
    if (i > 4) cerr << -1.0*log(1.0-params[i]) << "\t";
  }
  cerr << endl;

  while (!CONVERGED && SUCCESS) {
    vector<double> newparams;
    vector<double> newderiv;
    double newllk;
    SUCCESS = optimize_iteration(FULL, tolerance, meth_prob_table,
                                 reset_points, subtree_sizes, all_states,
                                 params, llk, deriv, CONVERGED, newparams,
                                 newllk, newderiv);
    if (SUCCESS) {
      params = newparams;
      llk = newllk;
      deriv = newderiv;
    }
    ++iteration;
    cerr << endl << "Iteration " << iteration << "\t llk=" << llk << "\t";
    for (size_t i = 0; i < params.size(); ++i) {
      if (i < 4) cerr << params[i] << "\t";
      if (i > 4) cerr << -1.0*log(1.0-params[i]) << "\t";
    }
    cerr << endl;
  }
}


static void
fill_leaf_prob(const vector<Site> &sites,
               const size_t desert_size,
               vector<vector<double> > &meth_prob_table) {
  size_t nsites = meth_prob_table.size();
  size_t n_leaf = meth_prob_table[0].size();
  for (size_t i = 0; i < n_leaf; ++i) {
    size_t prev = 0;
    size_t next = 0;
    size_t count = 0;
    for (size_t j = 0; j < nsites; ++j) {
      if (meth_prob_table[j][i] >= 0) {
        prev = j;
        next = j;
      } else {
        if (j > next) {
          while (next < nsites-1 && meth_prob_table[next][i]<0) ++next;
        }
        size_t d1 = distance_between_sites(sites[j], sites[prev]);
        size_t d2 = distance_between_sites(sites[j], sites[next]);
        if (j > prev && j < next && d1 < desert_size && d2 < desert_size) {
          meth_prob_table[j][i] = (meth_prob_table[prev][i]*d1 +
                                   meth_prob_table[next][i]*d2)/(d1+d2);
          ++count;
        } else if (prev < j && d1 < desert_size) {
          meth_prob_table[j][i] = (0.5 + meth_prob_table[prev][i])/2;
          ++count;
        } else {
          meth_prob_table[j][i] = (0.5 + meth_prob_table[next][i])/2;
          ++count;
        }
      }
    }
    cerr << "Filled " << count << " missing sites in leaf" << i << endl;
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////         COMPLETE DATA         /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
separate_regions(const size_t desert_size,
                 vector<string> &states, vector<Site> &sites,
                 vector<size_t> &reset_points) {

  const size_t totalsites = sites.size();

  reset_points.push_back(0);
  Site last_site = sites[0];
  for (size_t i = 1; i < totalsites; ++i) {
    if (distance_between_sites(last_site, sites[i]) < desert_size) {
    } else {
      reset_points.push_back(i);
    }
    last_site = sites[i];
  }
  reset_points.push_back(totalsites);
}


static void
loglik_complete_tree(const vector<string> &states,
                     const vector<size_t> &reset_points,
                     const vector<size_t> &subtree_sizes,
                     const vector<double> &params,
                     double &llk, vector<Logscale> &param_derivs) {

  vector<double> start_log_probs;
  vector<vector<Logscale> > start_log_prob_derivs;
  vector<vector<double> > trans_log_probs;
  vector<vector<vector<Logscale> > > trans_log_prob_derivs;

  vector<string> all_states;
  size_t n_hme = pow(2, subtree_sizes.size());
  size_t n_nodes = subtree_sizes.size();
  for (size_t i = 0; i < n_hme; ++i) {
    all_states.push_back(dec2binary(n_nodes, i));
  }
  auxiliary(subtree_sizes, all_states, params,
            start_log_probs, start_log_prob_derivs,
            trans_log_probs, trans_log_prob_derivs);

  const size_t n_params = params.size();
  param_derivs= vector<Logscale>(n_params,Logscale(0.0, 0.0));
  llk = 0.0;
  for (size_t i = 0; i < reset_points.size()-1; ++i) {
    const size_t start = reset_points[i];
    const size_t end = reset_points[i+1];
    size_t cur_id = binary2dec(states[start]);
    double lp = start_log_probs[cur_id];
    llk += lp;

    for (size_t j = 0; j < n_params; ++j) {
      Logscale deriv_inc(start_log_prob_derivs[cur_id][j].logval-lp,
                         start_log_prob_derivs[cur_id][j].symbol);
      Logscale sum;
      log_sum_log_sign(param_derivs[j], deriv_inc, sum);
      param_derivs[j].logval = sum.logval;
      param_derivs[j].symbol = sum.symbol;
    }

    for (size_t pos = start+1; pos < end; ++pos) {
      size_t prev_id = cur_id;
      cur_id = binary2dec(states[pos]);
      lp = trans_log_probs[prev_id][cur_id];
      llk += lp;
      for (size_t j = 0; j < n_params; ++j) {
        Logscale deriv_inc(trans_log_prob_derivs[prev_id][cur_id][j].logval-lp,
                           trans_log_prob_derivs[prev_id][cur_id][j].symbol);
        Logscale sum;
        log_sum_log_sign(param_derivs[j], deriv_inc, sum);
        param_derivs[j].logval = sum.logval;
        param_derivs[j].symbol = sum.symbol;
      }
    }
  }
}


static void
leaf_to_tree_prob(const vector<size_t> &subtree_sizes,
                  const vector<vector<double> > &meth_prob_table,
                  const vector<double> &params,
                  vector<vector<double> > &tree_prob_table) {



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
    bool COMPLETE = false;
    OptionParser opt_parse(strip_path(argv[0]), "test funcitonality of "
                           "phylo-methylome segmentation",
                           "<newick> <meth-tab>");
    opt_parse.add_opt("minCpG", 'm', "minimum observed #CpGs in a block"
                      "(default: 50)", false, minCpG);
    opt_parse.add_opt("complete", 'c', "complete observations",
                      false, COMPLETE);
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

    if (COMPLETE) {
      vector<string> states;
      load_full_states(meth_table_file, subtree_sizes.size(), sites, states);
      vector<size_t> reset_points;

      separate_regions(desert_size, states, sites, reset_points);
      if (VERBOSE)
        cerr << "Total" << reset_points.size() -1 << "blocks" << endl;
      //test transition probabilities.
      double pi0 = 0.1;
      double rate0 = 0.5;
      double g0 = 0.8;
      double g1 = 0.9;

      vector<double> Ts;
      for (size_t i = 0; i < branches.size(); ++i) {
        Ts.push_back(1.0-1.0/exp(branches[i]));
      }

      const size_t n_nodes = subtree_sizes.size();
      vector<double> start_param(n_nodes+4, 0.0);
      start_param[0] = pi0; start_param[1] = rate0;
      start_param[2] = g0; start_param[3] = g1;
      for (size_t i = 1; i < n_nodes; ++i) {
        start_param[4+i] = Ts[i];
      }

      vector<Logscale> param_derivs;
      double llk;
      loglik_complete_tree(states, reset_points, subtree_sizes,
                           start_param, llk, param_derivs);
      cerr << llk << endl;
      for (size_t i = 0; i < param_derivs.size(); ++i)
        cerr << exp(param_derivs[i].logval)*param_derivs[i].symbol << "\t" ;

    }

    if (!COMPLETE) {
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


      /**************************************************************************/
      /*************************  SEPARATE BY DESERTS ***************************/
      if (VERBOSE)
        cerr << "Read in total " << sites.size() << " sites" << endl
             << "Separating sites" << endl;
      vector<size_t> reset_points;
      double g0_est, g1_est;
      separate_regions(VERBOSE, desert_size, minCpG, meth_prob_table,
                       sites, reset_points, g0_est, g1_est);

      if (VERBOSE)
        cerr << "After filtering: " << sites.size()
             << "in " << reset_points.size() -1 << "blocks" << endl;

      fill_leaf_prob(sites, desert_size, meth_prob_table);
      if (VERBOSE)
        cerr << "Filled probabilities at missing sites" << endl;


      /**************************************************************************/
      /*************************  STARTING PARAMETERS ***************************/
      const size_t n_nodes = subtree_sizes.size();
      const size_t n_hme = pow(2,n_nodes);
      vector<string> all_states;
      for (size_t i = 0; i < n_hme; ++i) {
        all_states.push_back(dec2binary(n_nodes, i));
      }

      //test transition probabilities.
      double pi0 = 0.1;
      double rate0 = 0.39;
      double g0 = 0.72;
      double g1 = 0.78;
      double llk = 0;

      vector<double> Ts;
      for (size_t i = 0; i < branches.size(); ++i) {
        Ts.push_back(1.0-1.0/exp(branches[i]));
      }

      vector<double> start_param(n_nodes+4, 0.0);
      start_param[0] = pi0; start_param[1] = rate0;
      start_param[2] = g0; start_param[3] = g1;
      for (size_t i = 1; i < n_nodes; ++i) {
        start_param[4+i] = Ts[i];
      }

      /**************************************************************************/
      /*************************  OPTIMIZATION **********************************/
      bool FULL = false;
      bool DERIV = true;
      if (!DERIV) {
        vector<double> deriv;
        if (FULL) {
          llk = llk_by_forward_algorithm_deriv(reset_points, subtree_sizes,
                                               all_states, start_param,
                                               meth_prob_table, deriv);
        } else {
          llk = llk_by_rstr_forward_algorithm_deriv(reset_points, subtree_sizes,
                                                    all_states, start_param,
                                                    meth_prob_table, deriv);
        }
        cerr << "[Log-likelihood] " << llk << endl;
        for (size_t k = 0; k < n_nodes+4; ++k) {
          if (k!=4)
            cerr << deriv[k] << "\t";
          cerr << endl;
        }
      } else {
        if (FULL)
          cerr << "============== FULL TRANSITION===================" << endl;
        else
          cerr << "================ RESTRICTED ====================="<< endl;

        double tolerance = 1e-10;
        optimize(FULL, tolerance, meth_prob_table, reset_points, subtree_sizes,
                 all_states, start_param);
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