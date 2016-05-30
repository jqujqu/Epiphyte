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
  return (val == 0)? 0.0 : ( (val>0)? 1.0 : -1.0 );
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
  vector<pair<size_t, size_t> > counts;
  size_t meth_count = 0, unmeth_count = 0;
  while (iss >> meth_count) {
    if (iss >> unmeth_count) {
      counts.push_back(make_pair(meth_count, unmeth_count));
    } else {
      throw SMITHLABException("bad line format: " + iss.str());
    }
  }
  meth.push_back(counts);
}



void
parse_table_line(istringstream &iss,
                 vector<vector<double> > &meth) {
  double val;
  vector<double> meth_fields;
  while (iss >> val) {
    meth_fields.push_back(val);
  }
  meth.push_back(meth_fields);
  for (size_t i = 0; i < meth.back().size(); ++ i) {
    if ((meth.back()[i] < 0 && meth.back()[i]!=-1.0) || meth.back()[i] > 1.0 )
      throw SMITHLABException("bad line content: [" + iss.str() +
                              "], expecting hypomethylation probabilities (Missing values should be -1)");
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
read_params(const string &paramfile,
            double &root_unmeth_prob, double &rate0,
            double &g0, double &g1,
            PhyloTreePreorder &t) {
  std::ifstream in(paramfile.c_str());
  if (!in)
    throw SMITHLABException("cannot read: " + paramfile);

  string line;
  getline(in, line);
  istringstream iss(line);
  iss >> t;

  getline(in, line);
  istringstream issp_1(line);
  issp_1 >> root_unmeth_prob >> rate0;
  getline(in, line);
  istringstream issp_2(line);
  issp_2 >> g0 >> g1;
  if (root_unmeth_prob >= 1.0 || root_unmeth_prob <= 0.0 ||
      rate0 >= 1.0 || rate0 <= 0.0 || g0 >=1.0 || g0 <= 0.0 ||
      g1 >= 1.0 || g1 <= 0.0)
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
fill_leaf_prob(const bool VERBOSE,
               const vector<Site> &sites,
               const size_t desert_size,
               vector<vector<double> > &hypo_prob_table,
               double &pi0_est) {
  size_t nsites = hypo_prob_table.size();
  size_t n_leaf = hypo_prob_table[0].size();
  pi0_est = 0.0;
  for (size_t i = 0; i < n_leaf; ++i) {
    size_t prev = 0;
    size_t next = 0;
    size_t count = 0;
    double mean_hypo_prob = 0.0;
    size_t obs = 0;
    for (size_t j = 0; j < nsites; ++j) {
      if (hypo_prob_table[j][i] >= 0 ) {
        mean_hypo_prob += hypo_prob_table[j][i] ;
        ++obs;
      }
    }
    mean_hypo_prob = mean_hypo_prob/obs;
    pi0_est += mean_hypo_prob;
    if (VERBOSE) cerr << mean_hypo_prob << endl;

    for (size_t j = 0; j < nsites; ++j) {
      if (hypo_prob_table[j][i] >= 0) {
        prev = j;
        next = j;
      } else {
        if (j > next) {
          while (next < nsites-1 && hypo_prob_table[next][i]<0) ++next;
        }
        size_t d1 = distance_between_sites(sites[j], sites[prev]);
        size_t d2 = distance_between_sites(sites[j], sites[next]);
        if (j > prev && j < next && d1 < desert_size && d2 < desert_size) {
          hypo_prob_table[j][i] = (hypo_prob_table[prev][i]*d2 +
                                   hypo_prob_table[next][i]*d1)/(d1+d2);
          ++count;
        } else if (prev < j && d1 < desert_size) {
          hypo_prob_table[j][i] = (mean_hypo_prob*d1 + hypo_prob_table[prev][i]*(desert_size - d1))/desert_size;
          assert(hypo_prob_table[j][i] >= 0);
          ++count;
        } else if (j < next && d2 < desert_size){
          hypo_prob_table[j][i] = (mean_hypo_prob*d2 + hypo_prob_table[next][i]*(desert_size - d2))/desert_size;
          assert(hypo_prob_table[j][i] >= 0);
          ++count;
        } else {
          hypo_prob_table[j][i] = mean_hypo_prob;
        }
      }
    }
    if (VERBOSE) cerr << "Filled " << count << " missing sites in leaf" << i << endl;
  }
  pi0_est = pi0_est/n_leaf;
}



static void
separate_regions(const bool VERBOSE,
                 const size_t desert_size,
                 const size_t min_site_per_block,
                 vector<vector<double> > &meth, vector<Site> &sites,
                 vector<size_t> &reset_points,
                 double &g0_est, double &g1_est, double &pi0_est) {
  const size_t nleaf = meth[0].size();
  const size_t totalsites = sites.size();

  vector<vector<double> > meth_aug = meth;
  fill_leaf_prob(VERBOSE, sites, desert_size, meth_aug, pi0_est);

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

  size_t good = 0;
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
      meth_copy.push_back(meth_aug[i]);
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

  if (last_block_size < min_site_per_block) { //remove small block
    size_t start = reset_points.back();
    sites_copy.erase(sites_copy.begin()+start, sites_copy.end() );
    meth_copy.erase(meth_copy.begin()+start, meth_copy.end());
    reset_points.pop_back();
    last_block_size =
      reset_points.empty()? 0 : meth_copy.size()-reset_points.back();
  } else {
    reset_points.push_back(sites_copy.size());
  }

  meth.swap(meth_copy);
  sites.swap(sites_copy);

  double g00= 0.0, g01=0.0, g10=0.0, g11=0.0;
  for (size_t i = 0; i < nleaf; ++i) {
    for (size_t j = 0; j < meth.size()-1; ++j) {
      if (meth[j][i] >= 0 && meth[j+1][i] >= 0) {
        g00 += meth[j][i]*meth[j+1][i];
        g01 += meth[j][i]*(1.0-meth[j+1][i]);
        g10 += (1.0- meth[j][i])*meth[j+1][i];
        g11 += (1.0-meth[j][i])*(1.0-meth[j+1][i]);
      }
    }
  }
  g0_est = (g00 + 1.0)/(g01 + g00 + 1.0);
  g1_est = (g11 + 1.0)/(g11 + g10 + 1.0);
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


// P(T)[0,0]: hypo -> hypo prob
// P(T)[1,1]: hyper -> hyper prob
static void
temporal_trans_prob_mat(const double T, const double rate0,
                        vector<vector<double> > &time_transition_matrix) {
  assert(rate0 > 0 && rate0 < 1 && T < 1 && T > 0);
  time_transition_matrix = vector<vector<double> >(2, vector<double>(2, 0.0));
  time_transition_matrix[0][0] = 1.0 - rate0*T;
  time_transition_matrix[0][1] = rate0*T;
  time_transition_matrix[1][0] = (1.0 - rate0)*T;
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
  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      prev_anc_denom[prev][anc] = G[prev][0]*time_trans_mat[anc][0] +
        G[prev][1]*time_trans_mat[anc][1];
    }
  }

  combined_trans_mat =
    vector<vector<vector<double> > >(2, vector<vector<double> >(2,vector<double>(2,0.0)));
  for (size_t prev = 0; prev < 2; ++prev) {
    for (size_t anc = 0; anc < 2; ++anc) {
      for (size_t cur = 0; cur < 2; ++cur) {
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
    term = kronecker_delta(cur, 1)*2 - 1.0 - G[prev][cur]*(time_trans_mat[anc][1] - time_trans_mat[anc][0])/denom;
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
collect_transition_matrices(const double rate0, const double g0, const double g1,
                            const vector<double> &Ts,
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
///////////    transformation between states and states-id     /////////////////
////////////////////////////////////////////////////////////////////////////////

string
dec2binary(const size_t len, size_t n) {
  assert (n < pow(2, len));
  string result;
  do result.push_back( '0' + (n & 1) );
  while (n >>= 1);
  while (result.length() < len) {
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
////////////////////////////////////////////////////////////////////////////////
////////////////////      DATA AUGMENTATION        /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// uninitialized nodes have value -1.0
// assuming all leaves are initialized already
static void
init_internal_prob(const vector<size_t> &subtree_sizes,
                   const size_t node_id,
                   vector<double> &node_probs) {
  const double tol = 1e-5;
  size_t count = 1;
  // Recursion
  size_t n_desc_leaf = 0;
  double sum_desc_leaf_prob = 0;
  while (count < subtree_sizes[node_id]) {
    const size_t child_id = node_id + count;
    if (node_probs[child_id] < 0) {
      init_internal_prob(subtree_sizes, child_id, node_probs);
    }
    count += subtree_sizes[child_id];
  }

  double direction = 0.0;
  for (size_t i = 1; i < subtree_sizes[node_id]; ++i) {
    if (subtree_sizes[i+node_id]==1) {
      ++n_desc_leaf;
      sum_desc_leaf_prob += node_probs[i+node_id];
      direction += (node_probs[i+node_id] > 0.5) ? 1.0 : -1.0;
    }
  }

  if (node_id!=0 && abs(direction) < n_desc_leaf) {
    double sum_extra_leaf_hypo_prob = 0.0;
    size_t n_extra_leaves = 0;
    // children disagree, find majority state of leaf species outside this subtree
    for(size_t i = 0; i < subtree_sizes.size(); ++i) {
      if (subtree_sizes[i] == 1 &&
          (i < node_id || i >= node_id+subtree_sizes[node_id])) {
        sum_extra_leaf_hypo_prob += node_probs[i];
        ++n_extra_leaves;
      }
    }
    sum_desc_leaf_prob += sum_extra_leaf_hypo_prob/(n_extra_leaves + tol);
    ++n_desc_leaf;
  }

  node_probs[node_id] = min(max(tol, sum_desc_leaf_prob/n_desc_leaf), 1.0 - tol);
}


static void
leaf_to_tree_prob(const vector<size_t> &subtree_sizes,
                  const vector<vector<double> > &meth_prob_table,
                  vector<vector<double> > &tree_prob_table) {
  const double tol = 1e-5;
  const size_t n_sites = meth_prob_table.size();
  const size_t n_nodes = subtree_sizes.size();
  const size_t n_leaves = meth_prob_table[0].size();
  vector<size_t> leaves_preorder;
  subtree_sizes_to_leaves_preorder(subtree_sizes, leaves_preorder);

  // copy leaves first
  for (size_t i = 0; i < n_sites; ++i) {
    for (size_t j = 0; j < n_leaves; ++j) {
      tree_prob_table[i][leaves_preorder[j]] = min(max(tol, meth_prob_table[i][j]), 1.0-tol);
    }
  }
  // initialize internal probs
  for (size_t i = 0; i < n_sites; ++i) {
    init_internal_prob(subtree_sizes, 0, tree_prob_table[i]);
  }
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////         APPROX POSTERIOR        /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
posterior_start(const bool LEAF,
                const vector<size_t> &subtree_sizes,
                const vector<vector<double> >&G,
                const vector<vector<vector<double> > > &time_trans_mats,
                const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                const size_t pos,
                const size_t node_id,
                const size_t parent_id,
                vector<bool> &updated,
                vector<vector<double> > &meth_prob_table,
                double &diff) {
  //vector<bool> updated(subtree_sizes.size(), false);
  double p_cur = 0.0;
  double mar;
  if (subtree_sizes[node_id] == 1) { //leaf
    if (LEAF) {
      double p_parent_next = meth_prob_table[pos+1][parent_id];
      double p_parent_cur = meth_prob_table[pos][parent_id];
      double p_next = meth_prob_table[pos+1][node_id];

      for (size_t next = 0; next < 2; ++next) {
        double next_p = (next == 0) ? p_next : 1.0 - p_next;
        for (size_t par_next = 0; par_next <2; ++par_next) {
          double par_next_p = (par_next == 0) ? p_parent_next : 1.0- p_parent_next;
          for (size_t par_cur = 0; par_cur < 2; ++ par_cur) {
            double par_cur_p = (par_cur == 0) ? p_parent_cur : 1.0- p_parent_cur;
            mar = next_p*par_next_p*par_cur_p;
            double p_cur_u =  time_trans_mats[node_id][par_cur][0] *
              combined_trans_mats[node_id][0][par_next][next];
            double p_cur_m =  time_trans_mats[node_id][par_cur][1] *
              combined_trans_mats[node_id][1][par_next][next];
            double p_cur_mb = p_cur_u/(p_cur_u + p_cur_m);
            p_cur += p_cur_mb*mar;
          }
        }
      }
      diff += abs(meth_prob_table[pos][node_id] - p_cur);
      meth_prob_table[pos][node_id] = p_cur;
    }
  } else {
    //check if all children are updated
    vector<size_t> children;
    size_t count = 1;
    while (count < subtree_sizes[node_id]) {
      size_t child_id = node_id + count;
      children.push_back(child_id);
      if (!updated[child_id]) {
        posterior_start(LEAF, subtree_sizes, G, time_trans_mats, combined_trans_mats,
                        pos, child_id, node_id, updated, meth_prob_table, diff);
      }
      count += subtree_sizes[child_id];
    }

    if (node_id == 0) { //root node  //root has 3 children!
      double hypo_prob_next = meth_prob_table[pos+1][node_id];
      double hypo_prob_child_0 = meth_prob_table[pos][children[0]];
      double hypo_prob_child_1 = meth_prob_table[pos][children[1]];
      double hypo_prob_child_2 = meth_prob_table[pos][children[2]];

      for (size_t c0 = 0; c0 < 2; ++c0) {
        double c0_p = (c0 == 0) ? hypo_prob_child_0 : 1.0 - hypo_prob_child_0;
        for (size_t c1 = 0; c1 < 2; ++c1) {
          double c1_p = (c1 == 0) ? hypo_prob_child_1 : 1.0 - hypo_prob_child_1;
          for (size_t c2 = 0; c2 < 2; ++c2) {
            double c2_p = (c2 == 0) ? hypo_prob_child_2 : 1.0 - hypo_prob_child_2;
            for (size_t next = 0; next < 2; ++ next) {
              double next_p = (next == 0) ?  hypo_prob_next : 1.0 - hypo_prob_next;
              mar = c0_p*c1_p*c2_p*next_p;
              double p_cur_u = time_trans_mats[children[0]][0][c0]*
                time_trans_mats[children[1]][0][c1]*
                time_trans_mats[children[2]][0][c2]*G[0][next];
              double p_cur_m = time_trans_mats[children[0]][1][c0]*
                time_trans_mats[children[1]][1][c1]*
                time_trans_mats[children[2]][1][c2]*G[1][next];
              p_cur += p_cur_u/(p_cur_m + p_cur_u)*mar;
            }
          }
        }
      }
      diff += abs(meth_prob_table[pos][node_id] - p_cur);
      meth_prob_table[pos][node_id] = p_cur;
    } else { // internal node // two children
      double hypo_prob_child_0 = meth_prob_table[pos][children[0]];
      double hypo_prob_child_1 = meth_prob_table[pos][children[1]];
      double hypo_prob_next = meth_prob_table[pos+1][node_id];
      double hypo_prob_par_cur = meth_prob_table[pos][parent_id];
      double hypo_prob_par_next = meth_prob_table[pos+1][parent_id];

      for (size_t c0 = 0; c0 < 2; ++c0) {
        double c0_p = (c0 == 0) ? hypo_prob_child_0 : 1.0 - hypo_prob_child_0;
        for (size_t c1 = 0; c1 < 2; ++c1) {
          double c1_p = (c1 == 0) ? hypo_prob_child_1 : 1.0 - hypo_prob_child_1;
          for (size_t next = 0; next < 2; ++next) {
            double next_p = (next == 0) ? hypo_prob_next : 1.0 - hypo_prob_next;
            for (size_t par_cur = 0; par_cur < 2; ++par_cur) {
              double par_cur_p =
                (par_cur == 0) ? hypo_prob_par_cur : 1.0 - hypo_prob_par_cur;
              for (size_t par_next = 0; par_next < 2; ++par_next) {
                double par_next_p =
                  (par_next == 0) ? hypo_prob_par_next : 1.0 - hypo_prob_par_next;
                mar = c0_p*c1_p*next_p*par_cur_p*par_next_p;
                double p_cur_u = time_trans_mats[node_id][par_cur][0]*
                  time_trans_mats[children[0]][0][c0]*
                  time_trans_mats[children[1]][0][c1]*
                  combined_trans_mats[node_id][0][par_next][next];
                double p_cur_m = time_trans_mats[node_id][par_cur][1]*
                  time_trans_mats[children[0]][1][c0]*
                  time_trans_mats[children[1]][1][c1]*
                  combined_trans_mats[node_id][1][par_next][next];
                p_cur += p_cur_u/(p_cur_u + p_cur_m)*mar;
              }
            }
          }
        }
      }
    }
    diff += abs(meth_prob_table[pos][node_id] - p_cur);
    meth_prob_table[pos][node_id] = p_cur;
  }
  updated[node_id] = true;
}


static void
posterior_middle(const bool LEAF,
                 const vector<size_t> &subtree_sizes,
                 const vector<vector<double> >&G,
                 const vector<vector<vector<double> > > &time_trans_mats,
                 const vector<vector<vector<vector<double> > > > &combined_trans_mats,
                 const size_t pos,
                 const size_t node_id,
                 const size_t parent_id,
                 vector<bool> &updated,
                 vector<vector<double> > &meth_prob_table,
                 double &diff) {

  double mar;
  double p_cur = 0.0;
  if (subtree_sizes[node_id] == 1) { // leaf
    if (LEAF) {
      double p_prev = meth_prob_table[pos-1][parent_id];
      double p_parent_next = meth_prob_table[pos+1][parent_id];
      double p_parent_cur = meth_prob_table[pos][parent_id];
      double p_next = meth_prob_table[pos+1][node_id];

      for (size_t prev = 0; prev < 2; ++prev) {
        double prev_p = (prev == 0) ? p_prev: 1.0 - p_prev;
        for (size_t next = 0; next < 2; ++next) {
          double next_p = (next == 0) ? p_next : 1.0 - p_next;
          for (size_t par_next = 0; par_next <2; ++par_next) {
            double par_next_p = (par_next == 0) ? p_parent_next : 1.0- p_parent_next;
            for (size_t par_cur = 0; par_cur < 2; ++ par_cur) {
              double par_cur_p = (par_cur == 0) ? p_parent_cur : 1.0- p_parent_cur;
              mar = prev_p*next_p*par_next_p*par_cur_p;
              double p_cur_u =  combined_trans_mats[node_id][prev][par_cur][0] *
                combined_trans_mats[node_id][0][par_next][next];
              double p_cur_m =  combined_trans_mats[node_id][prev][par_cur][1] *
                combined_trans_mats[node_id][1][par_next][next];
              double p_cur_mb = p_cur_u/(p_cur_u + p_cur_m);
              p_cur += p_cur_mb*mar;
            }
          }
        }
      }
      diff += abs(meth_prob_table[pos][node_id] - p_cur);
      meth_prob_table[pos][node_id] = p_cur;
    }
  } else {
    //check if all children are updated
    vector<size_t> children;
    size_t count = 1;
    while (count < subtree_sizes[node_id]) {
      size_t child_id = node_id + count;
      children.push_back(child_id);
      if (!updated[child_id]) {
        posterior_middle(LEAF, subtree_sizes, G, time_trans_mats, combined_trans_mats,
                         pos, child_id, node_id, updated, meth_prob_table, diff);
      }
      count += subtree_sizes[child_id];
    }
    if (node_id == 0) { //root node  //root has 3 children!
      double hypo_prob_prev = meth_prob_table[pos-1][node_id];
      double hypo_prob_next = meth_prob_table[pos+1][node_id];
      double hypo_prob_child_0 = meth_prob_table[pos][children[0]];
      double hypo_prob_child_1 = meth_prob_table[pos][children[1]];
      double hypo_prob_child_2 = meth_prob_table[pos][children[2]];
      double hypo_prob_prev_child_0 = meth_prob_table[pos-1][children[0]];
      double hypo_prob_prev_child_1 = meth_prob_table[pos-1][children[1]];
      double hypo_prob_prev_child_2 = meth_prob_table[pos-1][children[2]];

      for (size_t prev = 0; prev < 2; ++prev) {
        double prev_p = (prev==0) ? hypo_prob_prev : 1.0 - hypo_prob_prev;
        for (size_t prev_c0 = 0; prev_c0 < 2; ++ prev_c0) {
          double prev_c0_p =
            (prev_c0==0) ? hypo_prob_prev_child_0 : 1.0 - hypo_prob_prev_child_0;
          for (size_t prev_c1 = 0; prev_c1 < 2; ++ prev_c1) {
            double prev_c1_p =
              (prev_c1==0) ? hypo_prob_prev_child_1 : 1.0 - hypo_prob_prev_child_1;
            for (size_t prev_c2 =0; prev_c2 < 2; ++ prev_c2) {
              double prev_c2_p =
                (prev_c2==0) ? hypo_prob_prev_child_2 : 1.0 - hypo_prob_prev_child_2;
              for (size_t c0 = 0; c0 < 2; ++c0) {
                double c0_p = (c0 == 0) ? hypo_prob_child_0 : 1.0 - hypo_prob_child_0;
                for (size_t c1 = 0; c1 < 2; ++c1) {
                  double c1_p = (c1 == 0) ? hypo_prob_child_1 : 1.0 - hypo_prob_child_1;
                  for (size_t c2 = 0; c2 < 2; ++c2) {
                    double c2_p = (c2 == 0) ? hypo_prob_child_2 : 1.0 - hypo_prob_child_2;
                    for (size_t next = 0; next < 2; ++ next) {
                      double next_p = (next == 0) ?  hypo_prob_next : 1.0 - hypo_prob_next;
                      mar = prev_p*prev_c0_p*prev_c1_p*prev_c2_p*c0_p*c1_p*c2_p*next_p;
                      double p_cur_u = combined_trans_mats[children[0]][prev][0][c0]*
                        combined_trans_mats[children[1]][prev][0][c1]*
                        combined_trans_mats[children[2]][prev][0][c2]*
                        G[0][next]*G[prev][0];
                      double p_cur_m = combined_trans_mats[children[0]][prev][1][c0]*
                        combined_trans_mats[children[1]][prev][1][c1]*
                        combined_trans_mats[children[2]][prev][1][c2]*
                        G[1][next]*G[prev][1];
                      p_cur += p_cur_u/(p_cur_m + p_cur_u)*mar;
                    }
                  }
                }
              }
            }
          }
        }
      }
      diff += abs(meth_prob_table[pos][node_id] - p_cur);
      meth_prob_table[pos][node_id] = p_cur;
    } else { // internal node // two children
      double hypo_prob_prev_child_0 = meth_prob_table[pos-1][children[0]];
      double hypo_prob_prev_child_1 = meth_prob_table[pos-1][children[1]];
      double hypo_prob_child_0 = meth_prob_table[pos][children[0]];
      double hypo_prob_child_1 = meth_prob_table[pos][children[1]];
      double hypo_prob_prev = meth_prob_table[pos-1][node_id];
      double hypo_prob_next = meth_prob_table[pos+1][node_id];
      double hypo_prob_par_cur = meth_prob_table[pos][parent_id];
      double hypo_prob_par_next = meth_prob_table[pos+1][parent_id];

      for (size_t prev_c0 = 0; prev_c0 < 2; ++ prev_c0) {
        double prev_c0_p = (prev_c0 == 0) ? hypo_prob_prev_child_0 : 1.0 - hypo_prob_prev_child_0;
        for (size_t prev_c1 = 0; prev_c1 < 2; ++ prev_c1) {
          double prev_c1_p = (prev_c1 == 0) ? hypo_prob_prev_child_1 : 1.0 - hypo_prob_prev_child_1;
          for (size_t prev = 0; prev < 2; ++ prev) {
            double prev_p = (prev == 0) ? hypo_prob_prev : 1.0 - hypo_prob_prev;
            for (size_t c0 = 0; c0 < 2; ++c0) {
              double c0_p = (c0 == 0) ? hypo_prob_child_0 : 1.0 - hypo_prob_child_0;
              for (size_t c1 = 0; c1 < 2; ++c1) {
                double c1_p = (c1 == 0) ? hypo_prob_child_1 : 1.0 - hypo_prob_child_1;
                for (size_t next = 0; next < 2; ++next) {
                  double next_p = (next == 0) ? hypo_prob_next : 1.0 - hypo_prob_next;
                  for (size_t par_cur = 0; par_cur < 2; ++par_cur) {
                    double par_cur_p = (par_cur == 0) ? hypo_prob_par_cur : 1.0 - hypo_prob_par_cur;
                    for (size_t par_next = 0; par_next < 2; ++par_next) {
                      double par_next_p = (par_next == 0) ? hypo_prob_par_next : 1.0 - hypo_prob_par_next;
                      mar = prev_c0_p*prev_c1_p*prev_p*c0_p*c1_p*next_p*par_cur_p*par_next_p;
                      double p_cur_u = combined_trans_mats[node_id][prev][par_cur][0]*
                        combined_trans_mats[children[0]][prev_c0][0][c0]*
                        combined_trans_mats[children[1]][prev_c1][0][c1]*
                        combined_trans_mats[node_id][0][par_next][next];
                      double p_cur_m =combined_trans_mats[node_id][prev][par_cur][1]*
                        combined_trans_mats[children[0]][prev_c0][1][c0]*
                        combined_trans_mats[children[1]][prev_c1][1][c1]*
                        combined_trans_mats[node_id][1][par_next][next];
                      p_cur += p_cur_u/(p_cur_u + p_cur_m)*mar;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    diff += abs(meth_prob_table[pos][node_id] - p_cur);
    meth_prob_table[pos][node_id] = p_cur;
  }
  updated[node_id] = true;
}



static void
posterior_end(const bool LEAF,
              const vector<size_t> &subtree_sizes,
              const vector<vector<double> >&G,
              const vector<vector<vector<double> > > &time_trans_mats,
              const vector<vector<vector<vector<double> > > > &combined_trans_mats,
              const size_t pos,
              const size_t node_id,
              const size_t parent_id,
              vector<bool> &updated,
              vector<vector<double> > &meth_prob_table,
              double &diff) {

  double mar;
  double p_cur = 0.0;

  if (subtree_sizes[node_id] == 1) { //leaf
    if (LEAF) {
      double p_parent_cur = meth_prob_table[pos][parent_id];
      double p_prev = meth_prob_table[pos-1][node_id];
      for (size_t prev = 0; prev < 2; ++prev) {
        double prev_p = (prev == 0) ? p_prev : 1.0 - p_prev;
        for (size_t par_cur = 0; par_cur < 2; ++ par_cur) {
          double par_cur_p = (par_cur == 0) ? p_parent_cur : 1.0- p_parent_cur;
          mar = prev_p*par_cur_p;
          double p_cur_u = combined_trans_mats[node_id][prev][par_cur][0];
          double p_cur_m = combined_trans_mats[node_id][prev][par_cur][1];
          double p_cur_mb = p_cur_u/(p_cur_u + p_cur_m);
          p_cur += p_cur_mb*mar;
        }
      }
      diff += abs(meth_prob_table[pos][node_id] - p_cur);
      meth_prob_table[pos][node_id] = p_cur;
    }
  } else {
    //check if all children are updated
    vector<size_t> children;
    size_t count = 1;
    while (count < subtree_sizes[node_id]) {
      size_t child_id = node_id + count;
      children.push_back(child_id);
      if (!updated[child_id]) {
        posterior_end(LEAF, subtree_sizes, G, time_trans_mats, combined_trans_mats,
                      pos, child_id, node_id, updated, meth_prob_table, diff);
      }
      count += subtree_sizes[child_id];
    }

    if (node_id == 0) { //root node  //root has 3 children!
      double hypo_prob_prev = meth_prob_table[pos-1][node_id];
      double hypo_prob_child_0 = meth_prob_table[pos][children[0]];
      double hypo_prob_child_1 = meth_prob_table[pos][children[1]];
      double hypo_prob_child_2 = meth_prob_table[pos][children[2]];
      double hypo_prob_prev_child_0 = meth_prob_table[pos-1][children[0]];
      double hypo_prob_prev_child_1 = meth_prob_table[pos-1][children[1]];
      double hypo_prob_prev_child_2 = meth_prob_table[pos-1][children[2]];

      for (size_t prev = 0; prev < 2; ++ prev) {
        double prev_p = (prev == 0) ? hypo_prob_prev : 1.0 - hypo_prob_prev;
        for (size_t c0 = 0; c0 < 2; ++c0) {
          double c0_p = (c0 == 0) ? hypo_prob_child_0 : 1.0 - hypo_prob_child_0;
          for (size_t c1 = 0; c1 < 2; ++c1) {
            double c1_p = (c1 == 0) ? hypo_prob_child_1 : 1.0 - hypo_prob_child_1;
            for (size_t c2 = 0; c2 < 2; ++c2) {
              double c2_p = (c2 == 0) ? hypo_prob_child_2 : 1.0 - hypo_prob_child_2;
              for (size_t prev_c0 = 0; prev_c0 < 2; ++prev_c0) {
                double prev_c0_p = (prev_c0 == 0) ? hypo_prob_prev_child_0 : 1.0 - hypo_prob_prev_child_0;
                for (size_t prev_c1 = 0; prev_c1 < 2; ++prev_c1) {
                  double prev_c1_p = (prev_c1 == 0) ? hypo_prob_prev_child_1 : 1.0 - hypo_prob_prev_child_1;
                  for (size_t prev_c2 = 0; prev_c2 < 2; ++prev_c2) {
                    double prev_c2_p = (prev_c2 == 0) ? hypo_prob_prev_child_2 : 1.0 - hypo_prob_prev_child_2;
                    mar = prev_p*c0_p*c1_p*c2_p*prev_c0_p*prev_c1_p*prev_c2_p;
                    double p_cur_u = combined_trans_mats[children[0]][prev_c0][0][c0]*
                      combined_trans_mats[children[1]][prev_c1][0][c1]*
                      combined_trans_mats[children[2]][prev_c2][0][c2]*G[prev][0];
                    double p_cur_m = combined_trans_mats[children[0]][prev_c0][1][c0]*
                      combined_trans_mats[children[1]][prev_c1][1][c1]*
                      combined_trans_mats[children[2]][prev_c2][1][c2]*G[prev][1];
                    p_cur += p_cur_u/(p_cur_m + p_cur_u)*mar;
                  }
                }
              }
            }
          }
        }
      }
      diff += abs(meth_prob_table[pos][node_id] - p_cur);
      meth_prob_table[pos][node_id] = p_cur;
    } else { // internal node // two children
      double hypo_prob_child_0 = meth_prob_table[pos][children[0]];
      double hypo_prob_child_1 = meth_prob_table[pos][children[1]];
      double hypo_prob_prev_child_0 = meth_prob_table[pos-1][children[0]];
      double hypo_prob_prev_child_1 = meth_prob_table[pos-1][children[1]];
      double hypo_prob_prev = meth_prob_table[pos-1][node_id];
      double hypo_prob_par_cur = meth_prob_table[pos][parent_id];

      for (size_t c0 = 0; c0 < 2; ++c0) {
        double c0_p = (c0 == 0) ? hypo_prob_child_0 : 1.0 - hypo_prob_child_0;
        for (size_t c1 = 0; c1 < 2; ++c1) {
          double c1_p = (c1 == 0) ? hypo_prob_child_1 : 1.0 - hypo_prob_child_1;
          for (size_t prev = 0; prev < 2; ++prev) {
            double prev_p = (prev == 0) ? hypo_prob_prev : 1.0 - hypo_prob_prev;
            for (size_t par_cur = 0; par_cur < 2; ++par_cur) {
              double par_cur_p =
                (par_cur == 0) ? hypo_prob_par_cur : 1.0 - hypo_prob_par_cur;
              for (size_t prev_c0 = 0; prev_c0 < 2; ++prev_c0) {
                double prev_c0_p =
                  (prev_c0 == 0) ? hypo_prob_prev_child_0 : 1.0 - hypo_prob_prev_child_0;
                for (size_t prev_c1 = 0; prev_c1 < 2; ++prev_c1) {
                  double  prev_c1_p = (prev_c1 == 0) ? hypo_prob_prev_child_1 : 1.0 - hypo_prob_prev_child_1;
                  mar = c0_p*c1_p*prev_p*par_cur_p*prev_c0_p*prev_c1_p;
                  double p_cur_u = combined_trans_mats[node_id][prev][par_cur][0]*
                    combined_trans_mats[children[0]][prev_c0][0][c0]*
                    combined_trans_mats[children[1]][prev_c1][0][c1];
                  double p_cur_m = combined_trans_mats[node_id][prev][par_cur][1]*
                    combined_trans_mats[children[0]][prev_c0][1][c0]*
                    combined_trans_mats[children[1]][prev_c1][1][c1];
                  p_cur += p_cur_u/(p_cur_u + p_cur_m)*mar;
                }
              }
            }
          }
        }
      }
    }
    diff += abs(meth_prob_table[pos][node_id] - p_cur);
    meth_prob_table[pos][node_id] = p_cur;
  }
  updated[node_id] = true;
}


static void
iterate_update(const vector<size_t> &subtree_sizes,
               const vector<size_t> &reset_points,
               const vector<vector<double> >&G,
               const vector<vector<vector<double> > > &time_trans_mats,
               const vector<vector<vector<vector<double> > > > &combined_trans_mats,
               const double tolerance,
               const size_t MAXITER,
               vector<vector<double> > &tree_hypo_prob_table) {

  size_t n_nodes = subtree_sizes.size();
  bool LEAF = false;
  for (size_t i = 0; i < reset_points.size()-1; ++i) {
    double diff = std::numeric_limits<double>::max();
    size_t start = reset_points[i];
    size_t end = reset_points[i+1];
    double tol = (LEAF) ? tolerance*n_nodes*(end-start) : tolerance*n_nodes*(end-start)/2;
    size_t iter = 0;
    while (diff > tol && iter < MAXITER) {
      diff = 0.0;
      for (size_t j = start; j < end; ++j) {
        vector<bool> updated (subtree_sizes.size(), false);
        if (j == start) {
          posterior_start(LEAF, subtree_sizes, G, time_trans_mats,
                          combined_trans_mats, j, 0,
                          0, updated, tree_hypo_prob_table, diff);
        } else if (j == end -1) {
          posterior_end(LEAF, subtree_sizes, G, time_trans_mats,
                        combined_trans_mats, j, 0,
                        0, updated, tree_hypo_prob_table, diff);
        } else {
          posterior_middle(LEAF, subtree_sizes, G, time_trans_mats,
                           combined_trans_mats, j, 0,
                           0, updated, tree_hypo_prob_table, diff);
        }
      }
      ++iter;
    }
  }
}


static void
approx_posterior(const vector<size_t> &subtree_sizes,
                 const vector<double> &params,
                 const vector<size_t> &reset_points,
                 const double tolerance,
                 const size_t MAXITER,
                 vector<vector<double> > &meth_prob_table) {

  const size_t n_nodes = subtree_sizes.size();
  vector<vector<vector<double> > > time_trans_mats;
  // collect probabilities and derivatives by branch
  vector<vector<vector<vector<double> > > > combined_trans_mats;
  double rate0 = params[1];
  double g0 = params[2];
  double g1 = params[3];
  vector<double> Ts(params.begin()+4, params.end());
  collect_transition_matrices(rate0, g0, g1, Ts, time_trans_mats,
                              combined_trans_mats);
  vector<vector<double> >G(2, vector<double>(2, 0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  //update iteratively
  iterate_update(subtree_sizes, reset_points, G, time_trans_mats,
                 combined_trans_mats, tolerance, MAXITER, meth_prob_table);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////         OPTIMIZATION          /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////     SLIM OPTIMIZATION         /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
get_posterior(const vector<vector<double> > &tree_prob_table,
              const size_t pos, const size_t node_id,
              const size_t state,
              double &prob) {
  prob = (state == 0)? tree_prob_table[pos][node_id] : 1.0 - tree_prob_table[pos][node_id];
}

static void
acc_triad_weight(const vector<size_t> &subtree_sizes,
                 const vector<vector<double> > &tree_prob_table,
                 const size_t pos,
                 const size_t node_id,
                 vector<vector<vector<vector<double> > > > &triad_weights) { // treesizex2x2x2
  assert(pos > 0);

  if (subtree_sizes[node_id] >1)  {
    size_t count =1;
    while (count < subtree_sizes[node_id]) {
      size_t child_id = node_id + count;
      // triad_weights
      for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
          for (size_t k = 0; k < 2; ++k) {
            double p_prev, p_cur, p_anc;
            get_posterior(tree_prob_table, pos-1, child_id, i, p_prev);
            get_posterior(tree_prob_table, pos, node_id, j, p_anc);
            get_posterior(tree_prob_table, pos, child_id, k, p_cur);
            triad_weights[child_id][i][j][k] += p_prev*p_cur*p_anc;
          }
        }
      }
      count += subtree_sizes[child_id];
    }
  }
}

static void
acc_start_weight(const vector<size_t> &subtree_sizes,
                 const vector<vector<double> > &tree_prob_table,
                 const size_t pos,
                 const size_t node_id,
                 vector<vector<vector<double> > > &start_weights) {// treesizex2x2
  assert(subtree_sizes[node_id] >1 );
  size_t count = 1;
  while (count < subtree_sizes[node_id]) {
    size_t child_id = node_id + count;
    // triad_weights
    double p_cur, p_anc;
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        get_posterior(tree_prob_table, pos, node_id, j, p_anc);
        get_posterior(tree_prob_table, pos, child_id, k, p_cur);
        start_weights[child_id][j][k] += p_cur*p_anc;
      }
    }
    count += subtree_sizes[child_id];
  }
}

static void
acc_root_weights(const vector<vector<double> > &meth_prob_table,
                 const size_t pos,
                 vector<vector<double> > &root_weights){ //2x2
  assert (pos > 0);
  // accumulate transition from position pos-1 to pos
  const size_t node_id = 0;
  double p_prev, p_cur;
  for (size_t i = 0; i < 2; ++i) {
    for (size_t k = 0; k < 2; ++k) {
      get_posterior(meth_prob_table, pos-1, node_id, i, p_prev);
      get_posterior(meth_prob_table, pos, node_id, k, p_cur);
      root_weights[i][k] += p_prev*p_cur;
    }
  }
}

static void
posterior_to_weights(const vector<size_t> &subtree_sizes,
                     const vector<double> &params,
                     const vector<size_t> &reset_points,
                     const vector<vector<double> > &tree_prob_table,
                     vector<vector<vector<vector<double> > > > &triad_weights, // treesizex2x2x2
                     vector<vector<vector<double> > > &start_weights,// treesizex2x2
                     vector<vector<double> > &root_weights){ //2x2

  const size_t n_nodes = subtree_sizes.size();
  vector<vector<double> > mat2x2(2,vector<double>(2,0.0));
  vector<vector<vector<double> > > mat2x2x2(2, mat2x2);
  root_weights = mat2x2;
  start_weights = vector<vector<vector<double> > >(n_nodes, mat2x2);
  triad_weights = vector<vector<vector<vector<double> > > >(n_nodes, mat2x2x2);

  for  (size_t block = 0; block < reset_points.size() -1; ++block) {
    // accumulate start_weights
    for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
      if (subtree_sizes[node_id]>1) {
        acc_start_weight(subtree_sizes, tree_prob_table,
                         reset_points[block], node_id, start_weights);
      }
    }
    for (size_t pos = reset_points[block]+1; pos < reset_points[block+1]; ++pos ) {
      //accumulate root_weights
      acc_root_weights(tree_prob_table, pos, root_weights);
      // accumulate triad_weight
      for (size_t node_id = 0; node_id < n_nodes; ++node_id) {
        if (subtree_sizes[node_id]>1)
          acc_triad_weight(subtree_sizes, tree_prob_table, pos,
                           node_id, triad_weights);
      }
    }
  }
}

static void
objective_branch(const vector<double> &params,
                 const vector<vector<vector<vector<double> > > > &triad_weights, // treesizex2x2x2
                 const vector<vector<vector<double> > > &start_weights,
                 const vector<vector<vector<double> > > &time_trans_mats,
                 const vector<vector<vector<vector<double> > > > &combined_trans_mats,// treesizex2x2
                 const vector<vector<vector<vector<double> > > > &combined_trans_mats_dT,
                 const size_t node_id,
                 double &F, double &deriv) {

  F = 0.0;
  deriv = 0.0;
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        F += triad_weights[node_id][i][j][k]*log(combined_trans_mats[node_id][i][j][k]);
        deriv += triad_weights[node_id][i][j][k]*(combined_trans_mats_dT[node_id][i][j][k]/combined_trans_mats[node_id][i][j][k]);
      }
    }
  }

  double rate0 = params[1];
  vector<vector<double> > time_trans_mats_dT(2, vector<double>(2, 0.0));
  time_trans_mats_dT[0][0] = -rate0;
  time_trans_mats_dT[0][1] = rate0;
  time_trans_mats_dT[1][0] = 1.0-rate0;
  time_trans_mats_dT[1][1] = rate0-1.0;

  for (size_t j = 0; j < 2; ++j) {
    for (size_t k = 0; k < 2; ++k) {
      F += start_weights[node_id][j][k]*log(time_trans_mats[node_id][j][k]);
      deriv += start_weights[node_id][j][k]*(time_trans_mats_dT[j][k]/time_trans_mats[node_id][j][k]);
    }
  }
}

static void
update_branch(const double TOL,
              const vector<size_t> &subtree_sizes,
              const vector<double> &params, const size_t node_id,
              const vector<vector<vector<vector<double> > > > &triad_weights,
              const vector<vector<vector<double> > > &start_weights,
              double &T){

  double rate0 = params[1];
  double g0 = params[2];
  double g1 = params[3];
  vector<double> Ts(params.begin()+4, params.end());
  vector<vector<vector<double> > > time_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats_drate;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg0;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg1;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dT;
  collect_transition_matrices_deriv(rate0, g0, g1, Ts, time_trans_mats, combined_trans_mats,
                                    combined_trans_mats_drate, combined_trans_mats_dg0,
                                    combined_trans_mats_dg1, combined_trans_mats_dT);
  vector<double> prev_params = params;
  vector<double> new_params = params;

  double prev_F, new_F;
  double prev_deriv, new_deriv;
  objective_branch(prev_params, triad_weights, start_weights, time_trans_mats,
                   combined_trans_mats, combined_trans_mats_dT,
                   node_id, prev_F, prev_deriv);
  cerr << prev_params[4 + node_id] << "\t" << prev_deriv<< "\t" << prev_F << "\t" << 1.0 << endl;

  double frac = 1.0;
  while (frac > TOL && (frac*sign(prev_deriv) + prev_params[4+node_id] <  TOL ||
                        frac*sign(prev_deriv) + prev_params[4+node_id] > 1.0-TOL ) ) {
    frac = frac/2;
  }

  bool CONVERGE = false;
  bool SUCCESS = false;
  double improve = 0.0;
  while (!CONVERGE || !SUCCESS) {
    new_params[4+node_id] = prev_params[4+node_id] + frac*sign(prev_deriv);
    // update trans mats and derivs with new-params
    vector<double> trial_Ts(new_params.begin()+4, new_params.end());
    collect_transition_matrices_deriv(rate0, g0, g1, trial_Ts, time_trans_mats, combined_trans_mats,
                                      combined_trans_mats_drate, combined_trans_mats_dg0,
                                      combined_trans_mats_dg1, combined_trans_mats_dT);

    objective_branch(new_params, triad_weights, start_weights, time_trans_mats,
                     combined_trans_mats, combined_trans_mats_dT,
                     node_id, new_F, new_deriv);

    if (new_F > prev_F) {
      SUCCESS= true;
      improve += new_F - prev_F;
      cerr << "SUCCEEDED\tImprov=" << improve << endl;
      prev_F = new_F;
      prev_deriv = new_deriv;
      prev_params = new_params;
    } else {
      frac = frac/2;
    }
    if (frac < TOL) CONVERGE = true;
  }
  T = prev_params[4+node_id];
}


static void
objective_rate(const vector<size_t> &subtree_sizes,
               const vector<double> &params,
               const vector<vector<vector<vector<double> > > > &triad_weights, // treesizex2x2x2
               const vector<vector<vector<double> > > &start_weights,
               const vector<vector<vector<double> > > &time_trans_mats,
               const vector<vector<vector<vector<double> > > > &combined_trans_mats, // treesizex2x2
               const vector<vector<vector<vector<double> > > > &combined_trans_mats_drate,
               double &F, double &deriv_rate) {

  F = 0.0;
  deriv_rate = 0.0;

  const size_t n_nodes = subtree_sizes.size();
  const vector<double> T(params.begin()+4, params.end());

  for (size_t node_id = 1; node_id < n_nodes; ++ node_id) {
    vector<vector<double> > time_trans_mats_drate(2, vector<double>(2,0.0));
    time_trans_mats_drate[0][0] = -T[node_id];
    time_trans_mats_drate[0][1] = T[node_id];
    time_trans_mats_drate[1][0] = -T[node_id];
    time_trans_mats_drate[1][1] = T[node_id];

    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        for (size_t k = 0; k < 2; ++k) {
          F += triad_weights[node_id][i][j][k]*log(combined_trans_mats[node_id][i][j][k]);
          deriv_rate += triad_weights[node_id][i][j][k]*(combined_trans_mats_drate[node_id][i][j][k]/combined_trans_mats[node_id][i][j][k]);
        }
      }
    }

    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        F += start_weights[node_id][j][k]*log(time_trans_mats[node_id][j][k]);
        deriv_rate += start_weights[node_id][j][k]*time_trans_mats_drate[j][k]/time_trans_mats[node_id][j][k];
      }
    }
  }
}

static void
update_rate(const double TOL,
            const vector<size_t> &subtree_sizes,
            const vector<double> &params,
            const vector<vector<vector<vector<double> > > > &triad_weights, // treesizex2x2x2
            const vector<vector<vector<double> > > &start_weights,
            double &new_rate){

  const double rate0 = params[1];
  const double g0 = params[2];
  const double g1 = params[3];
  const vector<double> Ts(params.begin()+4, params.end());
  vector<vector<vector<double> > > time_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats_drate;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg0;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg1;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dT;
  collect_transition_matrices_deriv(rate0, g0, g1, Ts, time_trans_mats, combined_trans_mats,
                                    combined_trans_mats_drate, combined_trans_mats_dg0,
                                    combined_trans_mats_dg1, combined_trans_mats_dT);
  vector<double> prev_params = params;
  vector<double> new_params = params;

  double prev_F = 0.0;
  double new_F = 0.0;
  double prev_deriv, new_deriv;

  objective_rate(subtree_sizes, prev_params, triad_weights,
                 start_weights, time_trans_mats,
                 combined_trans_mats, combined_trans_mats_drate,
                 prev_F, prev_deriv);

  double denom = abs(prev_deriv);
  double frac = 1.0;
  while (frac > TOL && (frac*sign(prev_deriv) + prev_params[1] > 1 - TOL ||
                        frac*sign(prev_deriv) + prev_params[1] < TOL)) {
    frac = frac/2;
  }

  bool CONVERGE = false;
  bool SUCCESS = false;
  double improve = 0.0;
  while (!CONVERGE || !SUCCESS) {
    new_params[1] = prev_params[1] + frac*(prev_deriv/denom); //rate0
    collect_transition_matrices_deriv(new_params[1], g0, g1, Ts,
                                      time_trans_mats, combined_trans_mats,
                                      combined_trans_mats_drate, combined_trans_mats_dg0,
                                      combined_trans_mats_dg1, combined_trans_mats_dT);
    objective_rate(subtree_sizes, new_params, triad_weights,
                   start_weights, time_trans_mats, combined_trans_mats,
                   combined_trans_mats_drate, new_F, new_deriv);

    if (new_F > prev_F) {
      SUCCESS= true;
      improve += new_F - prev_F;
      prev_F = new_F;
      prev_deriv = new_deriv;
      prev_params = new_params;
      denom = abs(new_deriv);
      cerr << "SUCCESS\tImprov=" << improve << "\trate=" << new_params[1] <<endl;
    } else {
      frac = frac/2;
    }
    if (frac < TOL ) CONVERGE = true;
  }
  new_rate = prev_params[1];
}

static void
objective_G(const vector<size_t> &subtree_sizes,
            const vector<double> &params,
            const vector<vector<vector<vector<double> > > > &triad_weights, // treesizex2x2x2
            const vector<vector<double> > &root_weights,
            const vector<vector<vector<vector<double> > > > &combined_trans_mats,
            const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg0,
            const vector<vector<vector<vector<double> > > > &combined_trans_mats_dg1,
            double &F, vector<double> &deriv_G){

  double g0 = params[2];
  double g1 = params[3];

  vector<vector<double> > G(2,vector<double>(2,0.0));
  G[0][0] = g0;
  G[0][1] = 1.0 - g0;
  G[1][0] = 1.0 - g1;
  G[1][1] = g1;

  F = 0.0;
  deriv_G = vector<double>(2, 0.0);
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++ node_id) {
    for (size_t i = 0; i < 2; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        for (size_t k = 0; k < 2; ++k) {
          F += triad_weights[node_id][i][j][k]*log(combined_trans_mats[node_id][i][j][k]);
        }
      }
    }
    for (size_t j = 0; j < 2; ++j) {
      for (size_t k = 0; k < 2; ++k) {
        deriv_G[0] += triad_weights[node_id][0][j][k]*(combined_trans_mats_dg0[node_id][0][j][k]/combined_trans_mats[node_id][0][j][k]);
        deriv_G[1] += triad_weights[node_id][1][j][k]*(combined_trans_mats_dg1[node_id][1][j][k]/combined_trans_mats[node_id][1][j][k]);
      }
    }
  }

  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      F += root_weights[i][j]*log(G[i][j]);
    }
  }

  deriv_G[0] +=  root_weights[0][0]/G[0][0] - root_weights[0][1]/G[0][1];
  deriv_G[1] +=  -1.0*root_weights[1][0]/G[1][0] + root_weights[1][1]/G[1][1];
}

static void
update_G(const double TOL,
         const vector<size_t> &subtree_sizes,
         const vector<double> &params,
         const vector<vector<vector<vector<double> > > > &triad_weights, // treesizex2x2x2
         const vector<vector<double> > &root_weights,
         double &new_g0, double &new_g1) {

  const double rate0 = params[1];
  const double g0 = params[2];
  const double g1 = params[3];
  const vector<double> Ts(params.begin()+4, params.end());
  vector<vector<vector<double> > > time_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats;
  vector<vector<vector<vector<double> > > > combined_trans_mats_drate;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg0;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dg1;
  vector<vector<vector<vector<double> > > > combined_trans_mats_dT;
  collect_transition_matrices_deriv(rate0, g0, g1, Ts, time_trans_mats, combined_trans_mats,
                                    combined_trans_mats_drate, combined_trans_mats_dg0,
                                    combined_trans_mats_dg1, combined_trans_mats_dT);
  vector<double> prev_params = params;
  vector<double> new_params = params;

  double prev_F =0.0;
  double new_F = 0.0;
  vector<double> prev_deriv, new_deriv;

  objective_G(subtree_sizes, prev_params, triad_weights,
              root_weights, combined_trans_mats,combined_trans_mats_dg0,
              combined_trans_mats_dg1, prev_F, prev_deriv);

  double denom = abs(prev_deriv[0]) + abs(prev_deriv[1]);
  double frac = 1.0;
  while (frac > TOL && (frac*(prev_deriv[0]/denom) + prev_params[2] > 1- TOL ||
                        frac*(prev_deriv[0]/denom) + prev_params[2] < TOL ||
                        frac*(prev_deriv[1]/denom) + prev_params[3] > 1- TOL ||
                        frac*(prev_deriv[1]/denom) + prev_params[3] < TOL ) ) {
    frac = frac/2;
  }

  bool CONVERGE = false;
  bool SUCCESS = false;
  double improve = 0.0;
  while (!CONVERGE || !SUCCESS) {
    new_params[2] = prev_params[2] + frac*(prev_deriv[0]/denom); //g0
    new_params[3] = prev_params[3] + frac*(prev_deriv[1]/denom); //g1
    collect_transition_matrices_deriv(rate0, new_params[2], new_params[3], Ts,
                                      time_trans_mats, combined_trans_mats,
                                      combined_trans_mats_drate, combined_trans_mats_dg0,
                                      combined_trans_mats_dg1, combined_trans_mats_dT);
    objective_G(subtree_sizes, new_params, triad_weights, root_weights,
                combined_trans_mats, combined_trans_mats_dg0,
                combined_trans_mats_dg1, new_F, new_deriv);
    if (new_F > prev_F) {
      SUCCESS= true;
      improve += new_F - prev_F;
      prev_F = new_F;
      prev_deriv = new_deriv;
      prev_params = new_params;
      denom = abs(new_deriv[0]) + abs(new_deriv[1]);

      cerr << "SUCCESS\tImprov=" << improve << "\tG=("
           << new_params[2] << ", " << new_params[3] << ")" << endl;
    } else {
      frac = frac/2;
    }
    if (frac < TOL) CONVERGE = true;
  }
  new_g0 = prev_params[2];
  new_g1 = prev_params[3];
}


static void
update_pi0(const vector<vector<double> > &tree_prob_table,
           const vector<size_t> &reset_points,
           double &pi0){
  pi0 = 0.0;
  for (size_t i = 0; i < reset_points.size()-1; ++i) {
    pi0 += tree_prob_table[reset_points[i]][0];
  }
  pi0 = pi0/(reset_points.size()-1);
}


static void
optimize_params(const vector<size_t> &subtree_sizes,
                const vector<vector<double> > &tree_prob_table,
                const vector<size_t> &reset_points,
                const vector<double> &params,
                vector<double> &newparams){

  const size_t n_nodes = subtree_sizes.size();
  vector<vector<double> > mat2x2(2,vector<double>(2,0.0));
  vector<vector<vector<double> > > mat2x2x2(2, mat2x2);
  vector<vector<double> > root_weights = mat2x2;
  vector<vector<vector<double> > > start_weights(n_nodes, mat2x2);
  vector<vector<vector<vector<double> > > > triad_weights(n_nodes, mat2x2x2);
  posterior_to_weights(subtree_sizes, params, reset_points,
                       tree_prob_table, triad_weights,
                       start_weights, root_weights);
  newparams = params;

  double TOL = 1e-4;
  double new_g0;
  double new_g1;
  update_G(TOL, subtree_sizes, params,
           triad_weights, root_weights, new_g0, new_g1);
  newparams[2] = new_g0;
  newparams[3] = new_g1;

  // update branches
  double T = 0.0;
  for (size_t node_id = 1; node_id < subtree_sizes.size(); ++ node_id) {
    update_branch(TOL, subtree_sizes, newparams, node_id,
                  triad_weights, start_weights, T);
    newparams[4+node_id] = T;
  }
  // update rate0
  double new_rate;
  update_rate(TOL, subtree_sizes, newparams, triad_weights, start_weights,
              new_rate);
  newparams[1] = new_rate;

  //update pi0
  double new_pi0;
  update_pi0(tree_prob_table, reset_points, new_pi0);
  newparams[0] = new_pi0;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static void
tree_prob_to_states(const vector<vector<double> > &tree_prob_table,
                    const double cutoff,
                    vector<string> &states) {
  const size_t n_nodes = tree_prob_table[0].size();
  for (size_t i = 0; i < tree_prob_table.size(); ++i) {
    string state;
    for (size_t j = 0; j < n_nodes; ++j ) {
      state += ((tree_prob_table[i][j] <= cutoff) ? '1' : '0' );
    }
    states.push_back(state);
  }
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
      state += ((tree_prob_table[i][j] <= 0.5) ? '1' : '0' );
      if (subtree_sizes[j] == 1)
        pme += ((tree_prob_table[i][j] <= 0.5) ? '1' : '0' );
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
build_domain(const size_t minCpG,
             const size_t desert_size,
             const vector<Site> &sites,
             const vector<string> &states,
             vector<GenomicRegion> &domains) {
  // first round collapsing
  GenomicRegion cpg;
  cpg.set_chrom(sites[0].chrom);
  cpg.set_start(sites[0].pos);
  cpg.set_end(sites[0].pos+1);
  cpg.set_name(states[0]);
  cpg.set_score(1.0);
  for (size_t i = 1; i < sites.size(); ++i) {
    size_t d = distance_between_sites(sites[i], sites[i-1]);
    if (d < desert_size && states[i] == states[i-1]) {
      cpg.set_score(1.0+ cpg.get_score());
      cpg.set_end(sites[i].pos);
    } else {
      domains.push_back(cpg);
      cpg.set_chrom(sites[i].chrom);
      cpg.set_start(sites[i].pos);
      cpg.set_end(sites[i].pos+1);
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
            domains[j].distance(merging_domains.back()) < desert_size ) {
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
    size_t minfragcpg = 5;
    double tolerance = 1e-4;
    size_t MAXITER = 10;
    string outfile;
    string paramfile;

    // run mode flags
    bool VERBOSE = false;
    bool SINGLE = false;

    OptionParser opt_parse(strip_path(argv[0]), "test funcitonality of "
                           "phylo-methylome segmentation",
                           "<newick> <hypoprob-tab>");
    opt_parse.add_opt("minCpG", 'm', "minimum observed #CpGs in a block"
                      "(default: 50)", false, minCpG);
    opt_parse.add_opt("maxiter", 'i', "maximum iteration"
                      "(default: 10)", false, MAXITER);
    opt_parse.add_opt("verbose", 'v', "print more run info (default: false)",
                      false, VERBOSE);
    opt_parse.add_opt("params", 'p', "given parameters", false, paramfile);
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("minfragCpG", 'f', "ignore fragments with fewer CpG sites"
                      "(default: 5)", false, minfragcpg);
    opt_parse.add_opt("single", 's', "also output states by sites (when -o is used)",
                      false, SINGLE);

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
    /******************** INITIALIZE PARAMETERS *******************************/
    double pi0 = 0.16181;
    double rate0 = 0.767393;
    double g0 = 0.99013;
    double g1 = 0.983651;

    bool PARAMFIX = false;
    if (!paramfile.empty()) {
      read_params(paramfile, pi0, rate0, g0, g1, t);
      if (VERBOSE) {
        cerr << "Read in tree branch lengths:" << t.tostring() << endl
             << "Root_unmeth_prob=" << pi0 << "\trate=" <<  rate0
             << "\tg0=" << g0 << "\tg1" << g1 << endl;
      }
      PARAMFIX = true;
      branches.clear();
      t.get_branch_lengths(branches);
    }
    vector<double> Ts;
    for (size_t i = 0; i < branches.size(); ++i) {
      Ts.push_back(1.0 - 1.0/exp(branches[i]));
    }

    const size_t n_nodes = subtree_sizes.size();
    vector<double> start_param(n_nodes+4, 0.0);
    start_param[0] = pi0; start_param[1] = rate0;
    start_param[2] = g0; start_param[3] = g1;
    for (size_t i = 1; i < n_nodes; ++i) {
      start_param[4+i] = Ts[i];
    }

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

    if (VERBOSE) {
      vector<size_t> species_in_order;
      subtree_sizes_to_leaves_preorder(subtree_sizes, species_in_order);
      for (size_t i = 0; i < species_in_order.size(); ++i) {
        cerr << meth_table_species[i] << "\t" << species_in_order[i] << endl;
      }
    }

    /**************************************************************************/
    /*************************  SEPARATE BY DESERTS ***************************/
    if (VERBOSE)
      cerr << "Read in total " << sites.size() << " sites" << endl
           << "Separating sites" << endl;
    vector<size_t> reset_points;
    double g0_est, g1_est, pi0_est;
    separate_regions(VERBOSE, desert_size, minCpG, meth_prob_table,
                     sites, reset_points, g0_est, g1_est, pi0_est);

    if (!PARAMFIX) {
      start_param[0] = pi0_est;
      start_param[2] = g0_est;
      start_param[3] = g1_est;
    }

    if (VERBOSE) {
      cerr << "After filtering: " << sites.size()
           << "in " << reset_points.size() -1 << "blocks" << endl;

      if (PARAMFIX) cerr << "Given parameter:\t";
      else cerr << "Start parameter:\t";
      for (size_t i = 0; i < start_param.size(); ++i) {
        if (i < 4) cerr << start_param[i] << "\t";
        if (i > 4) cerr << "[branch "<< i-4 <<"]=" << -log(1.0 - start_param[i]) << "\t";
      }
      cerr << endl;
    }

    vector<vector<double> > tree_prob_table(sites.size(), vector<double>(n_nodes, -1.0));
    leaf_to_tree_prob(subtree_sizes, meth_prob_table, tree_prob_table);
    if (VERBOSE)
      cerr << "Initialized internal probability" << endl;

    // make a copy of the initialized tree_prob_table
    vector<vector<double> > tree_prob_table_copy = tree_prob_table;

    /**************************************************************************/
    /******************** APPROXIMATE OPTIMIZATION ****************************/
    if (VERBOSE) cerr << "Mode: Approx posterior" << endl;
    const double tol = 1e-3;
    const size_t max_app_iter = 100;

    double diff = std::numeric_limits<double>::max();
    size_t iter = 0;
    while (iter < MAXITER && diff > tol) {
      cerr << "Iteration " << iter << "\t";
      for (size_t i = 0; i < start_param.size(); ++i)
        cerr << start_param[i] << "\t";
      cerr << endl;
      vector<double> newparams(n_nodes+4, 0.0);

      // start from the initial table in each iteration.
      //tree_prob_table = tree_prob_table_copy;
      approx_posterior(subtree_sizes, start_param, reset_points, tolerance,
                       max_app_iter, tree_prob_table);

      optimize_params(subtree_sizes, tree_prob_table,
                      reset_points, start_param, newparams);

      diff = 0.0;
      for (size_t i = 0; i < newparams.size(); ++i) {
        diff += abs(newparams[i]-start_param[i]);
      }
      start_param = newparams;
      ++iter;
    }

    if (diff <= tol && VERBOSE)
      cerr << "Converged at iteration " << iter << endl;

    for (size_t i = 5; i < start_param.size(); ++i) {
      branches[i-4] = -log(1.0 - start_param[i]);
    }
    t.set_branch_lengths(branches);

    if (VERBOSE) {
      cerr << endl << "[Results]\t";
      for (size_t i = 0; i < 4; ++i)
        cerr << start_param[i] << "\t";
      for (size_t i = 5; i < start_param.size(); ++i) {
        cerr << branches[i-4] << "\t";
      }
      cerr << endl;
      cerr << "Final pass of posterior approximation" << endl;
    }

    tree_prob_table = tree_prob_table_copy;
    approx_posterior(subtree_sizes, start_param, reset_points, tolerance,
                     max_app_iter, tree_prob_table);

    if (!outfile.empty()) {
      std::ofstream out(outfile.c_str());
      if (!out)
        throw SMITHLABException("bad output file: " + outfile);
      vector<string> states;
      vector<GenomicRegion> domains;
      double cutoff = 0.5;
      tree_prob_to_states(tree_prob_table, cutoff, states);
      build_domain(minfragcpg, desert_size, sites, states, domains);
      cerr << domains.size() << endl;

      out << "#" << t.Newick_format() << endl;
      for (size_t i = 0; i < domains.size(); ++i) {
        out << domains[i] << '\n';
      }

      if (SINGLE) {
        string outssfile = outfile + "_bysite";
        std::ofstream outss(outssfile.c_str());
        if (!outss)
          throw SMITHLABException("bad output file: " + outssfile);
        write_treeprob_states(subtree_sizes, sites, tree_prob_table, outss);
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
