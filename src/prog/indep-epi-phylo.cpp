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
#include "BetaBin.hpp"

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
    if (meth.back()[i] < 0 || meth.back()[i] > 1.0)
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

// Fix f/b parameters for now
// change this part later
// pair<size_t, size_t> > : <Coverge, M>
void
get_betabin_params(const vector<vector<pair<size_t, size_t> > > &meth_table,
                   const size_t col, double &fa, double &fb,
                   double &ba, double &bb) {
  fa = 0.443;
  fb = 4.314;
  ba = 2.576;
  bb = 0.491;
}


//meth_table dimemsion: nsites x nsample
//pair in meth_table: <coverage, meth>
static void
meth_to_posterior(const vector<vector<pair<size_t, size_t> > > &meth_table,
                  const size_t col, vector<double> &statepost,
                  double &p) {
  double fa, fb, ba, bb;
  get_betabin_params(meth_table, col, fa, fb, ba, bb);
  betabin FG(fa, fb), BG(ba, bb);
  const size_t nsites = meth_table.size();
  for (size_t i = 0; i < nsites; ++ i) {
    if (meth_table[i][col].first >= 1) {
      const double cov = meth_table[i][col].first;
      const double mr = meth_table[i][col].second;
      const pair<double, double> MU = make_pair(mr, cov - mr);
      const double fllk = FG.log_likelihood(MU);
      const double bllk = BG.log_likelihood(MU);
      const double sum = exp(fllk) + exp(bllk);
      statepost.push_back(exp(fllk)/sum);
      p += exp(fllk)/sum;
    } else {
      statepost.push_back(0.5);
      p += 0.5;
    }
  }
}

static void
meth_to_posterior(const vector<vector<double> > &meth_table,
                  const size_t col, vector<double> &statepost,
                  double &p) {
  const size_t nsites = meth_table.size();
  for (size_t i = 0; i < nsites; ++ i) {
    statepost.push_back(meth_table[i][col]);
    p += meth_table[i][col];
  }
}


static bool
has_same_species_order(const PhyloTreePreorder &the_tree,
                       const vector<string> &meth_table_species) {
  vector<string> leaf_names;
  the_tree.get_leaf_names(leaf_names);
  return leaf_names == meth_table_species;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////  !!! ADS: CODE ABOVE IS TEMPORARY ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



struct um_val {
  double u;
  double m;

  um_val() : u(0.0), m(0.0) {}
  um_val(const double u_, const double m_) : u(u_), m(m_) {}
  double dot(const vector<double> &v) const {return u*v[0] + m*v[1];}
};

double
dot(const um_val &a, const um_val &b) {
  return a.u*b.u + a.m*b.m;
}


struct um_mat {
  um_val u;
  um_val m;
  um_mat() : u(um_val()), m(um_val()) {}
  um_mat(const um_val u_, const um_val m_) : u(u_), m(m_) {}
};



static void
trans_prob_mat(const double rate0, const double rate1,
               const double interval,
               vector<vector<double> > &transition_matrix) {
  assert(rate0 > 0 && rate1 > 0 && interval > 0);
  const double h = 1.0/exp(interval*(rate0 + rate1));

  transition_matrix = vector<vector<double> >(2, vector<double>(2,0.0));

  transition_matrix[0][0] = (rate0*h + rate1)/(rate0 + rate1);
  transition_matrix[0][1] = 1.0 - transition_matrix[0][0];

  transition_matrix[1][1] = (rate0 + rate1*h)/(rate0 + rate1);
  transition_matrix[1][0] = 1.0 - transition_matrix[1][1];
}


static void
unitrate_trans_prob_mat(const double rate0,
                        const double exp_interval,
                        vector<vector<double> > &transition_matrix) {
  assert(rate0 > 0 && rate0 <1 && exp_interval > 1);
  const double rate1 = 1.0 - rate0;
  const double h = 1.0/exp_interval;

  transition_matrix = vector<vector<double> >(2, vector<double>(2,0.0));

  transition_matrix[0][0] = (rate0*h + rate1)/(rate0 + rate1);
  transition_matrix[0][1] = 1.0 - transition_matrix[0][0];

  transition_matrix[1][1] = (rate0 + rate1*h)/(rate0 + rate1);
  transition_matrix[1][0] = 1.0 - transition_matrix[1][1];
}


static void
unitrate_trans_prob_mat(const double rate0,
                        const double exp_interval,
                        um_mat &transition_matrix) {
  assert(rate0 > 0 && rate0 <1 && exp_interval > 1);
  const double rate1 = 1.0 - rate0;
  const double h = 1.0/exp_interval;

  transition_matrix = um_mat();

  transition_matrix.u.u = (rate0*h + rate1)/(rate0 + rate1);
  transition_matrix.u.m = 1.0 - transition_matrix.u.u;
  transition_matrix.m.m = (rate0 + rate1*h)/(rate0 + rate1);
  transition_matrix.m.u = 1.0 - transition_matrix.m.m;
}


// static void
// get_rate_helper_inc(const double u, const double m,
//                     const double u_dr, const double m_dr,
//                     const double unmeth_child_lik, const double meth_child_lik,
//                     const vector<vector<double> > &child_trans_mat,
//                     const double exp_branch,
//                     double &u_inc, double &m_inc){
//   u_inc =
//     (u_dr*child_trans_mat[0][0] + u*(1.0/exp_branch - 1.0) +
//      m_dr*child_trans_mat[0][1] + m*(1.0-1.0/exp_branch))/unmeth_child_lik;
//   m_inc =
//     (u_dr*child_trans_mat[1][0] + u*(1.0/exp_branch - 1.0) +
//      m_dr*child_trans_mat[1][1] + m*(1.0-1.0/exp_branch))/meth_child_lik;
// }


// static void
// update_bysite(const vector<size_t> &subtree_sizes,
//               const size_t tree_start,
//               const size_t child_start,
//               const vector<vector<double> > &Q,
//               const vector<vector<double> > &child_trans_mat,
//               const double u, const double m,
//               const double u_dr, const double m_dr,
//               const vector<double> &u_db, const vector<double> &m_db,
//               const vector<double> &exp_branches,
//               double &subtree_lik_given_root_unmeth,
//               double &subtree_lik_given_root_meth,
//               double &rate_helper_u,
//               double &rate_helper_m,
//               vector<double> &branch_helper_u,
//               vector<double> &branch_helper_m){
//   const double unmeth_child_lik =
//     child_trans_mat[0][0]*u + child_trans_mat[0][1]*m;
//   const double meth_child_lik =
//     child_trans_mat[1][0]*u + child_trans_mat[1][1]*m;

//   //likelihood
//   subtree_lik_given_root_unmeth *= unmeth_child_lik;
//   subtree_lik_given_root_meth *= meth_child_lik;

//   //deriv_rate helper
//   double exp_branch = exp_branches[child_start];
//   double u_inc, m_inc;
//   get_rate_helper_inc(u, m, u_dr, m_dr,
//                       unmeth_child_lik, meth_child_lik,
//                       child_trans_mat, exp_branch, u_inc, m_inc);
//   rate_helper_u += u_inc;
//   rate_helper_m += m_inc;

//   //deriv_branch helper
//   branch_helper_u[child_start - tree_start] =
//     (u*Q[0][0] + m*Q[0][1])/exp_branches[child_start]/unmeth_child_lik;
//   branch_helper_m[child_start - tree_start] =
//     (u*Q[1][0] + m*Q[1][1])/exp_branches[child_start]/meth_child_lik;

//   for (size_t i = 1; i < subtree_sizes[child_start]; ++i) {
//     const size_t subtree_branch_id = i + child_start - tree_start;
//     branch_helper_u[subtree_branch_id] =
//       (u_db[i]*child_trans_mat[0][0] + m_db[i]*child_trans_mat[0][1])/unmeth_child_lik;
//     branch_helper_m[subtree_branch_id] =
//       (u_db[i]*child_trans_mat[1][0] + m_db[i]*child_trans_mat[1][1])/meth_child_lik;
//   }
// }


// static void
// subtree_likelihood_deriv(const vector<size_t> &subtree_sizes,
//                          const double &rate,
//                          const vector<double> &expbranches,
//                          const size_t tree_start,
//                          const vector<vector<double> > &states,
//                          vector<double> &subtree_lik_given_root_unmeth,
//                          vector<double> &subtree_lik_given_root_meth,
//                          vector<double> &deriv_rate_given_root_unmeth,
//                          vector<double> &deriv_rate_given_root_meth,
//                          vector<vector<double> > &deriv_branches_given_root_unmeth,
//                          vector<vector<double> > &deriv_branches_given_root_meth) {
//   const size_t nsites = states[0].size();
//   vector<double> empty_dbranch(subtree_sizes[tree_start], 0.0); //for branches in the subtree
//   deriv_branches_given_root_unmeth =
//     vector<vector<double> >(nsites, empty_dbranch);
//   deriv_branches_given_root_meth =
//     vector<vector<double> >(nsites, empty_dbranch);

//   if (subtree_sizes[tree_start] == 1) {
//     for (size_t pos = 0; pos < nsites; ++pos) {
//       assert(states[tree_start][pos] >= 0);
//       subtree_lik_given_root_unmeth.push_back(states[tree_start][pos]);
//       subtree_lik_given_root_meth.push_back(1.0 - states[tree_start][pos]);
//     }
//     deriv_rate_given_root_unmeth = vector<double>(nsites, 0.0);
//     deriv_rate_given_root_meth = vector<double>(nsites, 0.0);
//   }
//   else {
//     vector<vector<double> > Q(2, vector<double>(2, 0.0));
//     Q[0][0] = 0.0 - rate;
//     Q[0][1] = rate;
//     Q[1][0] = 1.0 - rate;
//     Q[1][1] = rate - 1.0;

//     subtree_lik_given_root_unmeth = vector<double>(nsites, 1.0);
//     subtree_lik_given_root_meth = vector<double>(nsites, 1.0);
//     deriv_rate_given_root_unmeth = vector<double>(nsites, 0.0);
//     deriv_rate_given_root_meth = vector<double>(nsites, 0.0);

//     vector<vector<double> > branch_helper_u(nsites, empty_dbranch);
//     vector<vector<double> > branch_helper_m(nsites, empty_dbranch);

//     vector<double> rate_helper_u(nsites, 0.0);
//     vector<double> rate_helper_m(nsites, 0.0);

//     size_t child_start = 0;
//     for (size_t count = 1; count < subtree_sizes[tree_start];
//          count += subtree_sizes[child_start]) {
//       child_start = tree_start + count;
//       vector<vector<double> > child_trans_mat;
//       unitrate_trans_prob_mat(rate, expbranches[child_start], child_trans_mat);
//       vector<double> u, m, u_dr, m_dr;
//       vector<vector<double> > u_db, m_db;
//       subtree_likelihood_deriv(subtree_sizes, rate, expbranches, child_start,
//                                states, u, m, u_dr, m_dr, u_db, m_db);

//       for (size_t pos = 0; pos < nsites; ++pos) {
//         update_bysite(subtree_sizes, tree_start, child_start, Q, child_trans_mat,
//                       u[pos], m[pos], u_dr[pos], m_dr[pos],
//                       u_db[pos], m_db[pos], expbranches,
//                       subtree_lik_given_root_unmeth[pos],
//                       subtree_lik_given_root_meth[pos],
//                       rate_helper_u[pos], rate_helper_m[pos],
//                       branch_helper_u[pos], branch_helper_m[pos]);
//       }
//     }

//     for (size_t pos = 0; pos < nsites; ++pos) {
//       //deriv_rate
//       deriv_rate_given_root_unmeth[pos] =
//         rate_helper_u[pos]*subtree_lik_given_root_unmeth[pos];
//       deriv_rate_given_root_meth[pos] =
//         rate_helper_m[pos]*subtree_lik_given_root_meth[pos];

//       //deriv_branches
//       for (size_t i = 1; i < subtree_sizes[tree_start]; ++i) {
//         deriv_branches_given_root_unmeth[pos][i] =
//           branch_helper_u[pos][i]*subtree_lik_given_root_unmeth[pos];
//         deriv_branches_given_root_meth[pos][i]=
//           branch_helper_m[pos][i]*subtree_lik_given_root_meth[pos];
//       }
//     }
//   }
// }



static void
subtree_likelihood_site_deriv(const vector<size_t> &subtree_sizes,
                              const double &rate, const um_mat &Q,
                              const vector<double> &expbranches,
                              const vector<um_mat> trans_mat,
                              const size_t tree_start,
                              const vector<vector<double> > &states,
                              const size_t &pos,
                              um_val &likelihood, um_val &deriv_rate,
                              vector<um_val> &deriv_branch) {

  if (subtree_sizes[tree_start] == 1) {
    assert(states[tree_start][pos] >= 0);
    likelihood = um_val(states[tree_start][pos], 1.0 - states[tree_start][pos]);
    deriv_rate = um_val(0.0, 0.0);
  }
  else {
    um_val rate_helper = um_val(0.0, 0.0);
    vector<um_val> branch_helper(subtree_sizes[tree_start]);
    likelihood = um_val(1.0, 1.0);

    size_t i = 1;
    while (i < subtree_sizes[tree_start]) {
      const size_t child_idx = tree_start + i;

      um_val child_rate_deriv, child_lik;
      vector<um_val> child_branch_deriv;
      subtree_likelihood_site_deriv(subtree_sizes, rate, Q, expbranches,
                                    trans_mat, child_idx, states, pos, child_lik,
                                    child_rate_deriv, child_branch_deriv);

      const vector<um_mat>::const_iterator P(trans_mat.begin() + child_idx);
      const um_val child_lik_trans(dot(child_lik, P->u), dot(child_lik, P->m));

      // conditional likelihoods
      likelihood.u *= child_lik_trans.u;
      likelihood.m *= child_lik_trans.m;

      // compute the helpers for rate derivatives
      const double expbranch = expbranches[child_idx];
      const um_val expbranch_helper(1.0/expbranch - 1.0, 1.0 - 1.0/expbranch);
      const double common = dot(expbranch_helper, child_lik);
      rate_helper.u += (dot(child_rate_deriv, P->u) + common)/child_lik_trans.u;
      rate_helper.m += (dot(child_rate_deriv, P->m) + common)/child_lik_trans.m;

      // compute the helpers for branch derivatives
      branch_helper[i].u = (dot(child_lik, Q.u)/expbranch)/child_lik_trans.u;
      branch_helper[i].m = (dot(child_lik, Q.m)/expbranch)/child_lik_trans.m;
      for (size_t j = 1; j < subtree_sizes[child_idx]; ++j) {
        branch_helper[i + j].u = dot(child_branch_deriv[j], P->u)/child_lik_trans.u;
        branch_helper[i + j].m = dot(child_branch_deriv[j], P->m)/child_lik_trans.m;
      }
      i += subtree_sizes[child_idx];
    }

    //deriv_rate
    deriv_rate.u = rate_helper.u*likelihood.u;
    deriv_rate.m = rate_helper.m*likelihood.m;

    //deriv_branches
    deriv_branch = vector<um_val>(subtree_sizes[tree_start]);
    for (size_t i = 1; i < subtree_sizes[tree_start]; ++i) {
      deriv_branch[i].u = branch_helper[i].u*likelihood.u;
      deriv_branch[i].m = branch_helper[i].m*likelihood.m;
    }
  }
}


static void
subtree_likelihood_deriv(const vector<size_t> &subtree_sizes,
                         const double &rate, const um_mat &Q,
                         const vector<double> &expbranches,
                         const vector<um_mat> trans_mat,
                         const size_t tree_start,
                         const vector<vector<double> > &states,
                         vector<um_val> &likelihood,
                         vector<um_val> &deriv_rate,
                         vector<vector<um_val> > &deriv_branch) {

  const size_t n_sites = states.front().size();
  const size_t n_branches = states.size();

  if (subtree_sizes[tree_start] == 1) {
    for (size_t i = 0; i < n_sites; ++i) {
      assert(states[tree_start][i] >= 0);
      likelihood.push_back(um_val(states[tree_start][i], 1.0 - states[tree_start][i]));
      deriv_rate.push_back(um_val(0.0, 0.0));
    }
  }
  else {
    vector<um_val> rate_helper(n_sites);
    vector<vector<um_val> > branch_helper(n_sites, vector<um_val>(subtree_sizes[tree_start]));
    likelihood = vector<um_val>(n_sites, um_val(1.0, 1.0));

    size_t i = 1;
    while (i < subtree_sizes[tree_start]) {
      const size_t child_idx = tree_start + i;

      vector<um_val> child_rate_deriv, child_lik;
      vector<vector<um_val> > child_branch_deriv;
      subtree_likelihood_deriv(subtree_sizes, rate, Q, expbranches,
                               trans_mat, child_idx, states, child_lik,
                               child_rate_deriv, child_branch_deriv);

      const vector<um_mat>::const_iterator P(trans_mat.begin() + child_idx);
      const double expbranch = expbranches[child_idx];
      const um_val expbranch_helper(1.0/expbranch - 1.0, 1.0 - 1.0/expbranch);

      for (size_t pos = 0; pos < n_sites; ++pos) {

        const um_val child_lik_trans(dot(child_lik[pos], P->u),
                                     dot(child_lik[pos], P->m));

        // conditional likelihoods
        likelihood[pos].u *= child_lik_trans.u;
        likelihood[pos].m *= child_lik_trans.m;

        // compute the helpers for rate derivatives
        const double common = dot(expbranch_helper, child_lik[pos]);
        rate_helper[pos].u +=
          (dot(child_rate_deriv[pos], P->u) + common)/child_lik_trans.u;
        rate_helper[pos].m +=
          (dot(child_rate_deriv[pos], P->m) + common)/child_lik_trans.m;

        // compute the helpers for branch derivatives
        branch_helper[pos][i].u =
          (dot(child_lik[pos], Q.u)/expbranch)/child_lik_trans.u;
        branch_helper[pos][i].m =
          (dot(child_lik[pos], Q.m)/expbranch)/child_lik_trans.m;
        for (size_t j = 1; j < subtree_sizes[child_idx]; ++j) {
          branch_helper[pos][i + j].u =
            dot(child_branch_deriv[pos][j], P->u)/child_lik_trans.u;
          branch_helper[pos][i + j].m =
            dot(child_branch_deriv[pos][j], P->m)/child_lik_trans.m;
        }
      }
      i += subtree_sizes[child_idx];
    }

    deriv_branch.resize(n_sites, vector<um_val>(n_branches));
    deriv_rate.resize(n_sites);
    for (size_t pos = 0; pos < n_sites; ++pos) {
      //deriv_rate
      deriv_rate[pos].u = rate_helper[pos].u*likelihood[pos].u;
      deriv_rate[pos].m = rate_helper[pos].m*likelihood[pos].m;

      //deriv_branches
      deriv_branch[pos] = vector<um_val>(subtree_sizes[tree_start]);
      for (size_t i = 1; i < subtree_sizes[tree_start]; ++i) {
        deriv_branch[pos][i].u = branch_helper[pos][i].u*likelihood[pos].u;
        deriv_branch[pos][i].m = branch_helper[pos][i].m*likelihood[pos].m;
      }
    }
  }
}



// static void
// subtree_likelihood_site_deriv(const vector<size_t> &subtree_sizes,
//                               const double &rate,
//                               const vector<double> &expbranches,
//                               const size_t tree_start,
//                               const vector<vector<double> > &states,
//                               const size_t &pos,
//                               double &subtree_lik_given_root_unmeth,
//                               double &subtree_lik_given_root_meth,
//                               double &deriv_rate_given_root_unmeth,
//                               double &deriv_rate_given_root_meth,
//                               vector<double> &deriv_branches_given_root_unmeth,
//                               vector<double> &deriv_branches_given_root_meth) {

//   double deriv_rate_helper_u = 0.0;
//   double deriv_rate_helper_m = 0.0;

//   vector<double> deriv_branches_helper_u(subtree_sizes[tree_start], 0.0);
//   vector<double> deriv_branches_helper_m(subtree_sizes[tree_start], 0.0);

//   if (subtree_sizes[tree_start] == 1) {
//     assert(states[tree_start][pos] >= 0);
//     subtree_lik_given_root_unmeth = states[tree_start][pos];
//     subtree_lik_given_root_meth = 1.0 - states[tree_start][pos];
//     deriv_rate_given_root_unmeth = 0.0;
//     deriv_rate_given_root_meth = 0.0;
//   } else {
//     size_t count = 1;
//     subtree_lik_given_root_unmeth = 1.0;
//     subtree_lik_given_root_meth = 1.0;

//     while (count < subtree_sizes[tree_start]) {
//       size_t child_start = tree_start + count;
//       double u, m, u_dr, m_dr;
//       vector<double> u_db, m_db;
//       subtree_likelihood_site_deriv(subtree_sizes, rate, expbranches,
//                                     child_start, states, pos, u, m, u_dr, m_dr,
//                                     u_db, m_db);
//       const double expbranch = expbranches[child_start];
//       vector<vector<double> > child_trans_mat;
//       unitrate_trans_prob_mat(rate, expbranch, child_trans_mat);
//       const double unmeth_child_lik =
//         child_trans_mat[0][0]*u + child_trans_mat[0][1]*m;
//       const double meth_child_lik =
//         child_trans_mat[1][0]*u + child_trans_mat[1][1]*m;

//       //likelihood
//       subtree_lik_given_root_unmeth =
//         subtree_lik_given_root_unmeth*unmeth_child_lik;
//       subtree_lik_given_root_meth =
//         subtree_lik_given_root_meth*meth_child_lik;

//       //deriv_rate
//       deriv_rate_helper_u +=
//         (u_dr*child_trans_mat[0][0] + u*(1.0/expbranch - 1.0) +
//          m_dr*child_trans_mat[0][1] + m*(1.0-1.0/expbranch))/unmeth_child_lik;
//       deriv_rate_helper_m +=
//         (u_dr*child_trans_mat[1][0] + u*(1.0/expbranch - 1.0) +
//          m_dr*child_trans_mat[1][1] + m*(1.0-1.0/expbranch))/meth_child_lik;

//       vector<vector<double> > Q(2, vector<double>(2, 0.0));
//       Q[0][0] = 0.0 - rate;
//       Q[0][1] = rate;
//       Q[1][0] = 1.0 - rate;
//       Q[1][1] = rate - 1.0;

//       deriv_branches_helper_u[count] =
//         (u*Q[0][0] + m*Q[0][1])/expbranch/unmeth_child_lik;
//       deriv_branches_helper_m[count] =
//         (u*Q[1][0] + m*Q[1][1])/expbranch/meth_child_lik;

//       for (size_t i = 1; i < subtree_sizes[child_start]; ++i) {
//           deriv_branches_helper_u[i + count] = 1.0/unmeth_child_lik*
//             (u_db[i]*child_trans_mat[0][0] + m_db[i]*child_trans_mat[0][1]);
//           deriv_branches_helper_m[i + count] = 1.0/meth_child_lik*
//             (u_db[i]*child_trans_mat[1][0] + m_db[i]*child_trans_mat[1][1]);
//       }
//       count += subtree_sizes[child_start];
//     }

//     //deriv_rate
//     deriv_rate_given_root_unmeth =
//       deriv_rate_helper_u*subtree_lik_given_root_unmeth;
//     deriv_rate_given_root_meth =
//       deriv_rate_helper_m*subtree_lik_given_root_meth;

//     //deriv_branches
//     deriv_branches_given_root_unmeth = vector<double>(subtree_sizes[tree_start], 0.0);
//     deriv_branches_given_root_meth = vector<double>(subtree_sizes[tree_start], 0.0);
//     for (size_t i = 1; i < subtree_sizes[tree_start]; ++i) {
//       deriv_branches_given_root_unmeth[i] =
//         deriv_branches_helper_u[i]*subtree_lik_given_root_unmeth;
//       deriv_branches_given_root_meth[i]=
//         deriv_branches_helper_m[i]*subtree_lik_given_root_meth;
//     }
//   }
// }


static void
subtree_likelihood(const vector<size_t> &subtree_sizes,
                   const double &rate,
                   const vector<double> &expbranches,
                   const size_t tree_start,
                   const vector<vector<double> > &states,
                   vector<double> &subtree_lik_given_root_unmeth,
                   vector<double> &subtree_lik_given_root_meth) {

  size_t nsites = states[0].size();
  subtree_lik_given_root_unmeth = vector<double>(nsites, 1.0);
  subtree_lik_given_root_meth = vector<double>(nsites, 1.0);

  if (subtree_sizes[tree_start] == 1) {
    for (size_t pos = 0; pos < nsites; ++ pos) {
      assert (states[tree_start][pos] >= 0);
      subtree_lik_given_root_unmeth[pos] = states[tree_start][pos];
      subtree_lik_given_root_meth[pos] = 1.0 - states[tree_start][pos];
    }
  } else {
    size_t count = 1;
    while (count < subtree_sizes[tree_start]) {
      const size_t child_start = tree_start + count;
      vector<double> u, m;
      subtree_likelihood(subtree_sizes, rate, expbranches,
                         child_start, states, u, m);
      vector<vector<double> > child_trans_mat;
      unitrate_trans_prob_mat(rate, expbranches[child_start], child_trans_mat);
      for (size_t pos = 0; pos < nsites; ++pos) {
        double unmeth_child_lik =
          child_trans_mat[0][0]*u[pos] + child_trans_mat[0][1]*m[pos];
        double meth_child_lik =
          child_trans_mat[1][0]*u[pos] + child_trans_mat[1][1]*m[pos];
        subtree_lik_given_root_unmeth[pos] =
          subtree_lik_given_root_unmeth[pos]*unmeth_child_lik;
        subtree_lik_given_root_meth[pos] =
          subtree_lik_given_root_meth[pos]*meth_child_lik;
      }
      count += subtree_sizes[child_start];
    }
  }
}


static void
subtree_likelihood_site(const vector<size_t> &subtree_sizes,
                        const double &rate,
                        const vector<double> &expbranches,
                        const size_t tree_start,
                        const vector<vector<double> > &states,
                        const size_t &pos,
                        double &subtree_lik_given_root_unmeth,
                        double &subtree_lik_given_root_meth) {

  if (subtree_sizes[tree_start] == 1) {
    assert (states[tree_start][pos] >= 0);
    subtree_lik_given_root_unmeth = states[tree_start][pos];
    subtree_lik_given_root_meth = 1.0 - states[tree_start][pos];
  }
  else {
    subtree_lik_given_root_unmeth = 1.0;
    subtree_lik_given_root_meth = 1.0;
    size_t count = 1;
    while (count < subtree_sizes[tree_start]) {
      const size_t child_start = tree_start + count;
      double u, m;
      subtree_likelihood_site(subtree_sizes, rate, expbranches,
                              child_start, states, pos, u, m);
      vector<vector<double> > child_trans_mat;
      unitrate_trans_prob_mat(rate, expbranches[child_start], child_trans_mat);
      const double unmeth_child_lik =
        child_trans_mat[0][0]*u + child_trans_mat[0][1]*m;
      const double meth_child_lik =
        child_trans_mat[1][0]*u + child_trans_mat[1][1]*m;
      subtree_lik_given_root_unmeth =
        subtree_lik_given_root_unmeth*unmeth_child_lik;
      subtree_lik_given_root_meth =
        subtree_lik_given_root_meth*meth_child_lik;
      count += subtree_sizes[child_start];
    }
  }
}


// ADS: this function is for obtaining derivatives for the purpose of
// estimating parameters??? Needs name change.
static double
tree_loglikelihood_deriv(const vector<size_t> &subtree_sizes,
                         const double root_unmeth_prob,
                         const double rate,
                         const vector<double> &branches,
                         const vector<vector<double> > &states,
                         double &deriv_root_dist,
                         double &deriv_rate,
                         vector<double> &deriv_branches) {
  double likelihood = 0.0;
  const size_t n_sites = states[0].size();
  const size_t n_branches = branches.size();

  vector<double> expbranches;
  for (size_t i = 0; i < n_branches; ++i)
    expbranches.push_back(exp(branches[i]));

  const um_mat Q(um_val(-rate, rate), um_val(1.0 - rate, rate - 1.0));
  const um_val root_prob(root_unmeth_prob, 1.0 - root_unmeth_prob);

  vector<um_mat> trans_mat(n_branches);
  for (size_t i = 1; i < n_branches; ++i)
    unitrate_trans_prob_mat(rate, expbranches[i], trans_mat[i]);

  deriv_root_dist = 0.0;
  deriv_rate = 0.0;
  deriv_branches = vector<double>(n_branches, 0.0);

  bool BYSITE = false;
  if (BYSITE) {
    for (size_t pos = 0; pos < n_sites; ++pos) {

      um_val deriv_rate_given_root;
      um_val likelihood_given_root;
      vector<um_val> deriv_branch_given_root;
      subtree_likelihood_site_deriv(subtree_sizes, rate, Q, expbranches,
                                    trans_mat, 0, states, pos,
                                    likelihood_given_root,
                                    deriv_rate_given_root,
                                    deriv_branch_given_root);

      const double lk = dot(likelihood_given_root, root_prob);
      likelihood += log(lk);

      deriv_root_dist += (likelihood_given_root.u - likelihood_given_root.m)/lk;
      deriv_rate += dot(deriv_rate_given_root, root_prob)/lk;
      for (size_t i = 1; i < n_branches; ++i)
        deriv_branches[i] += dot(deriv_branch_given_root[i], root_prob)/lk;
    }
  }
  else {

    vector<um_val> deriv_rate_given_root;
    vector<um_val> likelihood_given_root;
    vector<vector<um_val> > deriv_branch_given_root;
    subtree_likelihood_deriv(subtree_sizes, rate, Q, expbranches,
                             trans_mat, 0, states, likelihood_given_root,
                             deriv_rate_given_root, deriv_branch_given_root);

    for (size_t pos = 0; pos < n_sites; ++pos) {
      const double lk = dot(likelihood_given_root[pos], root_prob);
      likelihood += log(lk);
      deriv_root_dist += (likelihood_given_root[pos].u -
                          likelihood_given_root[pos].m)/lk;
      deriv_rate += dot(deriv_rate_given_root[pos], root_prob)/lk;
      for (size_t i = 1; i < n_branches; ++i)
        deriv_branches[i] += dot(deriv_branch_given_root[pos][i], root_prob)/lk;
    }

  }
  return likelihood;
}

////////////////////////////////////////////////////////////////////////////////
// Below are wrapper functions for derivatives of each type of parameter  //////
////////////////////////////////////////////////////////////////////////////////

double
tree_loglikelihood_deriv_rootdist(const vector<size_t> &subtree_sizes,
                                  const double root_unmeth_prob,
                                  const double &rate,
                                  const vector<double> &branches,
                                  const vector<vector<double> > &states,
                                  double &deriv_root_dist) {
  double deriv_rate;
  vector<double> deriv_branches;
  double llk = tree_loglikelihood_deriv(subtree_sizes, root_unmeth_prob, rate,
                                        branches, states, deriv_root_dist,
                                        deriv_rate, deriv_branches);
  return llk;
}

double
tree_loglikelihood_deriv_branch(const vector<size_t> &subtree_sizes,
                                const double root_unmeth_prob,
                                const double &rate,
                                const vector<double> &branches,
                                const vector<vector<double> > &states,
                                const size_t which_branch,
                                double &deriv_branch) {
  double deriv_root_dist;
  double deriv_rate;
  vector<double> deriv_branches;
  double llk = tree_loglikelihood_deriv(subtree_sizes, root_unmeth_prob, rate,
                                        branches, states, deriv_root_dist,
                                        deriv_rate, deriv_branches);
  deriv_branch = deriv_branches[which_branch];
  return llk;
}

double
tree_loglikelihood_deriv_rate(const vector<size_t> &subtree_sizes,
                              const double root_unmeth_prob,
                              const double &rate,
                              const vector<double> &branches,
                              const vector<vector<double> > &states,
                              double &deriv_rate) {
  double deriv_root_dist;
  vector<double> deriv_branches;
  double llk = tree_loglikelihood_deriv(subtree_sizes, root_unmeth_prob, rate,
                                        branches, states, deriv_root_dist,
                                        deriv_rate, deriv_branches);
  return llk;
}

////////////////////////////////////////////////////////////////////////////////

double
tree_loglikelihood(const vector<size_t> &subtree_sizes,
                   const double root_unmeth_prob,
                   const double &rate,
                   const vector<double> &branches,
                   const vector<vector<double> > &states) {
  double llk = 0.0;
  size_t nsites = states[0].size();

  vector<double> expbranches;
  for (size_t i = 0; i < branches.size(); ++i) {
    expbranches.push_back(exp(branches[i]));
  }

  bool BYSITE = false;
  if (BYSITE) {
    for (size_t pos = 0; pos < nsites; ++pos) {
      double u, m;
      subtree_likelihood_site(subtree_sizes, rate, expbranches,
                              0, states, pos, u, m);
      llk += log(u*root_unmeth_prob + m*(1.0 - root_unmeth_prob));
    }
  } else {
    vector<double> u, m;
    subtree_likelihood(subtree_sizes, rate, expbranches,
                       0, states, u, m);
    for  (size_t pos = 0; pos < nsites; ++pos) {
      llk += log(u[pos]*root_unmeth_prob + m[pos]*(1.0 - root_unmeth_prob));
    }
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


static void
subtree_pattern_log_prob(const double safecutoff,
                         const double safetol,
                         const vector<size_t> &subtree_sizes,
                         const double rate,
                         const vector<double> &branches,
                         const size_t tree_start,
                         const vector<vector<double> > &states,
                         const size_t pos,
                         vector<double> &subtree_log_probs) {
  if (subtree_sizes[tree_start] == 1) {
    assert (states[tree_start][pos] >= 0);
    double p = states[tree_start][pos];
    if (p < safecutoff) p = safetol;
    if (p > 1.0 - safecutoff) p = 1.0-safetol;
    subtree_log_probs.push_back(log(p));
    subtree_log_probs.push_back(log(1.0 - p));
  } else {
    size_t count = 1;
    vector<vector<double> > pat_log_probs_by_subtree;
    vector<vector<vector<double> > > trans_mats;
    while (count < subtree_sizes[tree_start]) {
      const size_t child_start = tree_start + count;
      const double branch = branches[child_start];
      vector<double> child_log_probs;
      subtree_pattern_log_prob(safecutoff, safetol,
                               subtree_sizes, rate, branches, child_start,
                               states, pos, child_log_probs);
      vector<vector<double> > child_trans_mat;
      trans_prob_mat(rate, 1.0 - rate, branch, child_trans_mat);
      pat_log_probs_by_subtree.push_back(child_log_probs);

      trans_mats.push_back(child_trans_mat);
      count += subtree_sizes[child_start];
    }

    /* combine both subtree probabilities */

    // iterate over the assumed state of the parent
    for (size_t parent_state = 0; parent_state < 2; ++parent_state) {

      // iterate over assumed state of left subtree
      for (size_t i = 0; i < pat_log_probs_by_subtree[0].size(); ++i) {
        size_t left_start = tree_start + 1;
        string left_pattern = dec2binary(subtree_sizes[left_start], i);
        size_t left_root_state = (left_pattern[0] == '0') ? 0 : 1;

        // iterate over assumed state of right subtree
        for (size_t j = 0; j < pat_log_probs_by_subtree[1].size(); ++j) {
          size_t right_start = left_start + subtree_sizes[left_start];
          string right_pattern = dec2binary(subtree_sizes[right_start], j);
          size_t right_root_state = (right_pattern[0] == '0') ? 0 : 1;

          if (pat_log_probs_by_subtree.size()==3) {
            // iterate over assumed state of outgroup
            for (size_t k = 0; k < pat_log_probs_by_subtree[2].size(); ++k) {
              size_t outgroup_start = right_start + subtree_sizes[right_start];
              string outgroup_pattern =
                dec2binary(subtree_sizes[outgroup_start], k);
              size_t outgroup_root_state = (outgroup_pattern[0] == '0') ? 0 : 1;
              double raw_prob =
                log(trans_mats[0][parent_state][left_root_state]) +
                pat_log_probs_by_subtree[0][i] +
                log(trans_mats[1][parent_state][right_root_state]) +
                pat_log_probs_by_subtree[1][j] +
                log(trans_mats[2][parent_state][outgroup_root_state]) +
                pat_log_probs_by_subtree[2][k];

              subtree_log_probs.push_back(raw_prob);
            }
          } else {
            double raw_prob =
              log(trans_mats[0][parent_state][left_root_state]) +
              pat_log_probs_by_subtree[0][i] +
              log(trans_mats[1][parent_state][right_root_state]) +
              pat_log_probs_by_subtree[1][j];

            subtree_log_probs.push_back(raw_prob);
          }
        }
      }
    }
  }
}


double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


static void
tree_pattern_probs_best(const double safecutoff,
                        const double safetol,
                        const vector<size_t> &subtree_sizes,
                        const double root_unmeth_prob,
                        const double rate,
                        const vector<double> &branches,
                        const vector<vector<double> > &states,
                        vector<string> &best_patterns) {

  const size_t nsites = states[0].size();
  for (size_t pos = 0; pos < nsites; ++pos) {
    vector<double> tree_log_probs;
    subtree_pattern_log_prob(safecutoff, safetol,
                             subtree_sizes, rate, branches, 0, states, pos,
                             tree_log_probs);
    vector<double>::iterator result;
    result = std::max_element(tree_log_probs.begin(), tree_log_probs.end());
    size_t best = std::distance(tree_log_probs.begin(), result);
    string pattern = dec2binary(subtree_sizes[0], best);
    best_patterns.push_back(pattern);
  }
}


static void
tree_pattern_probs(const double safecutoff,
                   const double safetol,
                   const vector<size_t> &subtree_sizes,
                   const double root_unmeth_prob,
                   const double rate,
                   const vector<double> &branches,
                   vector<vector<double> > &states) {

  const size_t nsites = states[0].size();
  const size_t treesize = subtree_sizes[0];
  for (size_t pos = 0; pos < nsites; ++pos) {
    vector<double> tree_log_probs;
    subtree_pattern_log_prob(safecutoff, safetol,
                             subtree_sizes, rate, branches, 0, states, pos,
                             tree_log_probs);
    for (size_t i = 0; i < tree_log_probs.size(); ++i) {
      string pattern = dec2binary(subtree_sizes[0], i);
      if (pattern[i]=='0') {
        tree_log_probs[i] = log(root_unmeth_prob) + tree_log_probs[i];
      } else {
        tree_log_probs[i] = log(1.0 - root_unmeth_prob) + tree_log_probs[i];
      }
    }
    const double magic_positive_val = 1.0;
    vector<double> unmeth_log_probs_by_node(treesize, magic_positive_val);
    vector<double> meth_log_probs_by_node(treesize, magic_positive_val);
    for (size_t i = 0; i < tree_log_probs.size(); ++i) {
      string pattern = dec2binary(subtree_sizes[0], i);
      for (size_t j = 0; j < treesize; ++j) {
        if (pattern[j] == '0') {
          unmeth_log_probs_by_node[j] =
            log_sum_log(tree_log_probs[i], unmeth_log_probs_by_node[j]);
        } else {
          meth_log_probs_by_node[j] =
            log_sum_log(tree_log_probs[i], meth_log_probs_by_node[j]);
        }
      }
    }

    for (size_t i = 0; i < treesize; ++i) {
      if (subtree_sizes[i] >1 )
        states[i][pos] = exp(unmeth_log_probs_by_node[i] -
                             log_sum_log(unmeth_log_probs_by_node[i],
                                         meth_log_probs_by_node[i]));
    }
  }
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////// CODE FOR OPTIMIZING PARAMETERS BELOW /////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


inline double
sign(const double val) {
  if (val > 0) return 1.0;
  if (val < 0) return -1.0;
  return 0.0;
}


double
optimize_rootdist(const bool VERBOSE,
                  const double LIKTOL,
                  const double PARTOL,
                  const double STEPSIZE,
                  const vector<size_t> &subtree_sizes,
                  const double &rate, const vector<double> &branches,
                  const vector<vector<double> > &states,
                  double &root_unmeth_prob) {

  double deriv_root_u;
  double prev_llk =
    tree_loglikelihood_deriv_rootdist(subtree_sizes, root_unmeth_prob, rate,
                                      branches, states, deriv_root_u);
  double prev_val = root_unmeth_prob;
  double fac = 1.0; // modifies step size in case of failure to find
                    // a valid (i.e. better) next value
  double new_val;
  do {
    new_val = prev_val + STEPSIZE*sign(deriv_root_u)*fac;
    fac = fac*0.5;
  } while (new_val < PARTOL || new_val + PARTOL > 1.0);

  double new_llk =
    tree_loglikelihood_deriv_rootdist(subtree_sizes, new_val, rate, branches,
                                      states, deriv_root_u);

  while (new_llk < prev_llk && abs(new_val - prev_val) > PARTOL) {
    fac = fac*0.5;
    new_val = max(PARTOL,
                  min(1.0 - PARTOL,
                      prev_val + STEPSIZE*sign(deriv_root_u)*fac));
    new_llk =
      tree_loglikelihood_deriv_rootdist(subtree_sizes, new_val, rate,
                                        branches, states, deriv_root_u);
  }

  if (VERBOSE)
    cerr << "param = " << new_val << "\t"<< new_llk
         << "\tImprove = " << new_llk - prev_llk << endl;

  if (new_llk > prev_llk) {
    prev_val = new_val;
    prev_llk = new_llk;
  }

  root_unmeth_prob = prev_val;
  return prev_llk;
}


double
optimize_rate(const bool VERBOSE,
              const double LIKTOL,
              const double PARTOL,
              const double STEPSIZE,
              const vector<size_t> &subtree_sizes,
              const double root_unmeth_prob, const vector<double> &branches,
              const vector<vector<double> > &states, double &rate) {

  double deriv_rate;
  double prev_llk =
    tree_loglikelihood_deriv_rate(subtree_sizes, root_unmeth_prob, rate,
                                  branches, states, deriv_rate);
  double prev_rate = rate;
  double new_rate;
  double new_llk;
  double fac = 1.0;
  do {
    new_rate = prev_rate + STEPSIZE*sign(deriv_rate)*fac;
    fac = fac*0.5;
  } while (new_rate < PARTOL || new_rate > 1.0 - PARTOL);

  new_llk =
    tree_loglikelihood_deriv_rate(subtree_sizes, root_unmeth_prob, new_rate,
                                  branches, states, deriv_rate);
  while (new_llk < prev_llk && abs(new_rate - prev_rate) > PARTOL) {
    fac = fac*0.5;
    new_rate = max(PARTOL, min(prev_rate + STEPSIZE*sign(deriv_rate)*fac,
                                  1.0 - PARTOL));
    new_llk =
      tree_loglikelihood_deriv_rate(subtree_sizes, root_unmeth_prob,
                                    new_rate, branches, states, deriv_rate);
  }

  if (VERBOSE)
    cerr << "param = " << new_rate << "\t"<< new_llk
         << "\tImprove = " << new_llk - prev_llk << endl;

  if (new_llk > prev_llk) {
    prev_llk = new_llk;
    prev_rate = new_rate;
  }

  rate = prev_rate;
  return prev_llk;
}


double
optimize_branch(const bool VERBOSE,
                const double LIKTOL,
                const double PARTOL,
                const double STEPSIZE,
                const vector<size_t> &subtree_sizes,
                const double root_unmeth_prob,
                const vector<vector<double> > &states,
                const double &rate, const size_t which_branch,
                vector<double> &branches) {

  double deriv_branch;
  double prev_llk =
    tree_loglikelihood_deriv_branch(subtree_sizes, root_unmeth_prob, rate,
                                    branches, states, which_branch,
                                    deriv_branch);
  double prev_val = branches[which_branch];
  double new_val;
  double new_llk;
  double fac = 1.0;
  do {
    new_val = prev_val + STEPSIZE*sign(deriv_branch)*fac;
    fac = fac*0.5;
  } while (new_val < PARTOL);

  branches[which_branch] = new_val;
  new_llk =
    tree_loglikelihood_deriv_branch(subtree_sizes, root_unmeth_prob, rate,
                                    branches, states, which_branch,
                                    deriv_branch);
  while (new_llk < prev_llk && abs(new_val - prev_val) > PARTOL) {
    fac = fac*0.5;
    new_val = max(PARTOL,
                  prev_val + STEPSIZE*sign(deriv_branch)*fac);
    branches[which_branch] = new_val;
    new_llk =
      tree_loglikelihood_deriv_branch(subtree_sizes, root_unmeth_prob, rate,
                                      branches, states, which_branch,
                                      deriv_branch);
  }

  if (VERBOSE)
    cerr << "param = " << new_val << "\t"<< new_llk
         << "\tImprove = " << new_llk-prev_llk << endl;

  branches[which_branch] = prev_val;
  return prev_llk;
}

double 
optimize_branches(const bool VERBOSE,
                  const double LIKTOL,
                  const double PARTOL,
                  const double STEPSIZE,
                  const vector<size_t> &subtree_sizes,
                  const double root_unmeth_prob,
                  const vector<vector<double> > &states,
                  const double &rate,
                  vector<double> &branches) {
  double deriv_root_dist;
  double deriv_rate;
  vector<double> deriv_branches;
  double prev_llk = tree_loglikelihood_deriv(subtree_sizes, root_unmeth_prob, rate,
                                             branches, states, deriv_root_dist,
                                             deriv_rate, deriv_branches);
  double denom = 0.0;
  for (size_t i = 0; i < deriv_branches.size(); ++i) {
    denom += abs(deriv_branches[i]);
  }

  vector<double> prev_branches = branches;
  vector<double> new_branches = branches;
  double new_llk;
  double fac = 1.0;
  bool OUTRANGE = false;
  do {
    OUTRANGE = false;
    for (size_t i = 1; i < branches.size(); ++i) {
      new_branches[i] = prev_branches[i] + deriv_branches[i]*fac/denom;
      if (new_branches[i] < PARTOL)
        OUTRANGE = true;
    }
    fac = fac*0.5;
  } while (OUTRANGE);

  double new_deriv_root_dist;
  double new_deriv_rate;
  vector<double> new_deriv_branches;
  new_llk =
    tree_loglikelihood_deriv(subtree_sizes, root_unmeth_prob, rate,
                             new_branches, states, new_deriv_root_dist,
                             new_deriv_rate, new_deriv_branches);
  while (new_llk < prev_llk && fac > PARTOL) {
    cerr << fac << endl;
    fac = fac*0.5;
    for (size_t i = 1; i < branches.size(); ++i) {
      new_branches[i] = prev_branches[i] + deriv_branches[i]*fac/denom;
    }
    new_llk =
      tree_loglikelihood_deriv(subtree_sizes, root_unmeth_prob, rate,
                               new_branches, states, new_deriv_root_dist,
                               new_deriv_rate, new_deriv_branches);
  }

  if (new_llk > prev_llk) {
    branches = new_branches;
  }

  if (VERBOSE) {
    for (size_t i = 1; i < branches.size(); ++i) {
      cerr << "[branch_" << i <<"]\t"<< branches[i] << endl;
    }
    cerr << "[log-likelihood]\t" << new_llk << "\tImprove = " << new_llk-prev_llk << endl;
  }
  
  return (new_llk > prev_llk)? new_llk: prev_llk;
}


double
optimize(const bool VERBOSE, const size_t MAXITER, const double LIKTOL,
         const double PARTOL, const double stepsize,
         const vector<size_t> &subtree_sizes,
         double &root_unmeth_prob,
         double &lam,
         vector<double> &branches,
         vector<vector<double> > &states){

  double llk = tree_loglikelihood(subtree_sizes, root_unmeth_prob, lam,
                                  branches, states);
  if (VERBOSE) cerr << "BEGIN llk=" << llk << endl;
  double prev_llk;
  for (size_t iter = 0; iter < MAXITER; ++iter) {
    if (VERBOSE)
      cerr << "-----------------Iteration " << iter
           << "-----------------" << endl;
    prev_llk = llk;

    llk = optimize_rootdist(true, LIKTOL, PARTOL, stepsize, subtree_sizes,
                            lam, branches, states, root_unmeth_prob);
    if (VERBOSE)
      cerr << "[root_unmeth_prob]\t" << root_unmeth_prob
           << "\t[log-likelihood]\t" << llk << endl;

    llk = optimize_rate(true, LIKTOL, PARTOL, stepsize, subtree_sizes,
                        root_unmeth_prob, branches, states, lam);
    if (VERBOSE)
      cerr << "[rate]\t" << lam
           << "\t[log-likelihood]\t" << llk << endl;

    // for (size_t i = 1; i < branches.size(); ++i) {
    //   llk = optimize_branch(true, LIKTOL, PARTOL, stepsize, subtree_sizes,
    //                         root_unmeth_prob, states, lam, i, branches);
    //   if (VERBOSE)
    //     cerr << "[branch_" << i <<"]\t"<< branches[i]
    //          << "\t[log-likelihood]\t" << llk << endl;
    // }

    llk = optimize_branches(true, LIKTOL, PARTOL, stepsize,
                            subtree_sizes, root_unmeth_prob, states,
                            lam, branches);

    if (VERBOSE)
      cerr << "[Improved] "<< llk- prev_llk << "(TOL=" << LIKTOL << ")"<< endl;
    if (abs(llk - prev_llk) < LIKTOL) {
      if (VERBOSE)
        cerr << "Converged at iteration " << iter << endl;
      break;
    }
  }

  if (VERBOSE) {
    cerr << "MLE parameters:" << endl
         << "[root_unmeth_prob]\t" << root_unmeth_prob<< endl
         << "[rate]\t" << lam << endl;
    for (size_t i = 1; i < branches.size(); ++i)
      cerr << "[branch" << i << "]\t" << branches[i] << endl;
    cerr << "[Log-likelihood]" << "\t" << llk << endl;
  }
  return llk;
}


static void
write_states(std::ostream &out,
             const PhyloTreePreorder &t,
             const double lam,
             const double root_unmeth_prob,
             const double llk,
             const vector<Site> &sites,
             const vector<vector<double> > &states) {
  const size_t nsites = sites.size();
  const size_t n_species = states.size();
  out << "#" << t.Newick_format() << "\tRate = " << lam << ";"
      << "\tRoot unmeth prob = " << root_unmeth_prob
      << "\tlog-likelihood = " << llk << endl;
  for (size_t i = 0; i < nsites; ++i) {
    out << sites[i].chrom << ":" << sites[i].pos << "\t";
    for (size_t j = 0; j < n_species; ++j)
      out << ((states[j][i] > 0.5) ? "0" : "1");
    out << endl;
  }
}


static void
write_patterns(std::ostream &out,
               const PhyloTreePreorder &t,
               const double lam,
               const double root_unmeth_prob,
               const double llk,
               const vector<Site> &sites,
               const vector<string> &patterns) {
  const size_t nsites = sites.size();
  out << "#" << t.Newick_format() << "\tRate = " << lam << ";"
      << "\tRoot unmeth prob = " << root_unmeth_prob
      << "\tlog-likelihood = " << llk << endl;
  for (size_t i = 0; i < nsites; ++i) {
    out << sites[i].chrom << ":" << sites[i].pos << "\t"
        << patterns[i] << endl;
  }
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
//////////////////////          MAIN             ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int
main(int argc, const char **argv) {

  try {
    // run mode flags
    bool VERBOSE = false;
    bool COUNT = false;
    bool NODEMAP = false;

    double LIKTOL = 1e-2;
    size_t MAXITER = 20;
    double PARTOL = 1e-4;
    double STEP_SIZE = 0.1;

    string outfile;
    string paramfile;
    /************************* COMMAND LINE OPTIONS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "test functionality for "
                           "estimating phylo-epigenomic parameters ",
                           "assuming independent sites"
                           "<newick> <meth-tab>");
    opt_parse.add_opt("counts", 'c', "meth-tab contains read counts instead"
                      "of state probability (default: false)", false, COUNT);
    opt_parse.add_opt("params", 'p', "use parameters from file "
                      "and skip optimization", false, paramfile);
    opt_parse.add_opt("iteration", 'i', "max iteration (default: 10)",
                      false, MAXITER);
    opt_parse.add_opt("nodemap", 'n', "output  MAP states at each node"
                      "instead of most likely methylation patterns"
                      "(default: false)", false, NODEMAP);
    opt_parse.add_opt("out", 'o', "output file (default: stdout)",
                      false, outfile);
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
    /**************************************************************************/


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
    vector<vector<pair<size_t, size_t> > > meth_count_table;
    vector<vector<double> > meth_prob_table;
    if (COUNT) {
      load_meth_table(meth_table_file, sites, meth_count_table,
                      meth_table_species);
    }
    else {
      load_meth_table(meth_table_file, sites, meth_prob_table,
                      meth_table_species);
    }
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
    /**********  OPTIMIZATION *************************************************/
    const size_t nsites = sites.size();
    const size_t nleaf = meth_table_species.size();
    vector<vector<double> > states(subtree_sizes.size(),
                                   vector<double>(nsites, -1.0));
    double root_unmeth_prob = 0.0;
    for (size_t i = 0; i < meth_table_species.size(); ++i) {
      vector<double> statepost;
      if (COUNT) {
        meth_to_posterior(meth_count_table, i, statepost, root_unmeth_prob);
      } else {
        meth_to_posterior(meth_prob_table, i, statepost, root_unmeth_prob);
      }
      states[species_in_order[i]] = statepost;
    }
    root_unmeth_prob = root_unmeth_prob/(nleaf*nsites);


    for (size_t i = 0; i < 10; ++i) {
      for (size_t j = 0; j < meth_table_species.size(); ++j) {
        cerr << states[species_in_order[j]][i] << "\t";
      }
      cerr << endl;
    }


    double llk;
    double lam;
    if (paramfile.empty()) {
      // initialize rate between (0,1)
      lam = 0.4;
      llk = optimize(VERBOSE, MAXITER, LIKTOL, PARTOL, STEP_SIZE,
                     subtree_sizes, root_unmeth_prob, lam, branches, states);
      t.set_branch_lengths(branches);
    } else {
      read_params(VERBOSE, paramfile, root_unmeth_prob, lam, t);
      t.get_branch_lengths(branches);
      llk = tree_loglikelihood(subtree_sizes, root_unmeth_prob, lam,
                               branches, states);
    }

    /**************************************************************************/
    /**********  SET INTERNAL NODES STATES ************************************/
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    const double safecutoff = 0.3;
    const double safetol = 1e-4;
    if (!NODEMAP) {
      if (VERBOSE)
        cerr << "[COMPUTING MOST LIKELY PATTERNS]" << endl;
      vector<string> bestpatterns;
      tree_pattern_probs_best(safecutoff, safetol,
                              subtree_sizes, root_unmeth_prob, lam,
                              branches, states, bestpatterns);
      write_patterns(out, t, lam, root_unmeth_prob, llk, sites, bestpatterns);
    } else {
      if (VERBOSE)
        cerr << "[COMPUTING MAXIMUM A POSTERIORI STATES AT EACH NODE]" << endl;
      tree_pattern_probs(safecutoff, safetol,
                         subtree_sizes, root_unmeth_prob,
                         lam, branches, states);
      write_states(out, t, lam, root_unmeth_prob, llk, sites, states);
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


/* ADS:
   xx Move the stuff in main into their own functions
   -- Consider a wrapper for the optimization stuff
   xx Make variable names more sensible and typically longer
   xx No const references for primitive types
   -- Make sure ordering of parameters to functions makes sense and is
      always consistent
   -- What are the semantics of TOLERANCE? Should it be both a
      precision value for estimates and also used for convergence?
   xx Make sure convergence criteria is based on only one value
 */

/* TO-DO:
   xx Let size of deriv_branches be variable according to the subtree size.
 */
