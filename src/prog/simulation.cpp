/*****************************************************************************
 *  simulation: a program to simulate methylation states for species
 *    in a phylogenetic tree according to inheritance and
 *    auto-correlation of methylation states, and WGBS data sets for
 *    extant species.
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Jenny Qu
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
 *****************************************************************************/



#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <sstream>
#include <tr1/unordered_set>
#include <algorithm> //std::max

#include <gsl/gsl_randist.h>
#include <unistd.h>
#include <time.h>       /* time */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "PhyloTreePreorder.hpp"


using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::pair;
using std::make_pair;
using std::tr1::unordered_set;
using std::max;

static bool
is_G_line(const char *line){
  static const char *label = "G-MATRIX";
  static const size_t label_len = 8;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
}

static bool
is_Q_line(const char *line){
  static const char *label = "Q-MATRIX";
  static const size_t label_len = 8;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
 }

static bool
is_tree_line(const char *line){
  static const char *label = "NEWICKTREE";
  static const size_t label_len = 10;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
}

static bool
is_pi_line(const char *line) {
  static const char *label = "PI0";
  static const size_t label_len = 3;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
}

static bool
is_cov_line(const char *line){
  static const char *label = "COV-PARAM";
  static const size_t label_len = 9;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
}

static bool
is_fg_line(const char *line){
  static const char *label = "FG-BETA-PARAMS";
  static const size_t label_len = 14;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
}

static bool
is_bg_line(const char *line){
  static const char *label = "BG-BETA-PARAMS";
  static const size_t label_len = 14;
  for (size_t i = 0; i < label_len; ++i)
    if (line[i] != label[i])
      return false;
  return true;
}

void
parse_paramfile(const string param_file,
                double &pi0,
                vector<double> &G, vector<double> &Q,
                size_t &coverage, vector<double> &F_PARAM,
                vector<double> &B_PARAM,
                PhyloTreePreorder &t){
  G = vector<double>(2, 0);
  Q = vector<double>(2, 0);
  F_PARAM = vector<double>(2, 0);
  B_PARAM = vector<double>(2, 0);
  static const size_t buffer_size = 10000; // Magic
  std::ifstream in(param_file.c_str());
  if (!in)
    throw SMITHLABException("bad file: " + param_file);
  while (!in.eof()){
    char buffer[buffer_size];
    in.getline(buffer, buffer_size);
    if (in.gcount() == buffer_size - 1)
      throw SMITHLABException("Line too long in file: " + param_file);
    std::istringstream is(buffer);
    string tmp;
    if (is_tree_line(buffer)) {
      is >> tmp >> t;
    } else if (is_pi_line(buffer)) {
      is >> tmp >> pi0;
    } else if (is_G_line(buffer)) {
      is >> tmp >> G[0] >> G[1];
    } else if (is_Q_line(buffer)) {
      is >> tmp >> Q[0] >> Q[1];
    } else if (is_cov_line(buffer)) {
      is >> tmp >> coverage;
    } else if (is_fg_line(buffer)) {
      is >> tmp >> F_PARAM[0] >> F_PARAM[1];
    } else if (is_bg_line(buffer)) {
      is >> tmp >> B_PARAM[0] >> B_PARAM[1];
    }
  }
}


static void
load_cpgs(const string &cpgs_file,
          vector<SimpleGenomicRegion> &cpgs){

  vector<GenomicRegion> cpgs_in;
  ReadBEDFile(cpgs_file, cpgs_in);
  assert(check_sorted(cpgs_in));
  if (!check_sorted(cpgs_in))
    throw SMITHLABException("regions not sorted in file: " + cpgs_file);

  for (size_t i = 0; i < cpgs_in.size(); ++i) {
    cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
  }
}

bool bernoulli(const double p) {
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  long seed = rand();
  gsl_rng_set (r, seed);
  size_t result = gsl_ran_bernoulli(r, p);
  gsl_rng_free(r);
  return result==1;
}


//use boolean varaible for methylation state
// true=foreground=u=hypomethylated
// false=background=m=hypermethylated
void sim_root_state(const bool prev,
                    const vector<double> &G,
                    bool &new_state) {
  double p = prev? G[0]: G[1];
  new_state = bernoulli(p)? prev : !prev;
}

//NTP: Normalized Transition Probability
void sim_newstate(const bool prev, const bool ancestor,
                  const vector<double> &G,
                  const vector<double> &Q,
                  const double branch,
                  bool &new_state){
  double tol = 1e-10;
  double NTP = 0; //probability of getting foreground
  double t = Q[0]+Q[1];
  double h = exp(-t*branch);
  double p0, p1;
  if(ancestor){
    p0 = (Q[0]*h + Q[1])/t;
    p1 = Q[0]*(1-h)/t;
    if(p0+p1 !=1){
      assert (fabs(1-p0-p1) < tol);
      p0 =p0/(p0+p1);
      p1 = 1 - p0;
    }
  }else{
    p0 = Q[1]*(1-h)/t;
    p1 = (Q[0] + Q[1]*h)/t;
    if(p0+p1 !=1){
      assert (fabs(1-p0-p1) < tol);
      p0 =p0/(p0+p1);
      p1 = 1 - p0;
    }
  }
  double total;
  if(prev){
    total = G[0]*p0 + (1-G[0])*p1;
    NTP = G[0]*p0 /total;
  }else{
    total = (1-G[1])*p0 + G[1]*p1;
    NTP = (1-G[1])*p0/total;
  }
  new_state = bernoulli(NTP);
}

//simulate first positon in a segment
void sim_newstate(const bool ancestor,
                  const vector<double> &Q,
                  const double branch,
                  bool &new_state){
  double t = Q[0]+Q[1];
  double h = exp(-t*branch);
  double p0; //prob of transition into "u" state
  if(ancestor){
    p0 = (Q[0]*h + Q[1])/t;
  }else{
    p0 = Q[1]*(1-h)/t;
  }
  new_state = bernoulli(p0);
}


//simulate position in  preorder
void
simulate_position(const bool first,
                  const vector<bool> &prev_states,
                  const double pi0,
                  const vector<double> &G,
                  const vector<double> &Q,
                  const vector<size_t> &subtree_sizes,
                  const vector<size_t> &tree_parent_index,
                  const vector<double> &branches,
                  vector<bool> &states){
  states = vector<bool> (subtree_sizes.size(), false);
  bool new_state;
  //Begin simulation from the root
  if (first) {
    bool initstate = bernoulli(pi0);
    sim_root_state(initstate, G, new_state);
  } else {
    sim_root_state(prev_states[0], G, new_state);
  }
  states[0] = new_state;

  for (size_t i = 1; i < subtree_sizes.size(); ++i ) {
    double b = branches[i];
    bool ancestor_state = states[tree_parent_index[i]];
    if (first) {
      sim_newstate(ancestor_state, Q, b, new_state);
    } else {
      sim_newstate(prev_states[i], ancestor_state, G, Q, b, new_state);
    }
    states[i] = new_state;
  }
}


// HME: true = hypo = 0 = T;
//      false = hyper = 1 = C;
void
printHME(const vector<bool> &HME) {
  for (size_t i=0; i < HME.size(); ++i) {
    if(HME[i]) cerr << 0 ;
    else cerr << 1;
  }
  cerr <<endl;
}

void
simulate_counts(const size_t NB_r, const double NB_p,
                const vector<double> &beta_params,
                double &meth, size_t &N) {
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  long seed = rand();
  gsl_rng_set(r, seed);
  N = gsl_ran_negative_binomial (r, NB_p, NB_r);
  if (N == 0) N = 1;
  double pm = gsl_ran_beta (r, beta_params[0], beta_params[1]);
  meth = static_cast<double>(gsl_ran_binomial (r, pm, N))/N;
  gsl_rng_free(r);
}

bool
write_site(std::ostream &out,
           const string &chrom, const size_t &pos,
           const string &strand, const string &seq,
           const double &meth, const size_t &coverage) {
  return (out << chrom << "\t" << pos << "\t" << strand
          << "\t" << seq << "\t" << (coverage == 0 ? 0.0 : meth) << "\t"
          << coverage << '\n');
}

bool
write_states(std::ostream &out,
             const string &chrom, const size_t &pos,
             const string &strand, const vector<bool> &HME) {

  string seq="";
  for(size_t i =0; i < HME.size(); ++i){
    if(HME[i]) seq +="T";
    else seq += "C";
  }
  return (out << chrom << "\t" << pos << "\t" << strand
          << "\t" << seq << '\n');
}

// hypoprob = 1 <=> true, hypo, state="0", "T"
string
HME_to_hypoprob(const vector<size_t> subtree_sizes,
                const vector<bool> &HME) {
  string s;
  for (size_t i = 0; i < subtree_sizes.size(); ++i) {
    if (subtree_sizes[i] == 1) {
      s+="\t";
      s += (HME[i]) ? "1.0" : "0.0";
    }
  }
  return s;
}

bool
write_hypoprob(std::ostream &out,
               const string &chrom, const size_t &pos,
               const vector<size_t> &subtree_sizes,
               const vector<bool> HME) {
  return (out << chrom <<":" << pos << HME_to_hypoprob(subtree_sizes, HME)
          << '\n' );
}

static void
get_parent_index(const vector<size_t> &subtree_sizes,
                 vector<size_t> &tree_parent_index) {
  tree_parent_index = vector<size_t>(subtree_sizes.size(), 0);
  for (size_t i = 0; i < subtree_sizes.size(); ++i) {
    if (subtree_sizes[i] > 1) {
      size_t count =1;
      while (count < subtree_sizes[i]) {
        size_t child = i+count;
        tree_parent_index[child] = i;
        count += subtree_sizes[child];
      }
    }
  }
}


int main(int argc, const char **argv) {

  try {
    bool VERBOSE = false;
    bool DEBUG = false;
    string outfile, cpgs_file;
    size_t desertsize = 1000;
    bool SINGLE = false;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate methylomes "
                           "according to a phylogenetic tree",
                           "<parameter file>");
    opt_parse.add_opt("CpG", 'c', "methcount file or bed file",
                      true, cpgs_file);
    opt_parse.add_opt("singlesite", 's', "simulate sites independently",
                      false, SINGLE);
    opt_parse.add_opt("desert", 'd', "desert size (default 1000)",
                      false, desertsize);
    opt_parse.add_opt("output", 'o', "name of output file "
                      "(default: stdout)", false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false, VERBOSE);
    opt_parse.add_opt("debug", 'D', "print Debug info",
                      false, DEBUG);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string param_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in(param_file.c_str());
    if (!in)
      throw SMITHLABException("bad file: " + param_file);

    string tree_rep;
    vector<double> G, Q;
    size_t coverage;
    double pi0;
    vector<double> F_PARAM, B_PARAM;
    if(VERBOSE)
      cerr << "----Reading parameters----" << endl;
    PhyloTreePreorder t;
    parse_paramfile(param_file, pi0, G, Q, coverage, F_PARAM, B_PARAM, t);
    vector<size_t> subtree_sizes;
    t.get_subtree_sizes(subtree_sizes);
    vector<size_t> tree_parent_index;
    get_parent_index(subtree_sizes,tree_parent_index);
    vector<string> leaf_names;
    t.get_leaf_names(leaf_names);
    vector<double> branches;
    t.get_branch_lengths(branches);


    vector<SimpleGenomicRegion> cpgs;
    vector<pair<double, double> > meths;
    vector<size_t> reads;

    if(VERBOSE)
      cerr << "----Loading cpgs----" << endl;
    load_cpgs(cpgs_file, cpgs);


    // simulate position by position
    srand(time(NULL));
    vector<vector<bool> > HMEs;
    vector<bool> prev_HME, HME;
    bool first;

    for(size_t j = 0; j < cpgs.size(); ++j){
      first = (j == 0 ||
               (j > 0 && cpgs[j].distance(cpgs[j-1]) > desertsize));
      if (SINGLE) first = true;
      prev_HME = HME;
      simulate_position(first, prev_HME, pi0, G, Q, subtree_sizes,
                        tree_parent_index, branches, HME);
      HMEs.push_back(HME);
      if (DEBUG) printHME(HME);
    }

    if(VERBOSE){
       cerr << "Leaf nodes are:" ;
      for (size_t i = 0; i < leaf_names.size(); ++i)
        cerr << leaf_names[i] << "\t";
      cerr << endl;
    }

    vector<size_t> leafidx;
    for(size_t i = 0; i < subtree_sizes.size(); ++i){
      if (subtree_sizes[i] == 1)
        leafidx.push_back(i);
    }

    double meth;
    size_t N;
    string methstate;
    for(size_t i = 0; i < leaf_names.size(); ++i){
      //simulate counts and print
      string sim_outfile = outfile + '_' + leaf_names[i];
      if(VERBOSE)
        cerr <<" Writing to " << sim_outfile << endl;

      std::ofstream of;
      if (!outfile.empty()) of.open(sim_outfile.c_str());
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

      for(size_t j = 0; j < cpgs.size(); ++j){
        if(HMEs[j][leafidx[i]]) {
          simulate_counts(coverage, 0.5, F_PARAM, meth, N);
          methstate = "T";
        }
        else{
          simulate_counts(coverage, 0.5, B_PARAM, meth, N);
          methstate = "C";
        }
        write_site(out, cpgs[j].get_chrom(), cpgs[j].get_start(),
                   "+", methstate, meth, N);
      }
      of.close();
    }

    string state_outfile = outfile + "_treestates";
    std::ofstream of;
    if (!state_outfile.empty()) of.open(state_outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if(VERBOSE) cerr << "Writing to " << state_outfile << endl;
    out << "##" << t.Newick_format() << endl;
    for(size_t j = 0; j < cpgs.size(); ++j)
      write_states(out, cpgs[j].get_chrom(), cpgs[j].get_start(),
                   "+", HMEs[j]);
    of.close();

    string hypoprob_outfile = outfile + "_hypoprobs";
    if (!hypoprob_outfile.empty()) of.open(hypoprob_outfile.c_str());
    std::ostream hout(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if(VERBOSE) cerr << "Writing to " << hypoprob_outfile << endl;

    for (size_t i = 0; i < leaf_names.size(); ++i) {
      out << leaf_names[i] << "\t";
    }
    out << "\n";
    for(size_t j = 0; j < cpgs.size(); ++j)
      write_hypoprob(hout, cpgs[j].get_chrom(), cpgs[j].get_start(),
                     subtree_sizes, HMEs[j]);
    of.close();


  }catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
