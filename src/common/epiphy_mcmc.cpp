/* Copyright (C) 2015-16 University of Southern California and
 *                       Andrew D. Smith and Jenny Qu
 *
 * Authors: Jenny Qu and Andrew Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 */

#include "epiphy_mcmc.hpp"
#include "sufficient_statistics_helpers.hpp"

#include <vector>
#include <random>
#include <iostream>
#include <math.h>       /* pow */


using std::vector;
using std::cerr;
using std::endl;

////////////////////////////////////////////////////////////////////////////////
/////////////////////   Measure MCMC convergence    ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void
mcmc_stat::scale() {
  double tot = root_start_distr.first + root_start_distr.second;
  assert(tot > 0);
  root_start_distr.first = root_start_distr.first/tot;
  root_start_distr.second = root_start_distr.second/tot;

  tot = root_distr.uu + root_distr.um + root_distr.mu + root_distr.mm;
  assert(tot > 0);
  root_distr.div(tot);

  assert(start_distr.size() > 1 && triad_distr.size() > 1);
  for (size_t i = 1; i < start_distr.size(); ++i) {
    start_distr[i].to_probabilities();
    tot = (triad_distr[i].uuu + triad_distr[i].uum +
           triad_distr[i].umu + triad_distr[i].umm +
           triad_distr[i].muu + triad_distr[i].mum +
           triad_distr[i].mmu + triad_distr[i].mmm);
    assert(tot > 0);
    triad_distr[i].div(tot);
  }
}


static void
flatten_mcmcstat(const mcmc_stat &m,
                 vector<double> &v) {
  v.clear();
  v.push_back(m.root_start_distr.first);
  v.push_back(m.root_start_distr.second);

  vector<double> p;
  m.root_distr.flatten(p);
  v.insert(v.end(), p.begin(), p.end());

  assert(m.start_distr.size() > 1 && m.triad_distr.size() > 1);
  for (size_t i = 1; i < m.start_distr.size(); ++i) {
    m.start_distr[i].flatten(p);
    v.insert(v.end(), p.begin(), p.end());
  }

  for (size_t i = 1; i < m.triad_distr.size(); ++i) {
    m.triad_distr[i].flatten(p);
    v.insert(v.end(), p.begin(), p.end());
  }
}


void
within_mean(const vector<vector<vector<double> > > &y,
            vector<vector<double> > &y_within_mean) {
  const size_t n_chain = y.size();
  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();

  // #chains x #variables
  y_within_mean = vector<vector<double> >(n_chain, vector<double>(n_var, 0.0));
  for (size_t i = 0; i < n_chain; ++i) {
    for (size_t j = 0; j < n_var; ++j) {
      for (size_t k = 0; k < n_samp; ++k) {
        y_within_mean[i][j] += y[i][k][j];
      }
      y_within_mean[i][j] /= n_samp;
    }
  }
}


void
overall_mean(const vector<vector<vector<double> > > &y,
             vector<double> &y_overall_mean) {

  const size_t n_chain = y.size();
  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();
  y_overall_mean = vector<double>(n_var, 0.0); // #variables

  for (size_t i = 0; i < n_var; ++i) {
    for (size_t j = 0; j < n_chain; ++j) {
      for (size_t k = 0; k < n_samp; ++k) {
        y_overall_mean[i] += y[j][k][i];
      }
    }
    y_overall_mean[i] /= n_samp*n_chain;
  }
}


void
within_variance(const vector<vector<vector<double> > > &y,
                const vector<vector<double> > &y_within_mean,
                vector<double> &y_within_variance) {
  const size_t n_chain = y.size();
  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();

  y_within_variance = vector<double>(n_var);

  for (size_t i = 0; i < n_var; ++i) {
    for (size_t j = 0; j < n_chain; ++j) {
      for (size_t k = 0; k < n_samp; ++k) {
        y_within_variance[i] += pow(y[j][k][i] - y_within_mean[j][i], 2);
      }
    }
    y_within_variance[i] /= n_chain*(n_samp - 1);
  }
}


void
btwn_variance(const vector<vector<double> > &y_within_mean,
              const vector<double> &y_overall_mean,
              const size_t n_samp,
              vector<double> &y_btwn_variance) {
  const size_t n_chain = y_within_mean.size();
  const size_t n_var = y_within_mean[0].size();
  y_btwn_variance = vector<double>(n_var);
  for (size_t i = 0; i < n_var; ++i) {
    for (size_t j = 0; j < n_chain; ++j) {
      y_btwn_variance[i] += pow(y_within_mean[j][i] - y_overall_mean[i], 2);
    }
    y_btwn_variance[i] = y_btwn_variance[i]*n_samp/(n_chain - 1);
  }
}


void
EPSR(vector<vector<mcmc_stat> > &mcmcstats,
     vector<double> &epsr) {
  epsr.clear();

  // cerr << "[in EPSR]";

  vector<vector<vector<double> > > y; // #chains x #samples x #variables
  for (size_t i = 0; i < mcmcstats.size(); ++i) { // i-th chain
    vector<vector<double> > chain;
    for (size_t j = 0; j < mcmcstats[0].size(); ++j) { // j-th sample
      vector<double> stat;
      flatten_mcmcstat(mcmcstats[i][j], stat);
      chain.push_back(stat);
    }
    y.push_back(chain);
  }

  const size_t n_samp = y[0].size();
  const size_t n_var = y[0][0].size();

  vector<vector<double> > y_within_mean; // #chain x #var
  within_mean(y, y_within_mean);

  vector<double> y_overall_mean; //#var
  overall_mean(y, y_overall_mean);

  vector<double> y_within_variance; //
  within_variance(y, y_within_mean, y_within_variance);

  vector<double> y_btwn_variance;
  btwn_variance(y_within_mean, y_overall_mean, n_samp, y_btwn_variance);

  vector<double> vhat(n_var);
  for (size_t i = 0; i < n_var; ++i) {
    vhat[i] = (y_within_variance[i]*(n_samp - 1) + y_btwn_variance[i])/n_samp;
  }

  for (size_t i = 0; i < n_var; ++i) {
    epsr.push_back(pow(vhat[i]/y_btwn_variance[i], 0.5));
  }
}


void
sum(const vector<mcmc_stat> &mcmcstats,
    mcmc_stat &ave_mcmc_stat) {
  ave_mcmc_stat = mcmcstats[0];
  for (size_t i = 1; i < mcmcstats.size(); ++i) {
    ave_mcmc_stat.root_start_distr.first += mcmcstats[i].root_start_distr.first;
    ave_mcmc_stat.root_start_distr.second += mcmcstats[i].root_start_distr.second;
    ave_mcmc_stat.root_distr += mcmcstats[i].root_distr;
    for (size_t j = 0; j < mcmcstats[0].start_distr.size(); ++j) {
      ave_mcmc_stat.start_distr[j] += mcmcstats[i].start_distr[j];
      ave_mcmc_stat.triad_distr[j] += mcmcstats[i].triad_distr[j];
    }
  }
}
