/*    Copyright (C) 2016 University of Southern California and
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

#ifndef SUFFICIENT_STATISTICS_HELPERS_HPP
#define SUFFICIENT_STATISTICS_HELPERS_HPP

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////////  STRUCTS FOR HOLDING PAIRS AND TRIPLES OF STATES  ////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <cmath>
#include <sstream>

struct param_set;

struct pair_state {
  double uu, um; // think 2x2 matrix
  double mu, mm;

  pair_state(double _uu, double _um, double _mu, double _mm) :
    uu(_uu), um(_um), mu(_mu), mm(_mm) {}
  pair_state() : uu(0.0), um(0.0), mu(0.0), mm(0.0) {}

  // WARNING: no range checking (for i,j > 1)
  double operator()(int i, int j) const { // version for l-values
    return (i == 0) ? (j == 0 ? uu : um) : (j == 0 ? mu : mm);
  }
  double & operator()(int i, int j) { // version for r-values
    return (i == 0) ? (j == 0 ? uu : um) : (j == 0 ? mu : mm);
  }
  pair_state operator+(const pair_state &other) const {
    return pair_state(uu + other.uu, um + other.um,
                      mu + other.mu, mm + other.mm);
  }
  void div(const double x) {
    uu /= x; um /= x;
    mu /= x; mm /= x;
  }
  void make_logs() {
    uu = std::log(uu); um = std::log(um);
    mu = std::log(mu); mm = std::log(mm);
  }
  void flatten(std::vector<double> &p) const {
    p.clear();
    p.push_back(uu); p.push_back(um);
    p.push_back(mu); p.push_back(mm);
  }
  std::string tostring() const {
    std::ostringstream oss;
    oss << "[" << uu << ", " << um << "]\n"
        << "[" << mu << ", " << mm << "]";
    return oss.str();
  }
};

std::ostream &
operator<<(std::ostream &out, const pair_state &ps);

struct triple_state {
  double uuu, uum; // think: a 2x2 matrix
  double umu, umm;

  double muu, mum; // think: another 2x2 matrix
  double mmu, mmm;

  triple_state() : uuu(0.0), uum(0.0), umu(0.0), umm(0.0),
                   muu(0.0), mum(0.0), mmu(0.0), mmm(0.0) {};
  triple_state(double _uuu, double _uum, double _umu, double _umm,
               double _muu, double _mum, double _mmu, double _mmm) :
    uuu(_uuu), uum(_uum), umu(_umu), umm(_umm),
    muu(_muu), mum(_mum), mmu(_mmu), mmm(_mmm) {};

  // WARNING: no range checking
  double operator()(int i, int j, int k) const { // version for r-value (rhs)
    return  i == 0 ?
      (j == 0 ? (k == 0 ? uuu : uum) : (k == 0 ? umu : umm)) :
      (j == 0 ? (k == 0 ? muu : mum) : (k == 0 ? mmu : mmm));
  }
  double & operator()(int i, int j, int k) { // version for l-value (lhs)
    return  i == 0 ?
      (j == 0 ? (k == 0 ? uuu : uum) : (k == 0 ? umu : umm)) :
      (j == 0 ? (k == 0 ? muu : mum) : (k == 0 ? mmu : mmm));
  }
  triple_state operator+(const triple_state &other) const {
    return triple_state(uuu + other.uuu, uum + other.uum,
                        umu + other.umu, umm + other.umm,
                        muu + other.muu, mum + other.mum,
                        mmu + other.mmu, mmm + other.mmm);
  }
  void div(const double x) {
    uuu /= x; uum /= x;
    umu /= x; umm /= x;

    muu /= x; mum /= x;
    mmu /= x; mmm /= x;
  }
  void make_logs() {
    uuu = std::log(uuu); uum = std::log(uum);
    umu = std::log(umu); umm = std::log(umm);

    muu = std::log(muu); mum = std::log(mum);
    mmu = std::log(mmu); mmm = std::log(mmm);
  }
  void flatten(std::vector<double> &p) const {
    p.clear();
    p.push_back(uuu); p.push_back(uum);
    p.push_back(umu); p.push_back(umm);
    p.push_back(muu); p.push_back(mum);
    p.push_back(mmu); p.push_back(mmm);
  }
  std::string tostring() const {
    std::ostringstream oss;
    oss << "u [" << uuu << ", " << uum << "]\n"
        << "  [" << umu << ", " << umm << "]\n"
        << "m [" << muu << ", " << mum << "]\n"
        << "  [" << mmu << ", " << mmm << "]";
    return oss.str();
  }
};

std::ostream &
operator<<(std::ostream &out, const triple_state &ts);

void
get_transition_matrices(const param_set &ps,
                        std::vector<pair_state> &P,
                        std::vector<triple_state> &GP);

void
get_transition_matrices_deriv(const param_set &ps,
                              std::vector<pair_state> &P,
                              std::vector<triple_state> &GP,
                              std::vector<triple_state> &GP_drate,
                              std::vector<triple_state> &GP_dg0,
                              std::vector<triple_state> &GP_dg1,
                              std::vector<triple_state> &GP_dT);


#endif
