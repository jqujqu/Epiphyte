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

#ifndef PARAM_SET_HPP
#define PARAM_SET_HPP

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//////////////////// FOR HOLDING PARAMETER SETS ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

#include "PhyloTreePreorder.hpp"
#include "PhyloTree.hpp"

#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <cmath>
#include <iterator>

struct param_set {
  static double tolerance;  // tolerance for parameter estimate precision
  double pi0;
  double rate0;
  std::vector<double> g0;
  std::vector<double> g1;
  std::vector<double> T; // Definition: T[i] = 1.0 - 1.0/exp(branches[i]);

  param_set() {}

  param_set(double _pi0, double _rate0,
            const std::vector<double> &_g0, const std::vector<double> &_g1,
            const std::vector<double> &_T) :
    pi0(_pi0), rate0(_rate0), g0(_g0), g1(_g1), T(_T) {}

  param_set(double _pi0, double _rate0,
            const std::vector<double> &_g0, const std::vector<double> &_g1) :
    pi0(_pi0), rate0(_rate0), g0(_g0), g1(_g1) {}

  void assign_branches(PhyloTreePreorder &t) const;

  static double
  absolute_difference(const param_set &a, const param_set &b) {
    assert(a.T.size() == b.T.size());
    double d =  std::abs(a.pi0 - b.pi0) + std::abs(a.rate0 - b.rate0);
    for (size_t i = 0; i < a.T.size(); ++i) {
      d += std::abs(a.g0[i] - b.g0[i]) + std::abs(a.g1[i] - b.g1[i]);
      d += std::abs(a.T[i] - b.T[i]);
    }
    return d;
  }

  static double
  max_abs_difference(const param_set &a, const param_set &b) {
    assert(a.T.size() == b.T.size());
    double d = std::abs(a.pi0 - b.pi0);
    d = std::max(std::abs(a.rate0 - b.rate0), d);
    for (size_t i = 0; i < a.T.size(); ++i) {
      d = std::max(std::abs(a.g0[i] - b.g0[i]), d);
      d = std::max(std::abs(a.g1[i] - b.g1[i]), d);
      d = std::max(std::abs(a.T[i] - b.T[i]), d);
    }
    return d;
  }

  std::string tostring() const;

  void read(const std::string &paramfile, PhyloTreePreorder &t);

  void write(const PhyloTreePreorder t, const std::string &paramfile);

  bool is_valid() const { // strict inequalities to avoid invalid log value
    bool valid = pi0 < 1.0 && pi0 > 0.0 && rate0 < 1.0 && rate0 > 0.0;
    for (size_t i = 0; i < T.size(); ++i) {
      valid = valid && (g0[i] < 1.0 && g0[i] > 0.0 &&
                        g1[i] < 1.0 && g1[i] > 0.0);
      if (i > 0)
        valid = valid && T[i] < 1.0 && T[i] > 0.0; 
    }
    return valid;
  }
};

std::ostream &
operator<<(std::ostream &out, const param_set &ps);

#endif
