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

#include "param_set.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using std::vector;
using std::string;

double param_set::tolerance = 5e-4;  // tolerance for parameter estimate precision

std::ostream &
operator<<(std::ostream &out, const param_set &ps) {
  return out << ps.tostring();
}

void
param_set::read(const string &paramfile, PhyloTreePreorder &t) {

  std::ifstream in(paramfile.c_str());
  if (!in)
    throw std::runtime_error("cannot read: " + paramfile);

  // First read in the tree (Newick format). The reason this is done
  // here is to ensure the format of the file is correct. When
  // parameters are provided and not learned, the tree will be read as
  // usual, from the tree file, but also replicated here
  string dummy_label;
  in >> dummy_label >> t;
  in >> dummy_label >> pi0;
  in >> dummy_label >> rate0;

  // read in G parameters: The last pair of values will be repeated for the
  // remaining branches if number of pairs is less than tree-size
  for (size_t i = 0; i < t.get_size(); ++i) {
    double val0, val1;
    in >> dummy_label >> val0;
    in >> dummy_label >> val1;
    if (val0 > 0.0 && val1 > 0.0 && val0 < 1.0 && val1 > 0.0) {
      g0.push_back(val0);
      g1.push_back(val1);
    } else if (g0.size() > 0 && g1.size() > 0) {
      g0.push_back(g0.back());
      g1.push_back(g1.back());
    } else { // will fail assertion later
      g0.push_back(val0);
      g1.push_back(val1);
    }
  }

  // sync the transformed values for branches in parameter set
  vector<double> branches;
  t.get_branch_lengths(branches);
  T.clear();
  T.resize(branches.size(), 0.0);
  for (size_t i = 1; i < branches.size(); ++i)
    // ADS: need some way to check that this transformation has
    // happened and has been reversed when appropriate.
    T[i] = 1.0 - 1.0/exp(branches[i]);

  assert(is_valid());
}

void
param_set::write(PhyloTreePreorder t, const string &paramfile) {

  std::ofstream out(paramfile.c_str());
  if (!out)
    throw std::runtime_error("cannot write: " + paramfile);

  // first set the tree branches
  vector<double> branches(T.size(), 0.0);
  for (size_t i = 1; i < T.size(); ++i)
    branches[i] = -log(1.0 - T[i]);
  t.set_branch_lengths(branches);

  out << "tree\t" << t << std::endl;
  out << "pi0\t" << pi0 << std::endl;
  out << "rate0\t" << rate0 << std::endl;
  for (size_t i = 0; i < T.size(); ++i) {
    out << "g0[" << i << "]\t" << g0[i] << std::endl;
    out << "g1[" << i << "]\t" << g1[i] << std::endl;
  }
}

string
param_set::tostring() const {
  std::ostringstream oss;
  oss << "pi0=" << pi0 << ", "
      << "rate0=" << rate0 << ", ";
  for (size_t i = 0; i < T.size(); ++i) 
    oss << "g0[" << i << "]=" << g0[i] << ", "
        << "g1[" << i << "]=" << g1[i] << ", ";
  oss << "T=(";
  for (size_t i = 0; i < T.size() - 1; ++i)
    oss << -log(1.0 - T[i]) << ',';
  oss << -log(1.0 - T.back()) << ')';
  return oss.str();
}


void
param_set::assign_branches(PhyloTreePreorder &t) const {
  vector<double> branches(T.size());
  for (size_t i = 0; i < T.size(); ++i)
    branches[i] = -log(1.0 - T[i]);
  t.set_branch_lengths(branches);
}
