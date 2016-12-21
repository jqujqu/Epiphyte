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

double param_set::tolerance = 1e-4; // tolerance for parameter estimate precision

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
  in >> t;
  in >> pi0;
  in >> rate0;
  in >> g0;
  in >> g1;

  // sync the transformed values for brances in parameter set
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
