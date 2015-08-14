/*    Copyright (C) 2015 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Jenny Qu
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

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "PhyloTreePreorder.hpp"

using std::string;
using std::vector;


void
PhyloTreePreorder::get_subtree_sizes_implementation(const PhyloTree::PTNode &node,
                                                    vector<size_t> &subtree_sizes) const {
  subtree_sizes.push_back(1);

  if (node.has_children()) {
    for (size_t i = 0; i < node.child.size(); ++i) {
      vector<size_t> child_subtree_sizes;
      get_subtree_sizes_implementation(node.child[i], child_subtree_sizes);
      subtree_sizes.front() += child_subtree_sizes.size();
      copy(child_subtree_sizes.begin(), child_subtree_sizes.end(),
           back_inserter(subtree_sizes));
    }
  }
}


void
PhyloTreePreorder::get_branch_lengths_implementation(const PhyloTree::PTNode &node,
                                                     vector<double> &branch_lengths) const {
  branch_lengths.push_back(node.branch_length);
  for (size_t i = 0; i < node.child.size(); ++i)
    get_branch_lengths_implementation(node.child[i], branch_lengths);
}


void
PhyloTreePreorder::get_leaf_names_implementation(const PhyloTree::PTNode &node,
                                                 vector<string> &leafnames) const {
  if (node.child.size()==0) {
    leafnames.push_back(node.name);
  } else {
    for (size_t i = 0; i < node.child.size(); ++i)
      get_leaf_names_implementation(node.child[i], leafnames);
  }
}

void
PhyloTreePreorder::set_branch_lengths_implementation(PhyloTree::PTNode &node,
                                                     vector<double> &branch_lengths) {
  node.branch_length = branch_lengths[0];
  branch_lengths.erase(branch_lengths.begin());
  for (size_t i = 0; i < node.child.size(); ++i)
    set_branch_lengths_implementation(node.child[i], branch_lengths);
}
