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

#ifndef PHYLOTREE_PREORDER_HPP
#define PHYLOTREE_PREORDER_HPP

#include <string>
#include <vector>

#include "PhyloTree.hpp"

using std::string;
using std::vector;

class PhyloTreePreorder : public PhyloTree {
public:
  void
  get_subtree_sizes(vector<size_t> &subtree_sizes) const {
    get_subtree_sizes(root, subtree_sizes);
  }
  void
  get_branch_lengths(vector<double> &branch_lengths) const {
    get_branch_lengths(root, branch_lengths);
  }
  void
  get_leaf_names(vector<string> &names) const {
    get_leaf_names(root, names);
  };

  void
  set_branch_lengths(vector<double> branch_lengths) {
    set_branch_lengths(root, branch_lengths);
  }

private:
  void
  get_subtree_sizes(const PhyloTree::PTNode &node,
                    vector<size_t> &subtree_sizes) const;
  void
  get_branch_lengths(const PhyloTree::PTNode &node,
                     vector<double> &branch_lengths) const;
  void
  get_leaf_names(const PhyloTree::PTNode &node,
                 vector<string> &names) const;
  void
  set_branch_lengths(PhyloTree::PTNode &node,
                     vector<double> &branch_lengths);
};

#endif
