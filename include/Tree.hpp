/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef TREE_HPP
#define TREE_HPP

#include <cstddef>
#include <vector>

namespace rfs{

/** 
 * \class TreeNode
 * A generic tree node
 * \brief A generic tree node
 */
class TreeNode{

public:

  /** Constructor 
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  TreeNode(size_t n_children_exp = -1);
  
  /** Destructor, virtual becuase other classes will be derived from this class*/
  virtual ~TreeNode();

  /** 
   * Get a pointer to the parent node
   * \return pointer to parent, or NULL if this node has no parent
   */
  TreeNode* getParent();
  
  /**
   * Get pointers to a child
   * \param[in] child index
   * \return a pointer to the child
   */ 
  TreeNode* getChild(size_t idx);

  /**
   * Get the number of children that this node has
   * \return number of children
   */
  size_t getChildrenCount();

  /**
   * Add a new child to this node. 
   * \return pointer to new child node
   */
  TreeNode* addChild();

private:

  std::vector< boost::shared_ptr<TreeNode> > children_;
  TreeNode* parent_;
};

}

#endif
