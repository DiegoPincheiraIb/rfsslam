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
 * \brief A generic tree node, defined using curiously recurring template pattern (CRTP)
 * \tparam Derived the derived class
 */
template<class Derived>
class TreeNode{

public:

  /** Constructor 
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  TreeNode(size_t n_children_exp = 0);
  
  /** Destructor, virtual becuase other classes will be derived from this class*/
  virtual ~TreeNode();

  /** 
   * Get a pointer to the parent node
   * \return pointer to parent, or NULL if this node has no parent
   */
  Derived* getParent();
  
  /**
   * Get pointers to a child
   * \param[in] child index
   * \return a pointer to the child
   */ 
  Derived* getChild(size_t idx);

  /**
   * Get the number of children that this node has
   * \return number of children
   */
  size_t getChildrenCount();

  /**
   * Add a new child to this node. 
   * \return pointer to new child node
   */
  Derived* addChild();

private:

  std::vector< boost::shared_ptr<Derived> > children_;
  Derived* parent_;
};

}

// Implementation


#include <boost/shared_ptr.hpp>
#include "Tree.hpp"

namespace rfs{

  template<class Derived>
  TreeNode<Derived>::TreeNode(size_t n_children_exp):parent_(NULL){
    if(n_children_exp > 0)
      children_.resize(n_children_exp);
  }
  
  template<class Derived>
  TreeNode<Derived>::~TreeNode(){}
  
  template<class Derived>
  Derived* TreeNode<Derived>::getParent(){
    return parent_;
  }
  
  template<class Derived>
  Derived* TreeNode<Derived>::getChild(size_t idx){
    if(idx < children_.size())
      return children_[idx].get();
    else
      return NULL;
  }

  template<class Derived>
  size_t TreeNode<Derived>::getChildrenCount(){
    return children_.size();
  }

  template<class Derived>
  Derived* TreeNode<Derived>::addChild(){
    boost::shared_ptr<Derived> child(new Derived);
    child->parent_ = static_cast<Derived*>(this);
    children_.push_back(child);
    return child.get();
  }

}



#endif
