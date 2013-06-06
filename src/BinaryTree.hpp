// Binary Tree Class
// Keith Leung 2013

#ifndef BINARY_TREE_HPP
#define BINARY_TREE_HPP

#include <cstddef>


class BinaryTree;

/** 
 *  \class BinaryTreeNode
 *  \brief A class for the nodes of a binary tree
 */
class BinaryTreeNode
{
public:
  
  /** Default constructor */
  BinaryTreeNode();

  /** Copy constructor */
  BinaryTreeNode( const BinaryTreeNode* node);
  
  /** Default destructor */
  ~BinaryTreeNode();

  friend class BinaryTree;
  
  /** 
   * Set the direct children of this node
   * Each child will also have its parent updated to this node
   * \param child0 pointer to child node 0 (use NULL to not set)
   * \param child1 pointer to child node 1 (use NULL to not set)
   */
  void setChildren(BinaryTreeNode* child0=NULL, BinaryTreeNode* child1=NULL);

  /**
   * Get the direct children of this node
   * \param child0 pointer to child node 0 (use NULL to not get)
   * \param child1 pointer to child node 1 (use NULL to not get)
   */
  void getChildren(BinaryTreeNode* child0=NULL, BinaryTreeNode* child1=NULL); 

  /**
   * Remove the pointer to one or both children of this node
   * Note that memory does no get deallocated for the removed node[s]
   * \param child0 True if child0 is to be removed
   * \param child1 True if child1 is to be removed
   */
  void removeChildren(bool child0 = false, bool child1 = false);

  /**
   *  Get the point to the parent node
   *  \return pointer to parent or NULL if parent does not exist
   */
  BinaryTreeNode* getParent();

  /** 
   *  Get the number of descendents
   *  \return number of descendents
   */
  unsigned int getLeafCount();


private:

  unsigned int nLeaves_; /**< number of leaves*/
  BinaryTreeNode* children_ [2]; /**< pointers to children nodes */
  BinaryTreeNode* parent_; /**< pointer to parent node */
  BinaryTreeNode* next_; /**< pointer to the next leaf node */
  BinaryTreeNode* prev_; /**< pointer to the previous leaf node */

  /** 
   *  Update the number of descendents
   *  \return the updated count
   */
  unsigned int updateLeafCount();
  
};


/** 
 *  \class BinaryTree
 *  \brief A class for constructing and manipulating a binary tree
 *  \todo single leaf removal 
 *  \todo branch replacement 
 */
class BinaryTree
{

public:

  /** Default constructor */
  BinaryTree();
  
  /** Default destructor */
  ~BinaryTree();

  /** 
   *  Add a leaf to the tree while maintaining balance
   *  \param leaf pointer to the node to be added
   */
  void addLeaf( BinaryTreeNode* leaf);

  /** 
   *  Delete a node and all its children, with memory deallocation 
   *  \param node to delete. Defaults to the root node (i.e., delete the entire tree)
   *  \return parent of the deleted node, NULL if entire tree was deleted
   */
  BinaryTreeNode* deleteTree( BinaryTreeNode* startNode = NULL);

private:

  BinaryTreeNode* root_; /**< root node of the binary tree */
};

#endif
