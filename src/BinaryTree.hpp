#ifndef BINARY_TREE_HPP
#define BINARY_TREE_HPP

class BinaryTree;

/** 
 *  \class BinaryTreeNode
 *  \brief A class for the nodes of a binary tree
 *  \todo template on data type (i.e., map feature)
 *  \todo functions for next and previous pointers 
 */
class BinaryTreeNode
{
public:
  
  /** Default constructor */
  BinaryTreeNode();

  /** Copy constructor */
  BinaryTreeNode( const BinaryTreeNode& node);
  
  /** Default destructor */
  ~BinaryTreeNode();

  friend class BinaryTree;
  
  /** 
   *  Set the direct children of a node. 
   *  Each child will also have its parent updated to this node
   *  \param child1 pointer to child node 1
   *  \param child2 pointer to child node 2
   */
  void setChildren(BinaryTreeNode* child1, BinaryTreeNode* child2);

  /**
   *  Check if the node has children
   *  \return true or false
   */
  bool hasChildren();
  
  /** 
   *  Get the number of descendents / children (direct and indirect)
   *  \return number of descendents
   */
  int getChildrenCount();

  /**
   *  Get the pointer to child 1
   *  \return pointer to child 1 or NULL if child 1 does not exist
   */
  BinaryTreeNode* getChild1();

  /**
   *  Get the pointer to child 2
   *  \return pointer to child 2 or NULL if child 2 does not exist
   */
  BinaryTreeNode* getChild2();

  /**
   *  Set the pointer to child 1
   *  \param child the node to point child1 to
   */
  void setChild1(BinaryTreeNode* child );

  /**
   *  Set the pointer to child 2
   *  \param child the node to point child2 to
   */
  void setChild2(BinaryTreeNode* child );

  /**
   *  Get the point to the parent node
   *  \return pointer to parent or NULL if parent does not exist
   */
  BinaryTreeNode* getParent();


private:

  int nChildren_; /**< number of descendents / children (direct and indirect) */
  BinaryTreeNode* child1_; /**< pointer to child node 1 */
  BinaryTreeNode* child2_; /**< pointer to child node 2 */
  BinaryTreeNode* parent_; /**< pointer to parent node */
  BinaryTreeNode* next_; /**< pointer to the next leaf node */
  BinaryTreeNode* prev_; /**< pointer to the previous leaft node */
  
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
   *  Delete the tree and deallocate memory used by each node
   */
  void deleteTree();

private:

  BinaryTreeNode* root_; /**< root node of the binary tree */
};

#endif
