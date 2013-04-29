#include "BinaryTree.hpp"

#include <cstddef>

BinaryTreeNode::BinaryTreeNode(){
  child1_ = NULL;
  child2_ = NULL;
  parent_ = NULL;
  next_ = NULL;
  prev_ = NULL;
  nChildren_ = 0;
}

BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode& node):
  child1_( node.child1_), child2_( node.child1_), nChildren_(node.nChildren_),
  parent_( node.parent_), next_( node.next_), prev_(node.prev_)
{}

BinaryTreeNode::~BinaryTreeNode(){

}

void BinaryTreeNode::setChildren(BinaryTreeNode* child1, 
				 BinaryTreeNode* child2){
  child1_ = child1;
  child2_ = child2;
  nChildren_ = child1->nChildren_ + child2->nChildren_;
  child1->parent_ = this;
  child2->parent_ = this;
}

bool BinaryTreeNode::hasChildren(){
  if(child1_ != NULL || child2_ != NULL){
    return true;
  }
  return false;
}

int BinaryTreeNode::getChildrenCount(){
  return nChildren_;
}

BinaryTreeNode* BinaryTreeNode::getChild1(){
  return child1_;
}

BinaryTreeNode* BinaryTreeNode::getChild2(){
  return child2_;
}

void BinaryTreeNode::setChild1( BinaryTreeNode* child ){
  child1_ = child;
}

void BinaryTreeNode::setChild2( BinaryTreeNode* child ){
  child2_ = child;
}

BinaryTreeNode* BinaryTreeNode::getParent(){
  return parent_;
}

BinaryTree::BinaryTree(){
  root_ = new BinaryTreeNode();
}

BinaryTree::~BinaryTree(){
  deleteTree();
}

void BinaryTree::addLeaf( BinaryTreeNode* newLeaf){
  
  bool isNodeAdded = false;
  while(!isNodeAdded){
    BinaryTreeNode* node = root_; // current node being considered
    if (node->hasChildren()){
      node->nChildren_++;
      if( node->getChildrenCount() % 2 == 0){
	// even number, choose child 1
	node = node->getChild1();
      }else{
	// odd number, choose child 2
	node = node->getChild2();
      }
    }else{
      // current node is an existing leaf
      // make the current node a parent and the existing leaf as a child
      BinaryTreeNode* parent = node;
      BinaryTreeNode* existingLeaf = new BinaryTreeNode( *node );
      parent->setChildren(existingLeaf, newLeaf);
    }
  }
}

void BinaryTree::deleteTree(){
  
  bool isTreeEmpty = root_->hasChildren();
  BinaryTreeNode* node = root_; // the current node
  BinaryTreeNode* parent = node;
  while( !isTreeEmpty ){

    parent = node;
    if( parent->getChild1() != NULL ){
      node = parent->getChild1();
      parent->setChild1(NULL);
    }else if(parent->getChild2() != NULL ){
      node = parent->getChild2();
      parent->setChild2(NULL);
    }
    
    if(node->hasChildren() == false){
      BinaryTreeNode* node_temp = node;
      node = node->getParent();
      delete node_temp;

      if(node == root_ && root_->hasChildren() == false){
	isTreeEmpty = true;
      }
    }

  }

  delete root_;

}
