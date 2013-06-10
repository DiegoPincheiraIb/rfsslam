#include "BinaryTree.hpp"

BinaryTreeNode::BinaryTreeNode(){
  children_[0] = NULL;
  children_[1] = NULL;
  parent_ = NULL;
  next_ = NULL;
  prev_ = NULL;
  nLeaves_ = 0;
}

BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode* node){
  children_[0] = node->children_[0]; 
  children_[1] = node->children_[1]; 
  nLeaves_ = node->nLeaves_;
  parent_ = node->parent_; 
  next_ = node->next_; 
  prev_ = node->prev_;
}

BinaryTreeNode::~BinaryTreeNode(){

}

void BinaryTreeNode::setChildren(BinaryTreeNode* child0, 
				 BinaryTreeNode* child1){

  if( child0 != NULL){
    children_[0] = child0;
    children_[0]->parent_ = this;
  }
  if( child1 != NULL){
    children_[1] = child1;
    children_[1]->parent_ = this;
  }

  // update number of descendents for this node and its ascendents
  updateLeafCount();
  BinaryTreeNode* parent = parent_;
  while( parent != NULL ){
    parent->updateLeafCount();
    parent = parent->parent_;
  }

}

void BinaryTreeNode::getChildren(BinaryTreeNode* child0, 
				 BinaryTreeNode* child1){
  if( child0 != NULL){
    child0 = children_[0];
  }
  if( child1 != NULL){
    child1 = children_[1];
  }
}

void BinaryTreeNode::removeChildren(bool child0, bool child1){
  if(child0 == true)
    children_[0] = NULL;
  if(child1 == true)
    children_[1] = NULL;
  updateLeafCount();
}

BinaryTreeNode* BinaryTreeNode::getParent(){
  return parent_;
}

unsigned int BinaryTreeNode::getLeafCount(){
  return nLeaves_;
}

unsigned int BinaryTreeNode::updateLeafCount(){
  nLeaves_ = 0;
  if(children_[0] != NULL)
    nLeaves_ += children_[0]->nLeaves_;
  if(children_[1] != NULL)
    nLeaves_ += children_[1]->nLeaves_;  
  return nLeaves_;
}




////////// Binary Tree Implementation //////////

BinaryTree::BinaryTree(){
  root_ = new BinaryTreeNode();
}

BinaryTree::~BinaryTree(){
  deleteTree();
}

void BinaryTree::addLeaf( BinaryTreeNode* newLeaf){

  newLeaf->nLeaves_ = 1;

  BinaryTreeNode* node = root_; // current node being considered
  BinaryTreeNode* nodeParent = NULL; 
  unsigned int nLeaves0 = 0; // number of leaves under branch 0
  unsigned int nLeaves1 = 0; // number of leaves under branch 1

  // Traverse tree to find the node (with least children) to add the leaf to 
  while(node != NULL){
    
    nLeaves0 = 0; // number of leaves under branch 0
    nLeaves1 = 0; // number of leaves under branch 1

    if(node->children_[0] != NULL){
      nLeaves0 = node->children_[0]->getLeafCount();
    }
    if(node->children_[1] != NULL){
      nLeaves1 = node->children_[0]->getLeafCount();
    }

    nodeParent = node;

    if( nLeaves0 == 0 && nLeaves1 == 0){
      
      // current node is a leaf
      // split into a parent + curent leaf + new leaf
      
      BinaryTreeNode* currentLeaf = new BinaryTreeNode(node);
      nodeParent->setChildren(currentLeaf, newLeaf);
      node = NULL;
    }else if( nLeaves0 <= nLeaves1 ){ // go down branch 0
      
      node = node->children_[0];
      if(node == NULL){ // nodeParent is not a leaf but only has 1 child on branch 1
	nodeParent->setChildren( newLeaf, NULL );
      }

    }else{ // go down branch 1
      node = node->children_[1];
      if(node == NULL){ // nodeParent is not a leaf but only has 1 child on branch 0
	nodeParent->setChildren( NULL, newLeaf );
      }
    }

  }
    
}

BinaryTreeNode* BinaryTree::deleteTree( BinaryTreeNode* startNode){
  
  if(startNode == NULL){
    startNode = root_;
  }

  BinaryTreeNode* node = startNode; // the current node
  BinaryTreeNode* parent = NULL;
  BinaryTreeNode* child = NULL;

  while(1){

    if( node->children_[0] != NULL ){ // go down branch 0
      parent = node;
      node = node->children_[0];
    }else if(node->children_[1] != NULL ){ // go down branch 1
      parent = node;
      node = node->children_[1];
    }else{ // no children, delete node
      BinaryTreeNode* node_to_delete = node;
      
      // udpate the parent of node_to_delete
      node = node->parent_;
      if( node->children_[0] == node_to_delete )
	node->removeChildren(true, false);
      else if( node->children_[1] == node_to_delete )
	node->removeChildren(false, true); 
      node->updateLeafCount();
      delete node_to_delete;
      
      if( node_to_delete == startNode ){ // finished
	return node; // returns the parent of the node that we just deleted
      }
    }
  }
}
