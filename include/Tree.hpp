#ifndef TREE_HPP
#define TREE_HPP

#include <cstddef>
#include <list>

class Node{

public:

  /** Constructor 
   * \param[in] id The id number for this node. Default is -1
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  Node(int n_children_exp = -1):parent_(NULL), nChildren_(0){
    if(n_children_exp > 0)
      children_.resize(n_children_exp);
  }
  
  /** Destructor */
  ~Node(){
    Node* parent = NULL;
    Node* current = this;
    do{
      if(current->nChildren_ > 0){
	parent = current;
	current = parent->children_.front();
	parent->children_.pop_front();
	parent->nChildren_--;
      }else if(current != this){
	delete current;
	current = parent;
	parent = current->parent_;
      }
    }while(parent!=NULL || current->nChildren_ > 0);
  }

  /** 
   * Get a pointer to the parent node
   * \return pointer to parent, or NULL if this node has no parent
   */
  Node* getParent(){
    return parent_;
  }
  
  /**
   * Get pointers to the children
   * \param[out] children a vector containing pointers to children nodes
   * \return a pointers to a vector
   */ 
  void getChildren(std::vector<Node*> &children){
    children.clear();
    for( std::list<Node*>::iterator it = children_.begin(); it != children_.end(); it++ )
      children.push_back(*it);
  }

  /**
   * Get the number of children that this node has
   * \return number of children
   */
  int getChildrenCount(){
    return nChildren_;
  }

  /**
   * Add a child to this node. 
   * \note The parent node will take care of memory deallocation of children nodes
   * \param[in] child pointer to child node
   */
  void addChild(Node* child){
    child->parent_ = this;
    children_.push_front(child);
    nChildren_++;
  }

private:

  std::list<Node*> children_;
  int nChildren_;
  Node* parent_;
};

#endif
