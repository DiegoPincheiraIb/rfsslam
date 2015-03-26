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

#ifndef SPATIAL_INDEX_TREE_HPP
#define SPATIAL_INDEX_TREE_HPP

#include "SpatialIndexBox.hpp"
#include <iostream>

namespace rfs{

  /**
   * \class SpatialIndexTree
   * \tparam nDim number of dimensions
   * \tparam DataType type of object stored / associated with a node. This should be derived from the RandomVec class.
   * \brief Spatial indexing data structure (Octree when nDim = 3, Quadtree when nDim = 2)
   * \author Keith Leung
   */
  template <unsigned int nDim = 3, class DataType = RandomVec<3> >
  class SpatialIndexTree
  {
  public:

    typedef Box<nDim, DataType> TreeBox;
    typedef boost::shared_ptr<TreeBox> TreeBoxPtr;
    typedef Eigen::Matrix<double, nDim, 1> Pos;
    typedef boost::shared_ptr<DataType> DataPtr;

    /** 
     * Constructor
     * \param minBoxSize minimal geometric size of a box 
     */
    SpatialIndexTree(double minBoxSize = 1);
    
    /** Destructor */
    ~SpatialIndexTree();
    
    /**
     * Add data to the tree. This will use get() from DataType's base RandomVec class to obtain the position
     * from the first nDim elements.
     * \param[in] data Pointer to the data
     * \return box to which the data was added to
     */
    TreeBoxPtr addData(DataPtr data);

    /**
     * Remove data from the tree
     * \param[in] data Pointer to data to be removed
     */
    void removeData(DataPtr data);

    /**
     * Divide a box into smaller boxes
     * \param[in] b Pointer to box
     * \return True if successful
     */
    bool branch(TreeBoxPtr b);
    

  protected:

    /**
     * Search for the smaller box that contains a point
     * \param[in] p query point
     * \return pointer to box
     */
    TreeBoxPtr search(Pos p);

    TreeBoxPtr root_; /**< root of the tree */
    double const minBoxSize_; /**< smallest box size allowed */

  };

  /* Implementation */

  template<unsigned int nDim, class DataType>
  SpatialIndexTree<nDim, DataType>::SpatialIndexTree(double minBoxSize) : minBoxSize_(minBoxSize) {

    root_ = TreeBoxPtr( new TreeBox( minBoxSize*2, Pos::Ones() * -minBoxSize) );
    
    int nChildrenPerBox = pow(2,nDim);
    for(int i = 0; i < nChildrenPerBox; i++){
      
      Pos x = Pos::Zero();
      for(int d = 0; d < nDim; d++){
	x[d] = (i >> d) & 1 == 1; 
      }
      x *= minBoxSize;
      x += root_->getPos();

      root_->addChild( TreeBoxPtr(new TreeBox(minBoxSize, x)) );
     
    }
  }

  template<unsigned int nDim, class DataType>
  SpatialIndexTree<nDim, DataType>::~SpatialIndexTree(){}

  template<unsigned int nDim, class DataType>
  typename SpatialIndexTree<nDim, DataType>::TreeBoxPtr SpatialIndexTree<nDim, DataType>::addData(DataPtr data){
    
    // Data position
    Pos x = data->get().template head<nDim>();
    
    // Check if it fits inside root's bounding box
    if( root_->insideBox(x) ){ // within root's bounding box

      TreeBoxPtr b = search(x);
      if(b->size_ == minBoxSize_ || b->getDataSize() == 0 ){ // 
	b->addData(data);
      }else{ // check if adding data point will divide b
	Pos y = b->getData(0)->get().template head<nDim>(); // Get position of any arbitrary data point from b
	//check if the axis-aligned planes defining the center of the box also divide x and y
	Pos c = b->getPos( TreeBox::POS_CENTER );
	bool divide = false;
	for(unsigned int d = 0; d < nDim; d++){
	  if( (x[d] - c[d]) * (y[d] - c[d]) < 0 ){ // c divides x and y
	    divide = true;
	    break;
	  }
	} 
	if(!divide){
	  b->addData(data);
	}else{
	  
	}
      }
      
      
    }else{ // outside of root's bounding box
      
    }
    
  }

  template<unsigned int nDim, class DataType>
  void SpatialIndexTree<nDim, DataType>::removeData(SpatialIndexTree<nDim, DataType>::DataPtr data){
    
  }

  template<unsigned int nDim, class DataType>
  typename SpatialIndexTree<nDim, DataType>::TreeBoxPtr 
  SpatialIndexTree<nDim, DataType>::search(SpatialIndexTree<nDim, DataType>::Pos p){
    
    int nChildrenPerBox = pow(2,nDim);
    TreeBoxPtr b = root_; 
    unsigned int c = 0;
    
    while( b->getChildrenCount() != 0 ){
      for(int c = 0; c < nChildrenPerBox; c++){
	if (b->getChild()->insideBox(p)){
	  b = b->getChild();
	  break;
	}
      }
    }
    return b;

  }
  
  template<unsigned int nDim, class DataType>
  bool SpatialIndexTree<nDim, DataType>::branch(SpatialIndexTree<nDim, DataType>::TreeBoxPtr b){
    
    return true;
  }


}
#endif
