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

#ifndef SPATIAL_INDEX_BOX_HPP
#define SPATIAL_INDEX_BOX_HPP

#include <assert.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>

#include "RandomVec.hpp"
#include "Tree.hpp"

namespace rfs{

  /**
   * \class Box
   * \tparam nDim number of dimensions
   * \tparam DataType type of object stored / associated within a box 
   * \brief Box used for quadtree / octree
   */
  template <unsigned int nDim = 3, class DataType = RandomVec<3> >
  class Box : public TreeNode< Box<nDim, DataType> >
  {
  public:

    typedef Eigen::Matrix<double, nDim, 1> Pos; /**< Location of a point */
    typedef boost::shared_ptr<DataType> DataPtr; /**< Data pointer */

    enum POSITION{
      POS_LOWER, POS_UPPER, POS_CENTER
    };
    
    /** 
     * Constructor 
     * \param[in] length size of box
     * \param[in] pos Array specifying the location of the minimum corner.
     */
    Box(double length = 1, 
	Eigen::Matrix<double, nDim, 1> pos = Eigen::Matrix<double, nDim, 1>::Zero() );

    /** Destructor */
    ~Box();
    
    /** 
     * Get the position of the box 
     * \param[in] p which part of the box 
     */
    Eigen::Matrix<double, nDim, 1> getPos(POSITION p = POS_LOWER);

    /**
     * Add data to box
     * \param[in] data Data to add
     */
    void addData(DataPtr data);

    /**
     * Get the number of data points stored in this box
     */
    unsigned int getDataSize();

    /**
     * Get a data point 
     * \param[in] idx Data point index.
     * \return pointer to Data, or NULL pointer if index is invalid
     */
    DataPtr getData(unsigned int idx);

    double const size_; /**< size of the box */

  protected:
    
    /**
     * Check if a point is within a box
     * \param[in] p query point
     * \return True or False
     */
    bool isInside(Pos p);

    Pos bound_min_; /**< point defining the intersection of all minimum bounds */
    Pos bound_max_; /**< point defining the intersection of all maximum bounds */

    std::vector< DataPtr > data_; /**< data storage */

  };

  /* Implementation */

  template <unsigned int nDim, class DataType>
  Box<nDim, DataType>::Box(double length, Eigen::Matrix<double, nDim, 1> pos) : size_(length){

    assert(nDim > 0);

    bound_min_ = pos;
    bound_max_ = pos + Pos::Ones() * length;
    
  }

  template <unsigned int nDim, class DataType>
  Box<nDim, DataType>::~Box(){}
 
  template <unsigned int nDim, class DataType>
  typename Box<nDim, DataType>::Pos Box<nDim, DataType>::getPos(Box<nDim, DataType>::POSITION p){
    
    switch(p){
    case POS_LOWER: return bound_min_;
    case POS_UPPER: return bound_max_;
    case POS_CENTER: return (bound_max_ - bound_min_) / 2 + bound_min_;
    }
    return bound_min_;
  }

  template <unsigned int nDim, class DataType>
  bool Box<nDim, DataType>::isInside(Box<nDim, DataType>::Pos p){
    
    if( (bound_max_ - p).minCoeff() >= 0 && (bound_min_ - p).maxCoeff() <= 0 ){
      return true;
    }else{
      return false;
    }
  }

  template <unsigned int nDim, class DataType>
  void Box<nDim, DataType>::addData( Box<nDim, DataType>::DataPtr data){
    data_.push_back(data);
  }

  template <unsigned int nDim, class DataType>
  unsigned int Box<nDim, DataType>::getDataSize(){
    return data_.size();
  }

  template <unsigned int nDim, class DataType>
  typename Box<nDim, DataType>::DataPtr Box<nDim, DataType>::getData(unsigned int idx){
    
    if(data_.size() > idx){
      return data_[idx];
    }else{
      return DataPtr();
    }
  }
  
}

#endif

