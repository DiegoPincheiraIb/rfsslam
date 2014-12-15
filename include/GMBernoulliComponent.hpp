/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2014, Keith Leung, Felipe Inostroza
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
 
#ifndef GMMULTIBERNOULLI_HPP
#define GMMULTIBERNOULLI_HPP

#include "GaussianMixture.hpp"

namespace rfs{
/** 
 * \class GMMultiBernoulli
 * \brief A class representing a multi bernoulli random finite set, with a Gaussian mixture representing the spacial distribution of each Bernoulli component.
 *
 * This class represents a Multi Bernoulli Random finite set.
 * Each element has a Bernoulli Distribution according to GMBernoulliComponent.
 * 
 * \author Felipe Inostroza
 */
template< class Landmark > GMMultiBernoulli
{
public:
  
  typedef Landmark TLandmark;
  typedef Landmark* pLandmark;
  
  /*Defaulf Constructor*/
  GMMultiBernoulli();
  
  /*Copy Constructor*/
  GMMultiBernoulli(const GMMultiBernoulli& other);
  
  /*Destructor*/
  ~GMMultiBernoulli();
  
  /** 
   * Copy data from this GMMultiBernoulli  to another GMMultiBernoulli.
   * 
   * \param[out] other the other Multi Bernoulli Distribution to which data is copied to
   */
  void copyTo( GMMultiBernoulli *other); 
  


  std::vector<GMBernoulliComponent> tracks_;
  


}
////////// Implementation //////////
template< class Landmark >
GMMultiBernoulli<Landmark>::GMMultiBernoulli(){
  
}

template< class Landmark >
GMMultiBernoulli<Landmark>::GMMultiBernoulli(const GMMultiBernoulli& other){
  other.copyTo(this);
}

template< class Landmark >
GMMultiBernoulli<Landmark>::~GMMultiBernoulli(){
}

template< class Landmark >
void GMMultiBernoulli<Landmark>::copyTo(GMMultiBernoulli *other){

  other->tracks_.resize(this->tracks_.size());
  other->checkWeightingTableSize(WeightingTableNRows_ , WeightingTableNCols_)
  for(int i=0 ; i<this->tracks_.size() ; i++)
    tracks_[i].copyTo(other->tracks[i]);
}



/** 
 * \class GMBernoulliComponent
 * \brief A class representing a bernoulli random finite set, with a Gaussian mixture representing the spacial distribution.
 *
 * This class represents a Bernoulli Random finite set wich has no elements with probability (1-p) and a single element with probability p.
 * If the element exists it has a spacial distribution according to a gaussian mixture.
 *
 * \author Felipe Inostroza
 */

template< class Landmark > GMBernoulliComponent: GaussianMixture<Landmark>
{
public:
  
  typedef Landmark TLandmark;
  typedef Landmark* pLandmark;
  
  /*Defaulf Constructor*/
  GMBernoulliComponent();
  
  /*Copy Constructor*/
  GMBernoulliComponent(const GMBernoulliComponent& other);
  
  /*Destructor*/
  ~GMBernoulliComponent();
  
   /** 
   * Copy data from this GMBernoulliComponent  to another GMBernoulliComponent.
   * 
   * \param[out] other the other Gaussian mixture to which data is copied to
   */
  void copyTo( GMBernoulliComponent *other); 
  
  /**
   * Get probability of existance
   * \return Probability of existance of the Bernoulli Component
   **/
  double getP();
  /**
   * Set probability of existance
   * \param[in] p Probability of existance of the Bernoulli Component
   **/
  void setP(double p);
  
protected:



  double p_; /**< Probability of existence of the Bernoulli RFS. */
  double p_prev_; /**< Previous probability of existence of the Bernoulli RFS, used during particle weighting.*/
  
};
////////// Implementation //////////



template< class Landmark >
GMBernoulliComponent<Landmark>::GMBernoulliComponent(){
  p_=0;
  p_prev_=0;
}

template< class Landmark >
GMBernoulliComponent<Landmark>::GMBernoulliComponent(const GMBernoulliComponent& other){
  other.copyTo(this);
}

template< class Landmark >
GMBernoulliComponent<Landmark>::~GMBernoulliComponent(){
}

template< class Landmark >
void GMBernoulliComponent<Landmark>::copyTo(GMBernoulliComponent *other){
  other->p_=p_;
  other->p_prev_=p_prev_;
  GaussianMixture<Landmark>::copyTo(other);
}
template< class Landmark >
double GMBernoulliComponent<Landmark>::getP(){
  return p_;
}
template< class Landmark >
void GMBernoulliComponent<Landmark>::setP(double p){
  p_=p;
}
}


#endif
