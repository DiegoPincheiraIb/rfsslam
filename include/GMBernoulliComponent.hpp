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

#include <vector>
#include "GaussianMixture.hpp"

namespace rfs {

/** 
 * \class GMBernoulliComponent
 * \brief A class representing a bernoulli random finite set, with a Gaussian mixture representing the spacial distribution.
 *
 * This class represents a Bernoulli Random finite set wich has no elements with probability (1-p) and a single element with probability p.
 * If the element exists it has a spacial distribution according to a gaussian mixture.
 *
 * \author Felipe Inostroza
 */

template<class Landmark>
  class GMBernoulliComponent : public GaussianMixture<Landmark>
  {
  public:

    typedef Landmark TLandmark;
    typedef Landmark* pLandmark;

    /*Default Constructor*/
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
    void copyTo(GMBernoulliComponent *other) const;

    /**
     * Assignment operator
     * \param[in] rhs the right-hand-side from which data is copied
     */
    GMBernoulliComponent& operator=( const GMBernoulliComponent& rhs );

    /**
     * Get probability of existence
     * \return Probability of existence of the Bernoulli Component
     **/
    double getP();
    /**
     * Set probability of existence
     * \param[in] p Probability of existence of the Bernoulli Component
     **/
    void setP(double p);
    /**
     * Getter for the maximum weight gaussian in the spatial distribution
     * @param[out] landmark Pointer to the maximum weight gaussian in the gaussian mixture.
     */
    void getMaxGaussian(pLandmark &landmark);


    double getPrevP() const;
    void setPrevP(double prevP);


  protected:

    double p_; /**< Probability of existence of the Bernoulli RFS. */
    double prevP_; /**< Previous probability of existence of the Bernoulli RFS, used during particle weighting.*/

  };
////////// Implementation //////////

template<class Landmark>
  GMBernoulliComponent<Landmark>::GMBernoulliComponent() {
    p_ = 0;
    prevP_ = 0;
  }

template<class Landmark>
  GMBernoulliComponent<Landmark>::GMBernoulliComponent(const GMBernoulliComponent& other) {
    other.copyTo(this);
  }

template<class Landmark>
  GMBernoulliComponent<Landmark>::~GMBernoulliComponent() {
  }

template<class Landmark>
  void GMBernoulliComponent<Landmark>::copyTo(GMBernoulliComponent *other) const {
    other->p_ = p_;
    other->prevP_ = prevP_;
    GaussianMixture<Landmark>::copyTo(other);
  }
template<class Landmark>
GMBernoulliComponent<Landmark>& GMBernoulliComponent<Landmark>::operator=(const GMBernoulliComponent &other) {
  other.copyTo(this);
  return *this;
}


template<class Landmark>
  double GMBernoulliComponent<Landmark>::getP() {
    return p_;
  }
template<class Landmark>
  void GMBernoulliComponent<Landmark>::setP(double p) {
    p_ = p;
  }


template<class Landmark>
  double rfs::GMBernoulliComponent<Landmark>::getPrevP() const {
    return prevP_;
  }
template<class Landmark>
  void rfs::GMBernoulliComponent<Landmark>::setPrevP(double prevP) {
    prevP_ = prevP;
  }

template<class Landmark>
  void rfs::GMBernoulliComponent<Landmark>::getMaxGaussian(pLandmark& landmark) {
    double maxW = 0;
    int maxWidx = 0;
    for (int i = 0; i < this->n_; i++) {
      if (this->gList_[i].landmark != NULL && this->gList_[i].weight > maxW) {
        maxW = this->gList_[i].weight;
        maxWidx = i;
      }
    }
    landmark = this->gList_[maxWidx].landmark;
    return;

  }

/**
 * \class GMMultiBernoulli
 * \brief A class representing a multi bernoulli random finite set, with a Gaussian mixture representing the spacial distribution of each Bernoulli component.
 *
 * This class represents a Multi Bernoulli Random finite set.
 * Each element has a Bernoulli Distribution according to GMBernoulliComponent.
 *
 * \author Felipe Inostroza
 */
template<class Landmark>
  class GMMultiBernoulli
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
    void copyTo(GMMultiBernoulli *other) const;

    /**
     * Prune the Gaussian mixtures to remove Gaussians with weights that are
     * less than a given threshold.
     * \note Gaussians may have difference indices after using this function due to sorting performed on gList_
     * \param[in] t weight threshold, below which Gaussians are removed.
     */
    void pruneGM(const double t);

    /**
     * Check all Gaussians in the mixtures and merge those that are within a certain Mahalanobis distance of each other
     * \param[in] t distance threshold (the default value squares to 0.1)
     * \param[in] f_inflation merged Gaussian covariance inflation factor (default value causes no inflation)
     *
     */
    void mergeGM(const double t = 0.31622776601, const double f_inflation = 1.0);

    /**
     * Prune tracks with probability of existence lower than the threshold
     * \note Track order may change
     * @param t probability of existence threshold
     */

    void pruneTracks(const double t);
    /**
     * sort tracks according to probability of existence
     */
    void sortByPE();
    /**
     * Compare method used for sorting tracks
     * @param a track
     * @param b track
     * @return true if track a has greater probability of existence than track b
     */
    static bool probExistenceCompare(GMBernoulliComponent<TLandmark> a, GMBernoulliComponent<TLandmark> b);
    std::vector< GMBernoulliComponent< TLandmark > > tracks_;

  };
////////// Implementation //////////
template<class Landmark>
  GMMultiBernoulli<Landmark>::GMMultiBernoulli() {

  }
template<class Landmark>
  void GMMultiBernoulli<Landmark>::sortByPE() {
    std::sort(tracks_.begin(), tracks_.end(), probExistenceCompare);
  }

template<class Landmark>
  GMMultiBernoulli<Landmark>::GMMultiBernoulli(const GMMultiBernoulli& other) {
    other.copyTo(this);
  }

template<class Landmark>
  GMMultiBernoulli<Landmark>::~GMMultiBernoulli() {
  }

template<class Landmark>
  void GMMultiBernoulli<Landmark>::copyTo(GMMultiBernoulli *other) const {

    other->tracks_.resize(this->tracks_.size());
    for (int i = 0; i < this->tracks_.size(); i++)
      tracks_[i].copyTo(&(other->tracks_[i]));
  }

template<class Landmark>
  void GMMultiBernoulli<Landmark>::pruneGM(const double t) {
    for (int i = 0; i < tracks_.size(); i++) {
      tracks_[i].prune(t);
    }

  }
template<class Landmark>
  void GMMultiBernoulli<Landmark>::mergeGM(const double t, const double f_inflation) {
    for (int i = 0; i < tracks_.size(); i++) {
      tracks_[i].merge(t, f_inflation);
    }
  }
template<class Landmark>
  void GMMultiBernoulli<Landmark>::pruneTracks(const double t) {
    int pruned = 0;
    for (int i = tracks_.size() - 1; i >= 0; i--) {
      if (tracks_[i].getP() < t) {
        tracks_[tracks_.size() - 1 - pruned].copyTo(&tracks_[i]);
        pruned++;
      }
    }
    tracks_.resize(tracks_.size() - pruned);

  }
template< class Landmark >
bool GMMultiBernoulli<Landmark>::probExistenceCompare(GMBernoulliComponent<TLandmark> a, GMBernoulliComponent<TLandmark> b){
  return a.getP() > b.getP();
}
} //namespace

#endif
