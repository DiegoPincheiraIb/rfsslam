#ifndef GAUSSIAN_MIXTURE_HPP
#define GAUSSIAN_MIXTURE_HPP

#include <iostream>
#include "Landmark.hpp"
#include <vector>

/** 
 * \class GaussianMixture
 * \brief A class for mixture of Gaussians
 *
 * This class represents a mixture of weighted Gaussians using a binary tree.
 * It is designed to work in a Rao-Blackwellized particle filter
 *
 * \author Keith Leung
 *
 * \todo Function for checking Gaussians that are close for merging
 * \todo Function for merging Gaussians
 */

template< class Landmark >
class GaussianMixture
{
public:

  typedef Landmark* pLandmark;

  struct Gaussian{
    pLandmark landmark;
    double weight;
  };
  
  /** Default constructor */
  GaussianMixture();
  
  /** Desctructor */
  ~GaussianMixture();

  /** 
   * Add Gaussian 
   * \param p pointer to Landmark to add
   * \param w weight of the new Gaussian
   * \return number of Gaussians in the mixture
   */
  unsigned int addGaussian( pLandmark p, double w = 1);

  /** 
   * Remove Gaussian 
   * \param idx index number of Gaussian \ Landmark
   * \return number of Gaussians in the mixture
   */
  unsigned int removeGaussian( unsigned int idx );

  /**
   * Set the weight of a Gaussian indicated by the given index
   * \param idx index
   * \param w weight 
   */ 
  void setWeight( unsigned int idx, double w );

  /**
   * Get the state and weight of a Gaussian
   * \param idx index
   * \param p overwritten by pointer to landmark
   * \param w overwritten by the weight 
   */
  void getGaussian( unsigned int idx, pLandmark &p, double &w);

  /**
   * Prune the Gaussian mixture to remove Gaussians with weights that are
   * less than a threshold
   * \param t weight threshold
   * \return number of Gaussians removed
   */
  unsigned int prune( const double t );

protected:

  int n_; /**< number of Gaussians / Landmarks */
  std::vector<Gaussian> gList_; /**< container for Gaussians */ 
  std::vector<unsigned int> gListFreeIdx_; /**< Tracks free spots in gList_ */
  
};

template< class Landmark >
GaussianMixture<Landmark>::GaussianMixture(){
  n_ = 0;
}

template< class Landmark >
GaussianMixture<Landmark>::~GaussianMixture(){
  // Decrease reference count for all landmarks
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::addGaussian( pLandmark p, double w ){
  
  Gaussian g;
  g.landmark = p;
  g.landmark->incNRef();
  g.weight = w;
  
  if( gListFreeIdx_.empty() ){ 
    gList_.push_back(g);
  }else{
    unsigned int idx = gListFreeIdx_.back();
    gListFreeIdx_.pop_back();
    gList_[idx] = g;
  }
  
  n_++;
  return n_;
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::removeGaussian( unsigned int idx ){
  try{
    unsigned int nRef = gList_[idx].landmark->decNRef();
    if( nRef == 0 ){
      delete gList_[idx].landmark; // no more references to g, so delete it
    }
    if( idx == gList_.size() - 1){
      gList_.pop_back();
    }else{
      gList_[idx].landmark = NULL;
      gList_[idx].weight = 0;
      gListFreeIdx_.push_back(idx);
    }
    n_--;
    return n_;
  }catch(...){
    std::cout << "Error in removing Gaussian\n";
  }
}

template< class Landmark >
void GaussianMixture<Landmark>::setWeight( unsigned int idx, double w ){
  try{
    gList_[idx].weight = w;
  }catch(...){
    std::cout << "Unable to set Gaussian weight\n";
  }
}

template< class Landmark >
void GaussianMixture<Landmark>::getGaussian( unsigned int idx, pLandmark &p, double &w){
  try{
    p = gList_[idx].landmark;
    w = gList_[idx].weight;
  }catch(...){
    std::cout << "Unable to get Gaussian\n";
  }
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::prune( const double t ){
  
  unsigned int nPruned = 0;
  
  for( int idx = 0; idx < gList_.size(); idx++ ){
    if (gList_[idx].landmark != NULL && gList_[idx].weight < t){
      removeGaussian(idx);
      nPruned++;
    }
  }

  return nPruned;
}

#endif
