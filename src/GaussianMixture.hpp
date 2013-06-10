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
 */

template< class Landmark >
class GaussianMixture
{
public:

  typedef Landmark* pLandmark;
  
  /** Default constructor */
  GaussianMixture();
  
  /** Desctructor */
  ~GaussianMixture();

  /** 
   * Add Gaussian 
   * \param g pointer to Landmark to add
   * \return number of Gaussians in the mixture
   */
  unsigned int addGaussian( pLandmark g);

  /** 
   * Remove Gaussian 
   * \param idx index number of Gaussian \ Landmark
   * \return number of Gaussians in the mixture
   */
  unsigned int removeGaussian( unsigned int idx );

protected:

  int n_; /**< number of Gaussians / Landmarks */
  std::vector<pLandmark> gList_; /**< contains pointers to Gaussians / landmarks */ 
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
unsigned int GaussianMixture<Landmark>::addGaussian( pLandmark g ){
  if( gListFreeIdx_.empty() ){ 
    gList_.push_back(g);
  }else{
    unsigned int idx = gListFreeIdx_.back();
    gListFreeIdx_.pop_back();
    gList_[idx] = g;
  }
  g->incNRef();
  n_++;
  return n_;
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::removeGaussian( unsigned int idx ){
  try{
    unsigned int nRef = gList_[idx]->decNRef();
    if( nRef == 0 ){
      delete gList_[idx]; // no more references to g, so delete it
    }
    if( idx == gList_.size() - 1){
      gList_.pop_back();
    }else{
      gList_[idx] = NULL;
      gListFreeIdx_.push_back(idx);
    }
    n_--;
    return n_;
  }catch(...){
    std::cout << "Error in removing Gaussian\n";
  }
}





#endif
