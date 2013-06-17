#ifndef GAUSSIAN_MIXTURE_HPP
#define GAUSSIAN_MIXTURE_HPP

#include <algorithm>
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

  struct Gaussian{
    pLandmark landmark;
    double weight;
  };
  
  /** Default constructor */
  GaussianMixture();
  
  /** Desctructor */
  ~GaussianMixture();

  /** 
   * Copy data from this Gaussian mixture to another.
   * Reference counts to all landmarks will be incremented by 1
   * \param other the other Gaussian mixture to which data is copied to
   */
  void copyTo( GaussianMixture *other);

  /** 
   * Add Gaussian 
   * \param p pointer to Landmark to add
   * \param w weight of the new Gaussian
   * \param allocateMem if true, assumes memory for Landmark has been allocated
   * and will not go out of scope or deleted by anything method outside this object.
   * If false, memory is allocated for a new Landmark and data from p is copied over.
   * \return number of Gaussians in the mixture
   */
  unsigned int addGaussian( pLandmark p, double w = 1, bool allocateMem = false);

  /** 
   * Remove a Gaussian and updates the number of Gaussians in the mixture 
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
   * Update an Gaussian in the mixture. 
   * If there is only 1 reference to the Landmark, it will be updated in place
   * If there are multiple references to the landmark, a new landmark will be
   * created, and the reference count to the original landmark will decrease by 1
   * \param idx index of Gaussian / landmark to update
   * \param p pointer to Landmark object with the updated data
   */
  void updateGaussian( unsigned int idx, pLandmark p);

  /**
   * Merge Gaussians that are within a certain Mahalanobis distance of each other
   * \param t Distance threshold
   * \param f_inflation Merged Gaussian inflation factor
   * \return number of merging operations
   */
  unsigned int merge(const double t = 0.1, const double f_inflation = 1.0);

  /**
   * Merge two Guassians if the second is within a Mahalanobis distance of the first. 
   * If merging occurs The resulting Gaussian will overwrite the first Gaussian.
   * The second one will be removed from the Gaussian mixture.
   * \param idx1 index of the first Gaussian
   * \param idx2 index of the second Gaussian
   * \param t distance threshold
   * \param f_inflation Merged Gaussian inflation factor
   * \return true if merging is successful
   */
  bool merge(unsigned int idx1, unsigned int idx2, 
	     const double t = 0.1, const double f_inflation = 1.0);

  /**
   * Prune the Gaussian mixture to remove Gaussians with weights that are
   * less than a threshold. 
   * \note Gaussians may have difference indices after using this function
   * \param t weight threshold
   * \return number of Gaussians removed
   */
  unsigned int prune( const double t );

protected:

  int n_; /**< number of Gaussians / Landmarks */
  std::vector<Gaussian> gList_; /**< container for Gaussians */ 

  /**					       
   * Sort the Gaussian mixture container from highest to lowest Gaussian weight 
   */
  void sortByWeight();
  
  /**
   * Comparison function for sorting Gaussian container
   */
  bool weightCompare(Gaussian a, Gaussian b);
  
};

////////// Implementation //////////

template< class Landmark >
GaussianMixture<Landmark>::GaussianMixture(){
  n_ = 0;
}

template< class Landmark >
GaussianMixture<Landmark>::~GaussianMixture(){

  for( int i = 0; i < gList_.size(); i++){
    removeGaussian(i);
  }

}

template< class Landmark >
void GaussianMixture<Landmark>::copyTo( GaussianMixture *other){
  other->n_ = n_;
  other->gList_ = gList_;
  for(int i = 0; i < gList_.size(); i++){
    gList_.landmark->incNRef();
  }
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::addGaussian( pLandmark p, double w, 
						     bool allocateMem){
  
  Gaussian g;

  if(allocateMem){
    typename Landmark::tState x;
    typename Landmark::tUncertainty S;
    p->get(x, S);
    p = new Landmark(x, S);
  }

  g.landmark = p;
  g.landmark->incNRef();
  g.weight = w;
  gList_.push_back(g);
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
    gList_[idx].landmark = NULL;
    gList_[idx].weight = 0;
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
void GaussianMixture<Landmark>::updateGaussian( unsigned int idx, pLandmark p){
  
  typename Landmark::tState x;
  typename Landmark::tUncertainty S;
  p->get(x, S);
  
  if (gList_[idx].landmark->getNRef() == 1){ // update in place
    gList_[idx].landmark->set(x, S);
  }else{
    gList_[idx].landmark->decNRef();
    gList_[idx].landmark = new Landmark(x, S);
  }
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::merge(const double t, 
					      const double f_inflation){
  unsigned int nMerged = 0;

  for(unsigned int i = 0; i < gList_.size(); i++){
    for(unsigned int j = i+1; j < gList_.size(); j++){
      if( merge(i, j, t, f_inflation) ){
	nMerged++;
      } 
    }
  }
  return nMerged;
}


template< class Landmark >
bool GaussianMixture<Landmark>::merge(unsigned int idx1, unsigned int idx2,
				      const double t, const double f_inflation){

  if (gList_[idx1].landmark == NULL ||
      gList_[idx2].landmark == NULL ){
    return false;
  }

  double w_m, w_1, w_2, d_mahalanobis; 
  Gaussian m;
  typename Landmark::tState x_1, x_2, x_m, d_12, d_1, d_2;
  typename Landmark::tUncertainty S_1, S_2, S_m;

  gList_[idx1].landmark->get(x_1, S_1);
  gList_[idx2].landmark->get(x_2, S_2);

  d_12 = x_1 - x_2;
  d_mahalanobis = d_12.transpose() * S_1.inverse() * d_12;
  if( d_mahalanobis > t ){
    return false;
  }
  d_mahalanobis = d_12.transpose() * S_2.inverse() * d_12;
  if( d_mahalanobis > t ){
    return false;
  }

  w_1 = gList_[idx1].weight;
  w_2 = gList_[idx2].weight;
  w_m = w_1  + w_2;

  x_m = (x_1 * w_1 + x_2 * w_2) / w_m;

  d_1 = x_m - x_1;
  d_2 = x_m - x_2;
  S_m = ( w_1 * ( S_1 + f_inflation * d_1 * d_1.transpose() ) +
          w_2 * ( S_2 + f_inflation * d_2 * d_2.transpose() ) ) / w_m;

  gList_[idx1].landmark->set(x_m, S_m);
  gList_[idx1].weight = w_m;

  gList_[idx1].landmark->set(x_m, S_m);
  gList_[idx1].weight = w_m;

  removeGaussian( idx2 ); // this also takes care of updating Gaussian count

  return true;
  
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::prune( const double t ){
  
  unsigned int nPruned = 0;
  if( gList_.size() <= 1 ){
    return nPruned;
  }

  sortByWeight(); // Sort from greatest to smallest weight
  
  // Binary search for the Gaussian with weight closest to and greater than t
  unsigned int min_idx = 0;
  unsigned int max_idx = gList_.size() - 1;
  unsigned int idx = (unsigned int)( (max_idx + min_idx) / 2 );
  unsigned int idx_old = idx + 1;
  double w = gList_[idx].weight;
  while(idx != idx_old){if( w > t ){
      max_idx = idx;
    }else if ( w <= t ){
      min_idx = idx;
    }
    idx_old = idx;
    idx = (unsigned int)( (max_idx + min_idx) / 2 );
    w = gList_[idx].weight;
  }
  while( w <= t && idx < gList_.size() ){
    idx++;
    w = gList_[idx].weight;
  }
  idx_old = idx; 
  while( idx < gList_.size() ){
    removeGaussian(idx); // this already takes care updating the Gaussian count
    idx++;
    nPruned++;
  }
  gList_.resize(idx_old); 

  return nPruned;
}

template< class Landmark >
void GaussianMixture<Landmark>::sortByWeight(){
  std::sort( gList_.begin(), gList_.end(), weightCompare );
}

template< class Landmark >
bool GaussianMixture<Landmark>::weightCompare(Gaussian a, Gaussian b){
  return a.weight > b.weight;
}

#endif
