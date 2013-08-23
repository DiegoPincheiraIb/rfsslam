#ifndef GAUSSIAN_MIXTURE_HPP
#define GAUSSIAN_MIXTURE_HPP

#include <algorithm>
#include <iostream>
#include "Landmark.hpp"
#include "RandomVecMathTools.hpp"
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

  typedef Landmark TLandmark;
  typedef Landmark* pLandmark;

  /** \brief A data structure representing a weighted Gaussian distribution in GaussianMixture */ 
  struct Gaussian{
    pLandmark landmark; /**< pointer to landmark */
    double weight; /**< weight of Gaussian */
    double weight_prev; /**< previous weight of Gaussian (used by RB-PHD-Filter) */
  };
  
  /** Default constructor */
  GaussianMixture();
  
  /** Desctructor */
  ~GaussianMixture();

  /** 
   * Copy data from this Gaussian mixture to another.
   * Reference counts to all landmarks will be incremented by 1
   * \param[in/out] other the other Gaussian mixture to which data is copied to
   */
  void copyTo( GaussianMixture *other);

  /** 
   * Add Gaussian 
   * \param[in] p pointer to Landmark to add
   * \param[in] w weight of the new Gaussian
   * \param[in] allocateMem if true, assumes memory for Landmark has been allocated
   * and will not go out of scope or get deleted by anything method outside this object.
   * If false, memory is allocated for a new Landmark and data from p is copied over.
   * \return number of Gaussians in the mixture
   */
  unsigned int addGaussian( pLandmark p, double w = 1, bool allocateMem = false);

  /** 
   * Remove a Gaussian and updates the number of Gaussians in the mixture 
   * \param[in] idx index number of Gaussian \ Landmark
   * \return number of Gaussians in the mixture
   */
  unsigned int removeGaussian( unsigned int idx );

  /**
   * Get the number of Gaussians in the mixture
   * \return count 
   */
  unsigned int getGaussianCount();

  /**
   * Set the weight of a Gaussian indicated by the given index
   * \param[in] idx index
   * \param[in] w weight 
   */ 
  void setWeight( unsigned int idx, double w );

  /**
   * Get the weight of a Gaussian indicated by the given index
   * \param[in] idx index
   * \return weight 
   */ 
  double getWeight( unsigned int idx );

  /**
   * Get the state and weight of a Gaussian
   * \param[in] idx index
   * \return pointer to landmark. NULL if the Gaussian does not exist.
   */
  pLandmark getGaussian( unsigned int idx );
  
  /**
   * Get the state and weight of a Gaussian
   * \param[in] idx index
   * \param[out] p overwritten by pointer to landmark
   */
  void getGaussian( unsigned int idx, pLandmark &p);

  /**
   * Get the state and weight of a Gaussian
   * \param[in] idx index
   * \param[out] p overwritten by pointer to landmark
   * \param[out] w overwritten by the weight 
   */
  void getGaussian( unsigned int idx, pLandmark &p, double &w);

  /**
   * Get the state and weight of a Gaussian
   * \param[in] idx index
   * \param[out] p overwritten by pointer to landmark
   * \param[out] w overwritten by the weight 
   * \param[out] w_prev overwritten by the previous weight
   */
  void getGaussian( unsigned int idx, pLandmark &p, double &w, double &w_prev);

  /**
   * Update an Gaussian in the mixture. 
   * If there is only 1 reference to the Landmark, it will be updated in place
   * If there are multiple references to the landmark, a new landmark will be
   * created, and the reference count to the original landmark will decrease by 1
   * \param[in] idx index of Gaussian / landmark to update
   * \param[in] lm Landmark object with the updated data
   * \param[in] w weight of the updated Gaussian. No change if negative.
   * \return true if idx is valid and Gaussian is updated
   */
  bool updateGaussian( unsigned int idx, Landmark &lm, double w = -1);

  /**
   * Merge Gaussians that are within a certain Mahalanobis distance of each other
   * \param[in] t Distance threshold
   * \param[in] f_inflation Merged Gaussian inflation factor
   * \return number of merging operations
   */
  unsigned int merge(const double t = 0.31622776601, const double f_inflation = 1.0);

  /**
   * Merge two Guassians if the second is within a Mahalanobis distance of the first. 
   * If merging occurs The resulting Gaussian will overwrite the first Gaussian.
   * The second one will be removed from the Gaussian mixture.
   * \param[in] idx1 index of the first Gaussian
   * \param[in] idx2 index of the second Gaussian
   * \param[in] t distance threshold
   * \param[in] f_inflation Merged Gaussian inflation factor
   * \return true if merging is successful
   */
  bool merge(unsigned int idx1, unsigned int idx2, 
	     const double t = 0.1, const double f_inflation = 1.0);

  /**
   * Prune the Gaussian mixture to remove Gaussians with weights that are
   * less than a threshold. 
   * \note Gaussians may have difference indices after using this function
   * \param[in] t weight threshold
   * \return number of Gaussians removed
   */
  unsigned int prune( const double t );

  /**					       
   * Sort the Gaussian mixture container from highest to lowest Gaussian weight 
   */
  void sortByWeight();

protected:

  int n_; /**< number of Gaussians / Landmarks */
  std::vector<Gaussian> gList_; /**< container for Gaussians */ 
  bool isSorted_; /**< flag to prevent unecessary sorting */
  
  /**
   * Comparison function for sorting Gaussian container
   */
  static bool weightCompare(Gaussian a, Gaussian b);

   /** 
   * Add Gaussian and overwrite an existing spot in the Gaussian container
   * \param[in] idx element in gList_ to overwrite
   * \param[in] p pointer to Landmark to add
   * \param[in] w weight of the new Gaussian
   * \param[in] allocateMem if true, assumes memory for Landmark has been allocated
   * and will not go out of scope or get deleted by anything method outside this object.
   * If false, memory is allocated for a new Landmark and data from p is copied over.
   * \return number of Gaussians in the mixture
   */
  unsigned int addGaussian( unsigned int idx, pLandmark p, double w = 1, bool allocateMem = false);

  
};

////////// Implementation //////////

template< class Landmark >
GaussianMixture<Landmark>::GaussianMixture(){
  n_ = 0;
  isSorted_ = false;
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
    gList_[i].landmark->incNRef();
  }
  other->isSorted_ = isSorted_;
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::addGaussian( pLandmark p, double w, 
						     bool allocateMem){
  isSorted_ = false;

  if(allocateMem){
    Landmark* pNew = new Landmark;
    *pNew = *p;
    p = pNew;
  }

  Gaussian g = {p, w, 0};
  g.landmark->incNRef();
  gList_.push_back(g);
  n_++;
  return n_;
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::addGaussian( unsigned int idx, pLandmark p, double w, 
						     bool allocateMem){
  isSorted_ = false;

  if(allocateMem){
    Landmark* pNew = new Landmark;
    *pNew = *p;
    p = pNew;
  }

  if (gList_[idx].landmark != NULL){
    removeGaussian( idx );
  }

  gList_[idx].landmark = p;
  gList_[idx].landmark->incNRef();
  gList_[idx].weight = w;
  gList_[idx].weight_prev = 0;
  n_++;
  return n_;
}

template< class Landmark > 
unsigned int GaussianMixture<Landmark>::removeGaussian( unsigned int idx ){

  isSorted_ = false;
  try{
    if(gList_[idx].landmark != NULL){
      unsigned int nRef = gList_[idx].landmark->decNRef();
      if( nRef == 0 ){
	delete gList_[idx].landmark; // no more references to g, so delete it
      }
      gList_[idx].landmark = NULL;
      gList_[idx].weight = 0;
      gList_[idx].weight_prev = 0;
      n_--;
    }
    return n_;
  }catch(...){
    std::cout << "Error in removing Gaussian\n";
  }
}


template< class Landmark >
unsigned int GaussianMixture<Landmark>::getGaussianCount(){
  return n_;
}

template< class Landmark >
void GaussianMixture<Landmark>::setWeight( unsigned int idx, double w ){
  isSorted_ = false;
  try{
    gList_[idx].weight_prev = gList_[idx].weight;
    gList_[idx].weight = w;
  }catch(...){
    std::cout << "Unable to set Gaussian weight\n";
  }
}

template< class Landmark >
double GaussianMixture<Landmark>::getWeight( unsigned int idx){
  try{
    return gList_[idx].weight;
  }catch(...){
    std::cout << "Unable to get Gaussian weight\n";
  }
}

template< class Landmark >
typename GaussianMixture<Landmark>::pLandmark GaussianMixture<Landmark>::getGaussian( unsigned int idx ){
  return (gList_[idx].landmark);
}

template< class Landmark >
void GaussianMixture<Landmark>::getGaussian( unsigned int idx, pLandmark &p){
 
  try{
    p = gList_[idx].landmark;
  }catch(...){
    std::cout << "Unable to get Gaussian\n";
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
void GaussianMixture<Landmark>::getGaussian( unsigned int idx, pLandmark &p, double &w, double &w_prev){
 
  try{
    p = gList_[idx].landmark;
    w = gList_[idx].weight;
    w_prev = gList_[idx].weight_prev;
  }catch(...){
    std::cout << "Unable to get Gaussian\n";
  }
}

template< class Landmark >
bool GaussianMixture<Landmark>::updateGaussian( unsigned int idx, Landmark &lm, double w){
  
  isSorted_ = false;

  if( idx > gList_.size() ) 
    return false;

  if (gList_[idx].landmark == NULL)
    return false;

  typename Landmark::Vec x;
  typename Landmark::Mat S;
  lm.get(x, S);

  if ( w < 0 )
    w = getWeight( idx );

  if (gList_[idx].landmark->getNRef() == 1){ // update in place
    gList_[idx].landmark->set(x, S);
    gList_[idx].weight_prev = gList_[idx].weight;
    gList_[idx].weight = w;
  }else{
    double w_prev = gList_[idx].weight;
    removeGaussian( idx );
    Landmark* l = new Landmark;
    l->set(x, S);
    addGaussian( idx, l, w);
    gList_[idx].weight_prev = w_prev;
  }

  return true;

}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::merge(const double t, 
					      const double f_inflation){
  isSorted_ = false;
  unsigned int nMerged = 0;
  int nGaussians = gList_.size();

  for(unsigned int i = 0; i < nGaussians; i++){
   
    for(unsigned int j = i+1; j < nGaussians; j++){

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
  isSorted_ = false;

  if (gList_[idx1].landmark == NULL ||
      gList_[idx2].landmark == NULL ){
    return false;
  }

  double w_m, w_1, w_2; 
  double d_mahalanobis_1, d_mahalanobis_2; 

  w_1 = gList_[idx1].weight;
  w_2 = gList_[idx2].weight;
  
  double t2 = t * t;
  d_mahalanobis_1 = gList_[idx1].landmark->mahalanobisDist2( *(gList_[idx2].landmark) );
  if( d_mahalanobis_1 > t2 ){
    d_mahalanobis_2 = gList_[idx2].landmark->mahalanobisDist2( *(gList_[idx1].landmark) );
    if( d_mahalanobis_2 > t2 ){
      return false;
    }
  }

  w_m = w_1 + w_2;
  
  if( w_m == 0 )
    return false;

  typename Landmark::Vec x_1, x_2, x_m, d_1, d_2;
  typename Landmark::Mat S_1, S_2, S_m;

  gList_[idx1].landmark->get(x_1, S_1);
  gList_[idx2].landmark->get(x_2, S_2);

  x_m = (x_1 * w_1 + x_2 * w_2) / w_m;

  d_1 = x_m - x_1;
  d_2 = x_m - x_2;
  S_m = ( w_1 * ( S_1 + f_inflation * d_1 * d_1.transpose() ) +
          w_2 * ( S_2 + f_inflation * d_2 * d_2.transpose() ) ) / w_m;


  //std::cout << "\n" << x_1 << "\n" << S_1 << "\n\n";
  //std::cout << "\n" << x_2 << "\n" << S_2 << "\n\n";
  //std::cout << "\n" << x_m << "\n" << S_m << "\n\n";

  gList_[idx1].landmark->set(x_m, S_m);
  gList_[idx1].weight = w_m;
  gList_[idx1].weight_prev = 0;

  removeGaussian( idx2 ); // this also takes care of updating Gaussian count

  return true;
  
}

template< class Landmark >
unsigned int GaussianMixture<Landmark>::prune( const double t ){

  isSorted_ = false;

  unsigned int nPruned = 0;
  if( gList_.size() < 1 ){
    return nPruned;
  }

  sortByWeight(); // Sort from greatest to smallest weight
  
  // Binary search for the Gaussian with weight closest to and greater than t
  unsigned int min_idx = 0;
  unsigned int max_idx = gList_.size() - 1;
  unsigned int idx = (unsigned int)( (max_idx + min_idx) / 2 );
  unsigned int idx_old = idx + 1;
  double w = gList_[idx].weight;

  while(idx != idx_old){
    if( w <= t ){ 
      max_idx = idx;
    }else if ( w > t ){
      min_idx = idx;
    }
    idx_old = idx;
    idx = (unsigned int)( (max_idx + min_idx) / 2 );
    w = gList_[idx].weight;
  }
  while( w >= t && idx < gList_.size() ){
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
  if( !isSorted_ ){
    std::sort( gList_.begin(), gList_.end(), weightCompare );
    isSorted_ = true;
  }
}

template< class Landmark >
bool GaussianMixture<Landmark>::weightCompare(Gaussian a, Gaussian b){
  return a.weight > b.weight;
}

#endif
