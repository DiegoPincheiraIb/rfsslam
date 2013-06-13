// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP

#include <math.h>

#define PI acos(-1)

/**
 * \class State
 * \brief An abstract base class for defining a state
 * \tparam State container
 * \author Keith Leung
 */
template<class StateType>
class State
{

public:

  typedef StateType tState;

  /** Default constructor */
  State(unsigned int nDim = 0)
    : nDim_(nDim){
    /** \todo Assert error here that nDim <= 0 */
  }

  /** Default destructor */
  ~State(){}

  /** 
   * Function for setting the pose state
   * \param x state to be set
   */
  void set( StateType x ){x_ = x;}
  
  /** 
   * Function for getting the pose state 
   * \param x state [overwritten]
   */
  void get( StateType &x ){x = x_;}

  /** 
   * Get the number of dimensions
   * \return number of dimensions
   */ 
  unsigned int getNDim(){ return nDim_; }

protected:

  StateType x_; /**< State */
  unsigned int nDim_; /**< Number of dimensions */

};

/**
 * \class StateWithUncertainty
 * \brief An abstract base class for defining the vehicle pose state
 * \tparam State container for the state
 * \author Keith Leung
 */
template<class StateType, class UncertaintyType>
class StateWithUncertainty : public State<StateType>
{

public:

  typedef UncertaintyType tUncertainty;

  /** Default constructor */
  StateWithUncertainty(unsigned int nDim = 0) : State<StateType>( nDim ){
    pSxInv_ = NULL;
  };

  /** Default destructor */
  ~StateWithUncertainty(){};

  /** 
   * Function for setting the pose state with uncertainty
   * \param x state to be set
   * \param Sx uncertainty to be set
   */
  void set( StateType x, UncertaintyType Sx){
    this->x_ = x;
    Sx_ = Sx;
    SxInv_ = Sx_.inverse();
    pSxInv_ = &SxInv_;
    Sx_det_ = Sx.determinant();
  }
  
  /** 
   * Function for getting the pose state with uncertianty
   * \param x state [overwritten]
   * \param Sx uncertainty [overwritten]
   */
  void get( StateType &x, UncertaintyType &Sx){
    x = this->x_;
    Sx = Sx_;
  }

  /** 
   * Abstract function for returning the sqaured Mahalanobis distance 
   * from this object's state
   * \param x the state to which we measure the distance to
   * \return mahalanobis distance squared
   */
  double mahalanobisDist2( StateType &x){
    
    StateType e = this->x_ - x;
    if(pSxInv_ == NULL){
      SxInv_ = Sx_.inverse();
      pSxInv_ = &SxInv_;
    }
    return (e.transpose() * SxInv_ * e);
  }

  /**
   * Function for returning the Mahalanobis distance from this object's state
   * \param x the state to which we measure the distance to
   * \return mahalanobis distance
   */
  double mahalanobisDist( StateType &x){
    double md2 = mahalanobisDist2( x );
    if( md2 >= 0)
      return sqrt( md2 );
    else
      return -1;
  }

  /**
   * Function for returning the Mahalanobis distance from this object's state
   * \param x object containing the state to which we measure the distance to
   * \return mahalanobis distance
   */
  double mahalanobisDist( State<StateType> &x ){   
    StateType s;
    x.get(s);
    return mahalanobisDist( s );
  }

  /** 
   * Evaluate the Gaussian likelihood of a state
   * \param x the state at which the likelihood will be evaluated
   * \return likelihood
   */
  double evalGaussianLikelihood( StateType x ){
    double md2 = mahalanobisDist2( x );
    return ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, this->nDim_) * Sx_det_ )  );
  }

  /** 
   * Evaluate the Gaussian likelihood of a state
   * \param x the state at which the likelihood will be evaluated
   * \return likelihood
   */
  double evalGaussianLikelihood( State<StateType> x ){ 
    StateType s;
    x.get(s);
    return evalGaussianLikelihood( s );
  }
  

protected:

  UncertaintyType Sx_; /**< Covariance */
  UncertaintyType SxInv_; /**< Inverse covariance */
  UncertaintyType *pSxInv_; /**< Pointer to inverse covariance */
  double Sx_det_; /** Determinant of Sx_ */

};

#endif
