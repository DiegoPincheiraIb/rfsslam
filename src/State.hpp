// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP

#include <iostream>
#include <math.h>

double const PI = acos(-1);

/**
 * \class State
 * \brief An abstract base class for defining a state
 * \tparam StateType an Eigen vector;
 * \author Keith Leung
 *
 * \todo consider reading nDim_ using the size() function in Eigen
 */
template<class StateType>
class State
{

public:

  typedef StateType tState;

  /** Default constructor */
  State(){
    nDim_ = x_.size();
    if( nDim_ == 0 ){
      std::cerr << "Error: Dimension must be greater than 0\n";
      exit(-1);
    }
    x_ = StateType::Zero();
  }

  /** Constructor */
  State( StateType x ){ 
    nDim_ = x_.size();
    if( nDim_ == 0 ){
      std::cerr << "Error: Dimension must be greater than 0\n";
      exit(-1);
    }
    set(x); 
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
 * \tparam StateType An Eigen vector
 * \tparam UncertaintyType An Eigen matrix
 * \author Keith Leung
 */
template<class StateType, class UncertaintyType>
class StateWithUncertainty : public State<StateType>
{

public:

  typedef UncertaintyType tUncertainty;

  /** Default constructor */
  StateWithUncertainty(){
    if( Sx_.rows() != Sx_.cols() ){
      std::cerr << "Error: UncertaintyType must be a square matrix \n";
      exit(-1);
    }
    if( Sx_.rows() != this->x_.size() ){
      std::cerr << "Error: StateType and UncertaintyType dimension mismatch \n";
      exit(-1);
    }

    Sx_ = UncertaintyType::Zero();
  }

  /** Default destructor */
  ~StateWithUncertainty(){};

  /** 
   * Function for setting the pose state with uncertainty
   * \param x state to be set
   * \param Sx uncertainty to be set
   */
  void set( StateType x, UncertaintyType Sx){
    State<StateType>::set(x);
    Sx_ = Sx;
    SxInv_ = Sx_.inverse();
    Sx_det_ = Sx.determinant();
  }
  
  /** 
   * Function for getting the pose state with uncertianty
   * \param x state [overwritten]
   * \param Sx uncertainty [overwritten]
   */
  void get( StateType &x, UncertaintyType &Sx){
    State<StateType>::get(x);
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
  

private:

  UncertaintyType Sx_; /**< Covariance */
  UncertaintyType SxInv_; /**< Inverse covariance */
  double Sx_det_; /** Determinant of Sx_ */

};

#endif
