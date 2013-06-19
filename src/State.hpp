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
 * \tparam VecType an Eigen vector;
 * \author Keith Leung
 */
template<class VecType>
class State
{

public:

  typedef VecType Vec;

  /** Default constructor */
  State(){
    nDim_ = x_.size();
    if( nDim_ == 0 ){
      std::cerr << "Error: Dimension must be greater than 0\n";
      exit(-1);
    }
    x_ = Vec::Zero();
  }

  /** Constructor */
  State( Vec &x ){ 
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
  void set( Vec &x ){x_ = x;}

  /** 
   * Function for setting the pose state, exact same as set
   * \param x state to be set
   */
  void setState( Vec &x ){x_ = x;}
  
  /** 
   * Function for getting the pose state 
   * \param x state [overwritten]
   */
  void get( Vec &x ){x = x_;}

  /** 
   * Function for getting the pose state, exact same as get 
   * \param x state [overwritten]
   */
  void getState( Vec &x ){x = x_;}

  /** 
   * Get the number of dimensions
   * \return number of dimensions
   */ 
  unsigned int getNDim(){ return nDim_; }

protected:

  Vec x_; /**< State */
  unsigned int nDim_; /**< Number of dimensions */

};

/**
 * \class StateWithUncertainty
 * \brief An abstract base class for defining the vehicle pose state
 * \tparam VecType An Eigen vector
 * \tparam MatType An Eigen matrix
 * \author Keith Leung
 */
template<class VecType, class MatType>
class StateWithUncertainty : public State<VecType>
{

public:

  typedef MatType Mat;

  /** Default constructor */
  StateWithUncertainty(){
    if( Sx_.rows() != Sx_.cols() ){
      std::cerr << "Error: MatType must be a square matrix \n";
      exit(-1);
    }
    if( Sx_.rows() != this->x_.size() ){
      std::cerr << "Error: VecType and MatType dimension mismatch \n";
      exit(-1);
    }

    Sx_ = MatType::Zero();
  }

  /** Default destructor */
  ~StateWithUncertainty(){};

  /** 
   * Function for setting the pose state with uncertainty
   * \param x state to be set
   * \param Sx uncertainty to be set
   */
  void set( VecType &x, MatType &Sx){
    State<VecType>::set(x);
    setCov(Sx);
  }
  
  /** 
   * Function for getting the pose state with uncertianty
   * \param x state [overwritten]
   * \param Sx uncertainty [overwritten]
   */
  void get( VecType &x, MatType &Sx){
    State<VecType>::get(x);
    getCov(Sx);
  }

  /** 
   * Function for setting the pose uncertainty
   * \param Sx uncertainty to be set
   */
  void setCov( MatType &Sx){
    Sx_ = Sx;
    SxInv_ = Sx_.inverse();
    Sx_det_ = Sx.determinant();
  }
  
  /** 
   * Function for getting the pose uncertianty
   * \param Sx uncertainty [overwritten]
   */
  void getCov( MatType &Sx){
    Sx = Sx_;
  }

  /** 
   * Abstract function for returning the sqaured Mahalanobis distance 
   * from this object's state
   * \param x the state to which we measure the distance to
   * \return mahalanobis distance squared
   */
  double mahalanobisDist2( VecType &x){
    
    VecType e = this->x_ - x;
    return (e.transpose() * SxInv_ * e);
  }

  /**
   * Function for returning the Mahalanobis distance from this object's state
   * \param x the state to which we measure the distance to
   * \return mahalanobis distance
   */
  double mahalanobisDist( VecType &x){
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
  double mahalanobisDist( State<VecType> &x ){   
    VecType s;
    x.get(s);
    return mahalanobisDist( s );
  }

  /** 
   * Evaluate the Gaussian likelihood of a state
   * \param x the state at which the likelihood will be evaluated
   * \return likelihood
   */
  double evalGaussianLikelihood( VecType &x ){
    double md2 = mahalanobisDist2( x );
    return ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, this->nDim_) * Sx_det_ )  );
  }

  /** 
   * Evaluate the Gaussian likelihood of a state
   * \param x the state at which the likelihood will be evaluated
   * \return likelihood
   */
  double evalGaussianLikelihood( State<VecType> &x ){ 
    VecType s;
    x.get(s);
    return evalGaussianLikelihood( s );
  }
  

private:

  MatType Sx_; /**< Covariance */
  MatType SxInv_; /**< Inverse covariance */
  double Sx_det_; /** Determinant of Sx_ */

};

#endif
