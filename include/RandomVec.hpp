// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
//#include <Eigen/StdVector>
#include <iostream>
#include <stdio.h>

double const PI = acos(-1);

/**
 * \class RandomVec
 * \brief An abstract base class for deriving pose and measurement classes
 * \tparam VecType An Eigen vector of dimension n
 * \tparam MatType An Eigen matrix or dimension n x n
 * \author Keith Leung
 */
template<class VecType, class MatType>
class RandomVec
{

public:

  typedef VecType Vec;
  typedef MatType Mat;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** Default constructor */
  RandomVec() : 
    pSx_L_(NULL), 
    pSx_inv_(NULL),
    pSx_det_(NULL)
  {

    if ( !dimCheck() ){
      exit(-1);
    }
    x_.setZero();
    Sx_.setZero();
    t_ = -1;
  }

  /** 
   * Constructor 
   * \param x vector
   * \param Sx covariance
   * \param t time
   */
  RandomVec(VecType x, MatType Sx, double t = -1) : 
    pSx_L_(NULL),
    pSx_inv_(NULL),
    pSx_det_(NULL)
  {
    if ( !dimCheck() ){
      exit(-1);
    }
    set(x);
    setCov(Sx);
    t_ = t;
  }

  /** 
   * Constructor 
   * \param x vector
   * \param t time
   */
  RandomVec(VecType x, double t = -1) : 
    pSx_L_(NULL),
    pSx_inv_(NULL),
    pSx_det_(NULL){
    if ( !dimCheck() ){
      exit(-1);
    }
    set(x);
    Sx_.setZero();
    t_ = t;
  }

  /** Default destructor */
  ~RandomVec(){
    if( pSx_L_ != NULL)
      delete pSx_L_;
    if( pSx_inv_ != NULL)
      delete pSx_inv_;
    if( pSx_det_ != NULL)
      delete pSx_det_;
  };

  /** 
   * Set the vector
   * \param[in] x vector to be set
   */
  void set( VecType &x ){x_ = x;}

  /** 
   * Function for setting the pose uncertainty
   * \param[in] Sx uncertainty to be set
   */
  void setCov( MatType &Sx){
    Sx_ = Sx;
    if( pSx_L_ != NULL){
      delete pSx_L_;
      pSx_L_ = NULL;
    }
    if( pSx_inv_ != NULL){
      delete pSx_inv_;
      pSx_inv_ = NULL;
    }
    if( pSx_det_ != NULL){
      delete pSx_det_;
      pSx_det_ = NULL;
    }
  }

  /**
   * Set the time
   * \param[in] t time
   */
  void setTime( double t ){
    t_ = t;
  }

  /** 
   * Set the vector with a covariance matrix
   * \param[in] x vector to be set
   * \param[in] Sx covariance to be set
   */
  void set( VecType &x, MatType &Sx){
    set(x);
    setCov(Sx);
  }

  /** 
   * Set the vector with a time
   * \param[in] x vector to be set
   * \param[in] t time
   */
  void set( VecType &x, double t){
    set(x);
    t_ = t;
  }

  /** 
   * Set the vector with a covariance matrix, and time
   * \param[in] x vector to be set
   * \param[in] Sx covariance to be set
   * \param[in] t time
   */
  void set( VecType &x, MatType &Sx, double t){
    set(x);
    setCov(Sx);
    t_ = t;
  }


  /** 
   * Get the vector
   * \param[out] x vector
   */
  void get( VecType &x ){x = x_;}

  /** 
   * Getting the covariance matrix
   * \param[out] Sx uncertainty 
   */
  void getCov( MatType &Sx){
    Sx = Sx_;
  }

  /**
   * Get the lower triangular part of the Cholesky decomposition 
   * on the covariance matrx Sx_ 
   * \param[out] Sx_Chol_L The lower triangular part of the Choloesky decomposition
   */
  void getCovCholeskyDecompLower( MatType &Sx_Chol_L){
    if(pSx_L_ == NULL){
      pSx_L_ = new MatType;
      Eigen::LLT<MatType> cholesky( Sx_ );
      *pSx_L_ = cholesky.matrixL();
    }
    Sx_Chol_L = *pSx_L_;
  }

  /** 
   * Get the invserve covariance matrix 
   * \param[out] Sx_inv inverse covariance
   */ 
  void getCovInv( MatType &Sx_inv){
    if(pSx_inv_ == NULL){
      pSx_inv_ = new MatType;
      *pSx_inv_ = Sx_.inverse(); 
    }
    Sx_inv = *pSx_inv_;
  }

  /** 
   * Get the determinant of the covariance
   * \return determinant
   */
  double getCovDet(){
    if(pSx_det_ == NULL){
      pSx_det_ = new double;
      *pSx_det_ = Sx_.determinant();
    }
    return *pSx_det_;
  }

  /** 
   * Get the vector and covariance matrix
   * \param[out] x vector
   * \param[out] Sx uncertainty
   */
  void get( VecType &x, MatType &Sx){
    get(x);
    getCov(Sx);
  }

  /** 
   * Get the vector and time
   * \param[out] x vector
   * \param[out] t time
   */
  void get( VecType &x, double &t){
    get(x);
    t = t_;
  }

  /** 
   * Get the vector, covariance matrix, and time
   * \param[out] x vector
   * \param[out] Sx uncertainty
   * \param[out] t time
   */
  void get( VecType &x, MatType &Sx, double &t){
    get(x);
    getCov(Sx);
    t = t_;
  }

  /** 
   * Get an element of the vector
   * \param[in] n element index
   * \return element n
   */
  double get( int n ){ return x_(n);}

  /**
   * Get the time
   * \return time
   */
  double getTime(){
    return t_;
  }


  /** 
   * Get the dimension
   * \return dimension
   */ 
  unsigned int getNDim(){ return nDim_; }

  /**
   * Calculate the squared Mahalanobis distance to another random vector of the same type
   */
  double mahalanobisDist2( RandomVec<VecType, MatType> &to ){
    if(pSx_inv_ == NULL){
      pSx_inv_ = new MatType;
      *pSx_inv_ = Sx_.inverse(); 
    }
    e_ = to.x_ - x_;
    return (e_.transpose() * *pSx_inv_ * e_);
  }

  /**
   * Calculate the squared Mahalanobis distance to another random vector of the same type
   */
  double mahalanobisDist2( typename RandomVec<VecType, MatType>::Vec &to_x ){
    if(pSx_inv_ == NULL){
      pSx_inv_ = new MatType;
      *pSx_inv_ = Sx_.inverse(); 
    }
    e_ = to_x - x_;
    return (e_.transpose() * *pSx_inv_ * e_);
  }

  /**
   * Calculate likelihood
   */ 
  double evalGaussianLikelihood( RandomVec<VecType, MatType> &x_eval,
				 double* mDist2 = NULL){
    if(pSx_det_ == NULL){
      pSx_det_ = new double;
      *pSx_det_ = Sx_.determinant();
    }
    double md2 = mahalanobisDist2( x_eval );
    double l = ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim_ ) * *pSx_det_ ) );
    if( l != l) //If md2 is very large, l will become NAN;
      l = 0;
    if(mDist2 != NULL)
      *mDist2 = md2;
    return l;
  }

  /**
   * Calculate likelihood
   */ 
  double evalGaussianLikelihood( typename RandomVec<VecType, MatType>::Vec &x_eval,
				 double* mDist2 = NULL){
    if(pSx_det_ == NULL){
      pSx_det_ = new double;
      *pSx_det_ = Sx_.determinant();
    }
    double md2 = mahalanobisDist2( x_eval );
    double l = ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim_ ) * *pSx_det_ ) );
    if( l != l) //If md2 is very large, l will become NAN;
      l = 0;
    if(mDist2 != NULL)
      *mDist2 = md2;
    return l;
  }

private:

  VecType x_; /**< State */
  unsigned int nDim_; /**< Number of dimensions */
  MatType Sx_; /**< Covariance */
  MatType* pSx_inv_; /**< Inverse covariance */
  double* pSx_det_; /** Determinant of Sx_ */
  MatType* pSx_L_; /** Lower triangular part of Cholesky decomposition on Sx_ */
  double t_; /**< time */

  VecType e_; /**< temporary */

  /** Dimensionality check during initialization */
  bool dimCheck(){

    if( Sx_.rows() != Sx_.cols() ){
      std::cerr << "Error: MatType must be a square matrix \n";
      return false;
    }
    if( Sx_.rows() != x_.size() ){
      std::cerr << "Error: VecType and MatType dimension mismatch \n";
      return false;
    }
    nDim_ = x_.size();
    return true;
  }

};

#endif
