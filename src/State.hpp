// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
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
    pSxInv_(NULL),
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
    pSxInv_(NULL),
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
    pSxInv_(NULL),
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
    if( pSxInv_ != NULL)
      delete pSxInv_;
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
    if( pSxInv_ != NULL){
      delete pSxInv_;
      pSxInv_ = NULL;
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
  void set( double t ){
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
  void getInvCov( MatType &Sx_inv){
    if(pSxInv_ == NULL){
      pSxInv_ = new MatType;
      *pSxInv_ = Sx_.inverse(); 
    }
    Sx_inv = *pSxInv_;
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
   * Get the dimension
   * \return dimension
   */ 
  unsigned int getNDim(){ return nDim_; }

private:

  Vec x_; /**< State */
  unsigned int nDim_; /**< Number of dimensions */
  MatType Sx_; /**< Covariance */
  MatType* pSxInv_; /**< Inverse covariance */
  double* pSx_det_; /** Determinant of Sx_ */
  MatType* pSx_L_; /** Lower triangular part of Cholesky decomposition on Sx_ */
  double t_; /**< time */

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
