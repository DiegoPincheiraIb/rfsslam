// State class
// Keith Leung 2013

#ifndef STATE_HPP
#define STATE_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <iostream>
#include <math.h>
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
  RandomVec() : pSx_L_(NULL){

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
  RandomVec(VecType x, MatType Sx, double t = -1) : pSx_L_(NULL){
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
  RandomVec(VecType x, double t = -1) : pSx_L_(NULL){
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
    SxInv_ = Sx_.inverse();
    Sx_det_ = Sx.determinant();
    if( pSx_L_ != NULL){
      delete pSx_L_;
      pSx_L_ = NULL;
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

  /** 
   * Calculate the sqaured Mahalanobis distance 
   * from the current vector, scaled by the covariance
   * \param[in] x vector to which we measure the distance to
   * \return mahalanobis distance squared
   */
  double mahalanobisDist2( VecType &x){
    VecType e = x_ - x;
    return (e.transpose() * SxInv_ * e);
  }

  /**
   * Calculate the Mahalanobis distance from the current vector,
   * scaled by the covariance
   * \param[in] x vector to which we measure the distance to
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
   * Calculate the Mahalanobis distance from the current vector,
   * scaled by the covariance
   * \param[in] x vector to which we measure the distance to
   * \return mahalanobis distance
   */
  double mahalanobisDist( RandomVec<VecType, MatType> &x ){   
    VecType s;
    x.get(s);
    return mahalanobisDist( s );
  }

  /** 
   * Evaluate the Gaussian likelihood of a given evaluation point
   * the current vector as the mean, with its covariance
   * \param[in] x vector to the evaluation point
   * \return likelihood
   */
  double evalGaussianLikelihood( VecType &x ){
    double md2 = mahalanobisDist2( x );
    return ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim_) * Sx_det_ )  );
  }


  /** 
   * Evaluate the Gaussian likelihood of a given evaluation point
   * the current vector as the mean, with its covariance
   * \param[in] x vector to the evaluation point
   * \return likelihood
   */
  double evalGaussianLikelihood( RandomVec<VecType, MatType> &x ){ 
    VecType s;
    x.get(s);
    return evalGaussianLikelihood( s );
  }
  

private:

  Vec x_; /**< State */
  unsigned int nDim_; /**< Number of dimensions */
  MatType Sx_; /**< Covariance */
  MatType SxInv_; /**< Inverse covariance */
  double Sx_det_; /** Determinant of Sx_ */
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


template< class RandomVecDerivedClass >
class RandomVecMathTools
{
public:

  static void sample( RandomVecDerivedClass &s,
		      RandomVecDerivedClass &s_sample ){
    
    typename RandomVecDerivedClass::Vec x, indep_noise, e;
    typename RandomVecDerivedClass::Mat Sx, Sx_L;
    s.get(x, Sx);
    
    static boost::mt19937 rng_;
    static boost::normal_distribution<double> nd_;
    static boost::variate_generator< boost::mt19937, 
				     boost::normal_distribution<double> > 
      gen_(rng_, nd_);
    
    s.getCovCholeskyDecompLower(Sx_L);

    int n = Sx_L.cols();
    for(int i = 0; i < n; i++){
      indep_noise(i) = gen_();
    }
    e = Sx_L * indep_noise;
    x += e;
    s_sample.set( x, Sx );

  }

  static void sample( typename RandomVecDerivedClass::Vec &x,
		      typename RandomVecDerivedClass::Mat &Sx,
		      typename RandomVecDerivedClass::Mat &Sx_L,
		      RandomVecDerivedClass &s_sample ){
    
    typename RandomVecDerivedClass::Vec indep_noise, e;
    
    static boost::mt19937 rng_;
    static boost::normal_distribution<double> nd_;
    static boost::variate_generator< boost::mt19937, 
				     boost::normal_distribution<double> > 
      gen_(rng_, nd_);
    
    int n = Sx_L.cols();
    for(int i = 0; i < n; i++){
      indep_noise(i) = gen_();
    }
    e = Sx_L * indep_noise;
    x += e;
    s_sample.set( x, Sx );

  }

};


#endif
