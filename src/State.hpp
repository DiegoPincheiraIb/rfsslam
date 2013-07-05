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

/**
 * \class RandomVecMathTools
 * \brief A collection of methods for use on a RandomVec derived class.
 * All methods are static, therefore this class does not need to be
 * instantiated.
 * \tparam RandomVec derived class
 */
template< class RandomVecDerivedClass >
class RandomVecMathTools
{
public:

  /** 
   * Calculate the sqaured Mahalanobis distance from a random vector,
   * \param[in] x_fm the random vector from which we can calcuating the 
   * distance, and of which the covariance is taken for scaling
   * \param[in] x_to the random vector to which we are measuring.
   * The covariance of this is ignored
   * \return mahalanobis distance squared
   */
  static double mahalanobisDist2( RandomVecDerivedClass &x_fm,
				  RandomVecDerivedClass &x_to){
    typename RandomVecDerivedClass::Vec x1;
    x_to.get(x1);
    return mahalanobisDist2( x_fm, x1 );
  }

  /** 
   * Calculate the sqaured Mahalanobis distance from a random vector,
   * \param[in] x_fm the random vector from which we can calcuating the 
   * distance, and of which the covariance is taken for scaling
   * \param[in] x_to the random vector to which we are measuring.
   * \return mahalanobis distance squared
   */
  static double mahalanobisDist2( RandomVecDerivedClass &x_fm,
				  typename RandomVecDerivedClass::Vec &x_to){
    typename RandomVecDerivedClass::Vec x0;
    typename RandomVecDerivedClass::Mat Sx0Inv;
    x_fm.get(x0);
    x_fm.getInvCov( Sx0Inv );
    typename RandomVecDerivedClass::Vec e = x_to - x0;
    return (e.transpose() * Sx0Inv * e);
  }

  /** 
   * Calculate the Mahalanobis distance from a random vector,
   * \param[in] x_fm the random vector from which we can calcuating the 
   * distance, and of which the covariance is taken for scaling
   * \param[in] x_to the random vector to which we are measuring.
   * The covariance of this is ignored
   * \return mahalanobis distance
   */
  static double mahalanobisDist( RandomVecDerivedClass &x_fm,
				 RandomVecDerivedClass &x_to){
    return sqrt( mahalanobisDist2( x_fm, x_to ) );
  }

  /** 
   * Calculate the Mahalanobis distance from a random vector,
   * \param[in] x_fm the random vector from which we can calcuating the 
   * distance, and of which the covariance is taken for scaling
   * \param[in] x_to the random vector to which we are measuring.
   * \return mahalanobis distance
   */
  static double mahalanobisDist( RandomVecDerivedClass &x_fm,
				  typename RandomVecDerivedClass::Vec &x_to){
    return sqrt( mahalanobisDist2( x_fm, x_to ) );
  }

  
  /** 
   * Evaluate the Gaussian likelihood of a evaluation point
   * \param[in] gaussian The Gaussian distribution represented by a random vector
   * \param[in] x_eval evaluation point
   * \return likelihood
   */
  static double evalGaussianLikelihood( RandomVecDerivedClass &gaussian,
				 RandomVecDerivedClass &x_eval ){
    double nDim = gaussian.getNDim();
    double covDet = gaussian.getCovDet();
    double md2 = mahalanobisDist2( gaussian, x_eval );
    return ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim) * covDet ) );
  }

  /** 
   * Evaluate the Gaussian likelihood of a evaluation point
   * \param[in] gaussian The Gaussian distribution represented by a random vector
   * \param[in] x_eval evaluation point
   * \return likelihood
   */
  static double evalGaussianLikelihood( RandomVecDerivedClass &gaussian,
					typename RandomVecDerivedClass::Vec &x_eval ){
    double nDim = gaussian.getNDim();
    double covDet = gaussian.getCovDet();
    double md2 = mahalanobisDist2( gaussian, x_eval );
    return ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim) * covDet ) );
  }  

  /** 
   * Sample the random vector
   * \param[in] s The random vector with the mean and covariance
   * \param[out] s_sample The sampled vector. The covariance and time
   * of the s are copied.
   */
  static void sample( RandomVecDerivedClass &s,
		      RandomVecDerivedClass &s_sample ){
    
    typename RandomVecDerivedClass::Vec x, indep_noise, e;
    typename RandomVecDerivedClass::Mat Sx, Sx_L;
    double t;
    s.get(x, Sx, t);
    
    s.getCovCholeskyDecompLower(Sx_L);

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
    s_sample.set( x, Sx, t );

  }

  /** 
   * Sample a random vector
   * \param[in] x The mean of the random vector
   * \param[in] Sx The covariance of the random vector
   * \param[in] Sx_L The lower Cholesky decomposition of the covariance
   * \param[in] t Time
   * \param[out] s_sample The sampled vector. The covariance and time
   * of the s are copied.
   */
  static void sample( typename RandomVecDerivedClass::Vec &x,
		      typename RandomVecDerivedClass::Mat &Sx,
		      typename RandomVecDerivedClass::Mat &Sx_L,
		      double &t,
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
    s_sample.set( x, Sx, t );

  }


};

#endif
