#ifndef RANDOM_VEC_MATH_TOOLS
#define RANDOM_VEC_MATH_TOOLS

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <math.h>

double const PI = acos(-1);

/**
 * \class RandomVecMathTools
 * \brief A collection of methods for use on a RandomVec derived class.
 * All methods are static, therefore this class does not need to be
 * instantiated.
 * \tparam RandomVec derived class
 * \author Keith Leung
 */
template< class RandomVecDerived >
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
  static double mahalanobisDist2( RandomVecDerived &x_fm,
				  RandomVecDerived &x_to){
    typename RandomVecDerived::Vec x1;
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
  static double mahalanobisDist2( RandomVecDerived &x_fm,
				  typename RandomVecDerived::Vec &x_to){
    typename RandomVecDerived::Vec x0;
    typename RandomVecDerived::Mat Sx0Inv;
    x_fm.get(x0);
    x_fm.getCovInv( Sx0Inv );
    typename RandomVecDerived::Vec e = x_to - x0;
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
  static double mahalanobisDist( RandomVecDerived &x_fm,
				 RandomVecDerived &x_to){
    return sqrt( mahalanobisDist2( x_fm, x_to ) );
  }

  /** 
   * Calculate the Mahalanobis distance from a random vector,
   * \param[in] x_fm the random vector from which we can calcuating the 
   * distance, and of which the covariance is taken for scaling
   * \param[in] x_to the random vector to which we are measuring.
   * \return mahalanobis distance
   */
  static double mahalanobisDist( RandomVecDerived &x_fm,
				  typename RandomVecDerived::Vec &x_to){
    return sqrt( mahalanobisDist2( x_fm, x_to ) );
  }

  
  /** 
   * Evaluate the Gaussian likelihood of a evaluation point
   * \param[in] gaussian The Gaussian distribution represented by a random vector
   * \param[in] x_eval evaluation point
   * \param[out] mDist2 if not NULL, stores the squared mahalanobis distance used 
   * to calculate the likelihood
   * \return likelihood
   */
  static double evalGaussianLikelihood( RandomVecDerived &gaussian,
					RandomVecDerived &x_eval,
					double* mDist2 = NULL){
    double nDim = gaussian.getNDim();
    double covDet = gaussian.getCovDet();
    double md2 = mahalanobisDist2( gaussian, x_eval );
    double l = ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim) * covDet ) );
    //If md2 is very large, l will become NAN;
    if( l != l)
      l = 0;
    if(mDist2 != NULL)
      *mDist2 = md2;
    return l;
  }

  /** 
   * Evaluate the Gaussian likelihood of a evaluation point
   * \param[in] gaussian The Gaussian distribution represented by a random vector
   * \param[in] x_eval evaluation point
   * \param[out] mDist2 if not NULL, stores the squared mahalanobis distance used 
   * to calculate the likelihood
   * \return likelihood
   */
  static double evalGaussianLikelihood( RandomVecDerived &gaussian,
					typename RandomVecDerived::Vec &x_eval,
					double* mDist2 = NULL ){
    double nDim = gaussian.getNDim();
    double covDet = gaussian.getCovDet();
    double md2 = mahalanobisDist2( gaussian, x_eval );
    double l = ( exp(-0.5 * md2 ) / sqrt( pow( 2*PI, nDim) * covDet ) );
    // If md2 is very large, l will become NAN;
    if( l != l)
      l = 0;
    if(mDist2 != NULL)
      *mDist2 = md2;
    return l;
  }  

  /** 
   * Sample the random vector
   * \param[in] s The random vector with the mean and covariance
   * \param[out] s_sample The sampled vector. The covariance and time
   * of the s are copied.
   */
  static void sample( RandomVecDerived &s,
		      RandomVecDerived &s_sample ){
    
    typename RandomVecDerived::Vec x, indep_noise, e;
    typename RandomVecDerived::Mat Sx, Sx_L;
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
  static void sample( typename RandomVecDerived::Vec &x,
		      typename RandomVecDerived::Mat &Sx,
		      typename RandomVecDerived::Mat &Sx_L,
		      double &t,
		      RandomVecDerived &s_sample ){
    
    typename RandomVecDerived::Vec indep_noise, e;


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
