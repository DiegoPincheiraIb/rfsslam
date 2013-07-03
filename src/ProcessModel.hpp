// Classes related to motion model
// Keith Leung 2013

#ifndef PROCESSMODEL_HPP
#define PROCESSMODEL_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "Measurement.hpp"
#include "Pose.hpp"

/**
 * \class ProcessModel
 * \brief An abstract class for defining vehicle motion models
 * \tparam StateType StateWithUncertainty derived type for state 
 * \tparam InputType Measurement derived for process input
 * \author Keith Leung
 */
template<class StateType, class InputType>
class ProcessModel
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef StateType TState;
  typedef InputType TInput;

  /** Default constructor */
  ProcessModel() : nd_(0, 1), gen_(rng_, nd_) {}

  /** Constructor 
   * \param S additive zero mean while Gaussian noise for this model 
   */
  ProcessModel(typename StateType::Mat &S) : S_zmgn_(S), nd_(0, 1), gen_(rng_, nd_) {
    if( S_zmgn_ != StateType::Mat::Zero() ){
      Eigen::LLT<typename StateType::Mat> cholesky( S );
      L_ = cholesky.matrixL();
    }
  }

  /** Default destructor */
  ~ProcessModel(){};

  /** Set the additive zero mean Gaussian noise covariance matrix 
   *  \param[in] S additive zero mean while Gaussian noise for this model
   */
  void setNoise(typename StateType::Mat &S){
    S_zmgn_ = S;
    Eigen::LLT<typename StateType::Mat> cholesky( S );
    L_ = cholesky.matrixL();
  }

  /** 
   * Abstract function for determining pose at time-step k from pose at 
   * time-step k-1
   * This must be implemented in a derived class
   * \param s_k pose at current time-step k. This can be the same object as s_km (to update in place). [overwritten]
   * \param s_km pose at previous time-step k-1
   * \param input_k input to process model
   * \param dT size of time-step
   */
  virtual void step( StateType &s_k, StateType &s_km, 
		     InputType &input_k , double const dT = 0 ) = 0;		     

  /**
   * Sample the process noise to predict the pose at k from k-1
   * \note This function can be overwritten in derived classes for implementing
   * other user-defined sampling methods. 
   * \warning This function does not check that the noise covariance matrices
   * are valid (i.e., semi-positive definite)
   * \param s_k pose at current time-step k. This can be the same object as s_km (to update in place). [overwritten]
   * \param s_km pose at previous time-step k-1
   * \param input_k input to process model
   * \param dT size of time-step
   * \param useAdditiveWhiteGaussianNoise if true, the output includes 
   * the zero-mean additive white Gaussian noise specified in the constructor
   * \param useInputWhiteGaussianNoise if true, the output includes
   * the noise specified in the input, and assumes that it is zero-mean white
   * Gaussian noise.
   */
  virtual void sample( StateType &s_k, StateType &s_km, 
		       InputType &input_k, double const dT = 0,
		       bool useAdditiveWhiteGaussianNoise = true,
		       bool useInputWhiteGaussianNoise = false ){
    
    if(useInputWhiteGaussianNoise){

      InputType in = input_k;
      typename InputType::Vec u;
      typename InputType::Mat Su, Su_L;
      double t;
      in.get( u, Su, t );
      Eigen::LLT<typename InputType::Mat> cholesky( Su );
      Su_L = cholesky.matrixL();
    
      int n = Su_L.cols();
      typename InputType::Vec randomVecNormal, randomVecGaussian;
      for(int i = 0; i < n; i++){
	randomVecNormal(i) = randn();
      }
      randomVecGaussian = Su_L * randomVecNormal;
      u = u + randomVecGaussian;
      step( s_k, s_km, in, dT );

    }else{
    
      step( s_k, s_km, input_k, dT );

    }
    
    if( useAdditiveWhiteGaussianNoise && S_zmgn_ != StateType::Mat::Zero() ){
    
      int n = S_zmgn_.cols();
      typename StateType::Vec randomVecNormal, randomVecGaussian, x_k;
      for(int i = 0; i < n; i++){
	randomVecNormal(i) = randn();
      }
      randomVecGaussian = L_ * randomVecNormal;

      s_k.State<typename StateType::Vec>::get(x_k);
      x_k = x_k + randomVecGaussian;
      s_k.State<typename StateType::Vec>::set(x_k);
    }
  }

protected:
  
  /** Covariance matrix for zero mean white Gaussian noise */
  typename StateType::Mat S_zmgn_;

  /** Lower triangular part of Cholesky decomposition on S_zmgn */
  typename StateType::Mat L_;

  /** Generate a random number from a normal distribution */
  double randn(){
    return gen_();
  }

private:

  boost::mt19937 rng_;
  boost::normal_distribution<double> nd_;
  boost::variate_generator< boost::mt19937, boost::normal_distribution<double> > gen_;

};

/**
 * \class StaticProcessModel
 * \brief A template process model class with not inputs
 * \author Keith Leung
 */
template< class StateType >
class StaticProcessModel : public ProcessModel< StateType, NullInput>
{

public:

  StaticProcessModel(){
    inputNoiseDefined = false;
  }

  StaticProcessModel(typename StateType::Mat &S): 
    ProcessModel< StateType, NullInput>(S){
    inputNoiseDefined = true;
  }

  ~StaticProcessModel(){}

  /** 
   * Determine the pose at time-step k from pose at time-step k-1
   * \param[out] s_k pose at current time-step k
   * \param[in] s_km pose at previous time-step k-1
   * \param[in] input_k input to process model
   * \param[in] dT size of time-step
   */
  void step( StateType &s_k, StateType &s_km, 
	     NullInput &input_k , double const dT = 1 ){
    
    if( inputNoiseDefined ){
      typename StateType::Vec x;
      typename StateType::Mat S;
      s_km.get(x, S);
      S += (this->S_zmgn_ * dT * dT);
      s_k.set(x, S);
    }else{
      s_k = s_km;
    }

  }
  
  /** 
   * Step function to allow for no inputs to the process model
   * \param[out] s_k State at current time-step k.
   * \param[in] s_km State at previous time-step k-1
   * \param[in] dT size of time-step
   */		     
  void staticStep( StateType &s_k, StateType &s_km, double const dT = 1){
    NullInput input;
    step(s_k , s_km , input , dT);
  }
		     
private:

  bool inputNoiseDefined;

};


////////// Example 2d Odometry Motion Model //////////

/**
 * \class OdometryMotionModel2d
 * \brief A 2d odometry motion model with translational and rotational
 *        displacement
 * \author Keith Leung
 */
class OdometryMotionModel2d : public ProcessModel< Pose2d, Odometry2d >
{
public:

  /** Default constructor */
  OdometryMotionModel2d();

  /** Constructor with process noise input 
   * \param S additive zero-mean white Gaussian noise covariance matrix
   */
  OdometryMotionModel2d( Pose2d::Mat S );

  /** Default destructor */
  ~OdometryMotionModel2d();

   /** 
   * This overrides the virtual function in the parent class for
   * determining the pose at time-step k from pose at time-step k-1
   * \param s_k pose at current time-step k [overwritten]
   * \param s_km pose at previous time-step k-1
   * \param input_k input to process model
   * \param dT size of time-step (not used)
   */
  void step( Pose2d &s_k, Pose2d &s_km, Odometry2d &input_k, 
	     double const dT=0);

private:

  Pose2d::Vec x_k_i_;       /**< \f[ \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k} \f]*/
  Pose2d::Vec x_km_i_;      /**< \f[ \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k-1} \f]*/
  Eigen::Vector2d p_k_i_;   /**< \f[ \begin{bmatrix} x \\ y \end{bmatrix}_{k} \f]*/
  Eigen::Vector2d p_km_i_;  /**< \f[ \begin{bmatrix} x \\ y \end{bmatrix}_{k-1} \f]*/
  double theta_k_;          /**< \f[ \theta_k \f\] */
  double theta_km_;         /**< \f[ \theta_k \f\] */
  Eigen::Matrix2d C_k_i_;   /**< rotation matrix k */
  Eigen::Matrix2d C_km_i_;  /**< rotation matrix k-1 */
  
  Eigen::Vector3d u_k_km_;  /**< odometry input */
  Eigen::Matrix3d Sd_k_km_; /**< uncertainty of translation input */
  Eigen::Vector2d dp_k_km_; /**< translation input */
  double dtheta_k_km_;      /**< rotation input */
  Eigen::Matrix2d C_k_km_;  /**< rotation matrix from odometry input */
  
};

#endif
