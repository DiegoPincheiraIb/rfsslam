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
 * \tparam PoseType Pose derived type for state 
 * \tparam InputType Measurement derived for process input
 * \author Keith Leung
 *
 * \todo Add function for predicting covariance for step function
 */
template<class PoseType, class InputType>
class ProcessModel
{
public:

  typedef PoseType TPose;
  typedef InputType TInput;

  /** Default constructor */
  ProcessModel() : nd_(0, 1), gen_(rng_, nd_) {}

  /** Constructor 
   * \param S additive zero mean while Gaussian noise for this model 
   */
  ProcessModel(typename PoseType::Mat &S) : S_zmgn_(S), nd_(0, 1), gen_(rng_, nd_) {
    if( S_zmgn_ != PoseType::Mat::Zero() ){
      Eigen::LLT<typename PoseType::Mat> cholesky( S );
      L_ = cholesky.matrixL();
    }
  }

  /** Default destructor */
  ~ProcessModel(){};

  /** 
   * Abstract function for determining pose at time-step k from pose at 
   * time-step k-1
   * This must be implemented in a derived class
   * \param s_k pose at current time-step k [overwritten]
   * \param s_km pose at previous time-step k-1
   * \param input_k input to process model
   * \param dT size of time-step
   */
  virtual void step( PoseType &s_k, PoseType &s_km, 
		     InputType &input_k, double const dT = 0 ) = 0;

  /**
   * Sample the zero mean white Gaussian process noise
   * to predict the pose at k from k-1
   * \note For other noise models, such as noise on the input which
   * requires more complex noise propagation methods, override this 
   * function in a derived class
   * \param s_k pose at current time-step k [overwritten]
   * \param s_km pose at previous time-step k-1
   * \param input_k input to process model
   * \param dT size of time-step
   */
  virtual void sample( PoseType &s_k, PoseType &s_km, 
		       InputType &input_k, double const dT = 0 ){
    
    step( s_k, s_km, input_k, dT );

    if( S_zmgn_ != PoseType::Mat::Zero() ){
    
      boost::mt19937 rng;
      boost::normal_distribution<double> nd(0, 1);
      boost::variate_generator< boost::mt19937, boost::normal_distribution<double> > gen(rng, nd);
      int n = S_zmgn_.cols();
      typename PoseType::Vec randomVecUniform, randomVecGaussian, x_k;
      for(int i = 0; i < n; i++){
	randomVecUniform(i) = randn();
      }
      randomVecGaussian = L_ * randomVecUniform;

      s_k.get(x_k);
      x_k = x_k + randomVecGaussian;
      s_k.set(x_k);
    }
  }

protected:
  
  /** Covariance matrix for zero mean white Gaussian noise */
  typename PoseType::Mat S_zmgn_;

  /** Lower triangular part of Cholesky decomposition on S_zmgn */
  typename PoseType::Mat L_;

  /** Generate a random number from a normal distribution */
  double randn(){
    return gen_();
  }

private:

  boost::mt19937 rng_;
  boost::normal_distribution<double> nd_;
  boost::variate_generator< boost::mt19937, boost::normal_distribution<double> > gen_;

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
