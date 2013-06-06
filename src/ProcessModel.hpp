// Classes related to motion model
// Keith Leung 2013

#ifndef PROCESSMODEL_HPP
#define PROCESSMODEL_HPP

#include <Eigen/Core>
#include "Measurement.hpp"
#include "Pose.hpp"

/**
 * \class ProcessModel
 * \brief An abstract class for defining vehicle motion models
 * \tparam StateType Object type for state 
 * \tparam InputType Object type for process input information
 * \author Keith Leung
 *
 * \todo Add function for predicting covariance for step function
 */
template<class StateType, class InputType>
class ProcessModel
{
public:

  /** Default constructor */
  ProcessModel(){};

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
  virtual void step( StateType &s_k, StateType &s_km, 
		     InputType &input_k, double const dT = 0 ) = 0;

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

  Pose2d::vec x_k_i_;       /**< \f[ \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k} \f]*/
  Pose2d::vec x_km_i_;      /**< \f[ \begin{bmatrix} x \\ y \\ \theta \end{bmatrix}_{k-1} \f]*/
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
