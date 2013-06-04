// Classes related to motion model
// Keith Leung 2013

#ifndef PROCESSMODEL_HPP
#define PROCESSMODEL_HPP

#include <Eigen/Core>
#include "Pose.hpp"

/** 
 * \class Measurement
 * \brief An abstract class for the inputs to a process model
 * \tparam InputValType Object type representing process inputs
 * \tparam InputUncertaintyType Object type representing input uncertainties
 * \author Keith Leung
 */
template<class InputValType, class InputUncertaintyType>
class Measurement
{
public:

  /** Default constructor */
  Measurement(){ t_ = 0; };

  /** Default destructor */
  ~Measurement(){};

  /** 
   * Function for setting process input values
   * \param u input
   * \param Su input uncertainty
   * \param t input time
   */
  void set(InputValType u, InputUncertaintyType Su, double t = -1)
  {
    u_ = u;
    Su_ = Su;
    if(t >= 0)
      t_ = t;
  }

  /** 
   * Function for setting process input values
   * \param u input
   * \param t input time
   */
  void set(InputValType u, double t = -1)
  {
    u_ = u;
    if(t >= 0)
      t_ = t;
  }

  /** 
   * Function for getting process input values
   * \param u input [overwritten]
   * \param Su input uncertainty [overwritten]
   * \param t time of input [overwritten]
   */
  void get(InputValType &u, InputUncertaintyType &Su, double &t){
    u = u_;
    Su = Su_;
    t = t_;
  }

  /** 
   * Function for getting process input values
   * \param u input [overwritten]
   * \param t time of input [overwritten]
   */
  void get(InputValType &u, double &t){
    u = u_;
    t = t_;
  }


private:

  double t_; /**< time of the input */
  InputValType u_;  /**< Input */
  InputUncertaintyType Su_; /**< Input covariance */

};


/**
 * \class ProcessModel
 * \brief An abstract class for defining vehicle motion models
 * \tparam StateType Object type for state 
 * \tparam InputType Object type for process input information
 * \author Keith Leung
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



/********** Classes for an examples 2D motion model **********/

/**
 * \class Odometry2d
 * \brief A class for 2d odometry measurements for a 2d motion model
 * \author Keith Leung
 */
class Odometry2d : public Measurement< Eigen::Vector3d, Eigen::Matrix3d >
{
public:
  
  /** Default constructor */
  Odometry2d();
  
  /** 
   * Constructor, not necessary, but defined for convenience
   * \param dx_k_km x-displacement from frame k-1
   * \param dy_k_km y-displacement from frame k-1
   * \param dtheta_k_km rotational displacement from frame k-1
   * \param vardx_k_km variance in dx_k_km
   * \param vardy_k_km variance in dy_k_km
   * \param vardtheta_k_km variance in dtheta_k_km
   * \param t time of odometry reading
   */
  Odometry2d(double dx_k_km, double dy_k_km, double dtheta_k_km,
	     double vardx_k_km, double vardy_k_km, double vartheta_k_km,
	     double t);

  /** Destructor */
  ~Odometry2d();
};

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
