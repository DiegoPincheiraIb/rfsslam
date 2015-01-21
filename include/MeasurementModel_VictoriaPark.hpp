#ifndef MEASUREMENTMODEL_VICTORIAPARK_HPP
#define MEASUREMENTMODEL_VICTORIAPARK_HPP







#include <vector>
#include <rfsslam/MeasurementModel.hpp>
#include <rfsslam/MeasurementModel_RngBrg.hpp>
using namespace rfs;






////////// 2d Range-Bearing Measurement Model Extended for Circle Detection using a laser//////////

/**
 * \class MeasurementModel_VictoriaPark
 * A range and bearing measurement model circle detection using a Lidar
 *
 * \brief A 2d range-bearing measurement model for 2d Circle Features
 * \author Felipe Inostroza
 */



class MeasurementModel_VictoriaPark : public MeasurementModel <Pose2d, Landmark3d, Measurement3d>{

public:
  /** Default constructor */
  MeasurementModel_VictoriaPark();

  /**
    * Constructor that sets the uncertainty (covariance) of the measurement model, \f$\mathbf{R}\f$
    * \param covZ measurement covariance, \f$\mathbf{R}\f$
    * \param Slb laser bearing variance used in calculating the diameter variance
    */
  MeasurementModel_VictoriaPark(Eigen::Matrix3d &covZ, double Slb);


  /**
    * Constructor that sets the uncertainty (covariance) of the measurement model, \f$\mathbf{R}\f$
    * \param covP measurement covariance for the range bearing model, \f$\mathbf{R_p}\f$
    * \param covR measurement covariance for the diameter, \f$\mathbf{R}_R\f$
    * \param Slb laser bearing variance used in calculating the diameter variance
    */
  MeasurementModel_VictoriaPark(Eigen::Matrix2d &covP, double covR, double Slb);

  /**
   * Constructor that sets the uncertainty (covariance) of the measurement model,
   * \f[\mathbf{R} = \begin{bmatrix} \sigma_r^2 & 0 & 0 \\ 0 & \sigma_b^2 & 0 \\ 0 & 0 & \sigma_R^2 \end{bmatrix}\f]
   * radius, range and bearing are assumed to be uncorrelated
   * \param Sr Range variance \f$\sigma_r^2\f$
   * \param Sb Bearing variance \f$\sigma_b^2\f$
   * \param Sd flat diameter variance \f$\sigma_R^2\f$
   * \param Slb laser bearing variance used in calculating the diameter variance
   */
   MeasurementModel_VictoriaPark(double Sr, double Sb , double Sd, double Slb);

   /** Default destructor */
    ~MeasurementModel_VictoriaPark();



   /**
    * Set the zero-mean-white-Gaussian additive noise covariance matrix, \f$\mathbf{R}\f$
    * \param[in] R covariance matrix
    * \param Slb laser bearing variance used in calculating the diameter variance
    */
   void setNoise( Measurement3d::Mat &R , double Slb);

   /**
    * Obtain a measurement from a given robot pose and landmark position
    * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f]
    * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
    * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
    * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
    * \param[out] measurement \f$\mathbf{x}\f$, the measurement
    * \param[out] jacobian_wrt_lmk if not NULL, the pointed-to matrix is overwritten
    * by the Jacobian of the measurement model, \f$\mathbf{H}\f$, evaluated at \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
    * \return true if a valid measurement is produced
    */
   bool measure( const Pose2d &pose, const Landmark3d &landmark,
                 Measurement3d &measurement, Eigen::Matrix3d *jacobian_wrt_lmk = NULL, Eigen::Matrix3d *jacobian_wrt_pose = NULL);

   /**
    * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f]
    * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
    * \param[in] pose \f$\mathbf{x}\f$, robot pose (the uncertainty is not used here because the
    * RBPHDFilter represents robot pose estimates with particles)
    * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
    * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
    */
   void inverseMeasure(const Pose2d &pose, const Measurement3d &measurement, Landmark3d &landmark);


   /**
    * Abstract function of determining a landmark's probability of detection, and if the landmark is close to the sensing limit.
    * Through this we can indirectly specify sensing limits and other sensor characteristics
    * The probability of detection is necessary as a parameter is the PHD Filter. Indicating whether a landmark is close to the
    * sensing limit matters in the implementation for providing a better map estimate, as it reduces landmark disappearance
    * near the sensing limit due to the probability of detection mismatch.
    * \param[in] pose robot pose
    * \param[in] landmark landmark position
    * \param[out] isCloseToSensingLimit true if lahttps://www.google.com/search?client=ubuntu&channel=fs&q=23%2F30*6%2B1&ie=utf-8&oe=utf-8#channel=fs&q=twitter+%23shamelessplug&safe=offndmark is close to the sensing limit
    * \return probability of detection
    */
   double probabilityOfDetection( const Pose2d &pose,
                                  const Landmark3d &landmark,
                                  bool &isCloseToSensingLimit);

   double probabilityOfDetection2( const Pose2d &pose,
                                  const Landmark3d &landmark,
                                  bool &isCloseToSensingLimit);

   /**
    * Determine the clutter intensity in measurement space.
    * Uniform clutter intensity is assumed
    * \param[in] z measurement point at which clutter intensity will be determined
    * \param[in] nZ the cardinality of Z, of which z is a member.
    * \return clutter intensity
    */
   double clutterIntensity( Measurement3d &z,
                            int nZ );

   /**
    * Determine the clutter intensity integral in measurement space.
    * This is calculated based on the probablity of false alarm,
    * defined as p( NULL | measurement exists)
    * \param[in] nZ the cardinality of Z
    * \return clutter intensity
    */
   double clutterIntensityIntegral( int nZ = 0);




   void setLaserscan(const std::vector<double> &laserscan){

     laserscan_=laserscan;

     double FoVArea=0;
     for(int i=1;i<laserscan_.size();i++){
     
       FoVArea+=laserscan_[i]*laserscan_[i-1];
     }
     FoVArea+=laserscan_[0]*laserscan_[laserscan_.size()-1];
     FoVArea*=sin(PI/360)/2;
     
     clutterIntensity_=config.expectedClutterNumber_/FoVArea;
     
   }




   /** \brief Configuration for this 2d Circle  Measurement Model */
   struct Config{
     std::vector<double> probabilityOfDetection_; /**< Array containing probabilities of detection for a circle given the amount of points that should show up on the scan */
     double expectedClutterNumber_; /**< Expected number of clutter measurements, the intensity is asummed constant over the field of view of the robot */
     double rangeLimMax_; /**< sensing range limit, beyond which \f$ P_D = 0 \f$*/
     double rangeLimMin_; /**< sensing range limit, below which \f$ P_D = 0 \f$*/
     double bearingLimitMin_; /**< sensing angle limit, below which \f$ P_D = 0 \f$*/
     double bearingLimitMax_; /**< sensing angle limit, beyond which \f$ P_D = 0 \f$*/

   }config;



private:

   double Slb_;
   MeasurementModel_RngBrg rangeBearingModel;
   std::vector<double> laserscan_;
   double clutterIntensity_;

};







#endif

