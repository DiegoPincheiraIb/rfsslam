// Measurement model
// Felipe Inostroza, Keith Leung 2013

#ifndef MEASUREMENTMODEL_HPP
#define MEASUREMENTMODEL_HPP

#include "Measurement.hpp"
#include "Landmark.hpp"
#include "Pose.hpp"
#include "RandomVecMathTools.hpp"

/** 
 * \class MeasurementModel
 * \brief An abstract class for defining the measurement model
 * \author Felipe Inostroza, Keith Leung
 */
template<class PoseType, class LandmarkType, class MeasurementType>
class MeasurementModel
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef PoseType TPose;
  typedef LandmarkType TLandmark;
  typedef MeasurementType TMeasurement;
  
  /** Default constructor */
  MeasurementModel() : nd_(0, 1), gen_(rng_, nd_), R_( MeasurementType::Mat::Zero()) {}

  /** Default destructor */
  ~MeasurementModel(){}

  /** Set the zero-mean-white-Gaussian noise covariance matrix for this model
   *  \param[in] R covariance matrix
   */
  void setNoise( typename MeasurementType::Mat &R ){
    R_ = R;
    if( R_ != MeasurementType::Mat::Zero() ){
      Eigen::LLT<typename MeasurementType::Mat> cholesky(R_);
      L_ = cholesky.matrixL();
    }
  }

  /** Get the zero-mean-white-Gaussian noise covariance matrix for this model
   *  \param[out] R covariance matrix
   */
  void getNoise( typename MeasurementType::Mat &R ){
    R = R_;
  }
  

  /** 
   * Abstract function for predicting measurements using pose and landmark estimates
   * \note This must be implemented in a derived class
   * \param[in] pose robot pose from which the measurement is made
   * \param[in] landmark The measured landmark
   * \param[out] measurement The measurement
   * \param[out] jacobian If not NULL, the pointed-to matrix will be overwritten 
   * by the Jacobian of the measurement model at the point where the prediction is made
   * \return true if a measurement can be produced
   */
  virtual bool measure( PoseType &pose, LandmarkType &landmark, 
			MeasurementType &meaurement, 
			Eigen::Matrix<double , 
				      MeasurementType::Vec::RowsAtCompileTime ,
				      LandmarkType::Vec::RowsAtCompileTime > 
			*jacobian = NULL ) = 0;
  
  /**
   * Sample a measurement with noise
   * \param[in] pose robot pose 
   * \param[in] landmark The measured landmark
   * \param[out] measurement Sampled measurement (which does not contain uncertainty information)
   * \param[in] useAdditiveWhiteGaussianNoise include the zero-mean white 
   * Gaussian noise set for this model
   * \param[in] usePoseWhiteGaussianNoise include the noise set for the pose 
   * and interpret as zero-mean-white-Gaussian noise
   * \param[in] useLandmarkWhiteGaussianNoise include the noise set for 
   * landmark and interpret as zero-mean-white-Gaussian noise
   * \return true if sucessfully sampled a measurement
   */
  bool sample( PoseType &pose, LandmarkType &landmark, 
	       MeasurementType &measurement,
      	       bool useAdditiveWhiteGaussianNoise = true,		       
	       bool usePoseWhiteGaussianNoise = false,
	       bool useLandmarkWhiteGaussianNoise = false){
	
    typename MeasurementType::Vec z;
    typename MeasurementType::Vec noise;

    PoseType* pose_sample = &pose;
    LandmarkType* landmark_sample = &landmark;
    bool deallocatePose = false;
    bool deallocateLandmark = false;
    
    if(usePoseWhiteGaussianNoise){

      pose_sample = new PoseType;
      deallocatePose = true;      
      RandomVecMathTools<PoseType>::sample(pose, *pose_sample);
    }

    if(useLandmarkWhiteGaussianNoise){

      landmark_sample = new LandmarkType;
      deallocateLandmark = true;
      RandomVecMathTools<LandmarkType>::sample(landmark, *landmark_sample);
    }

    bool success = this->measure( *pose_sample, *landmark_sample, measurement);

    if(success){
      if(useAdditiveWhiteGaussianNoise){

	typename MeasurementType::Vec z_k;
	double t;
	measurement.get(z_k, t );
	RandomVecMathTools<MeasurementType>::sample(z_k, R_, L_, t, measurement);
      }
    }

    
    if(deallocatePose)
      delete pose_sample;
    if(deallocateLandmark)
      delete landmark_sample;

    return success;

  }

  /** 
   * Abstract function for the inverse measurement model
   * \note This must be implemented in a derived class
   * \param[in] pose robot pose 
   * \param[in] measurement measurement
   * \param[out] landmark position
   */
  virtual void inverseMeasure( PoseType &pose,
			       MeasurementType &measurement, 
			       LandmarkType &landmark ) = 0;

  /**
   * Abstract function of determining a landmark's probability of detection
   * \note If this is not reimplemented in a derived class, it will always
   * return a probability of detection of 1
   * \brief In addition to specifying a single probability of detection as the
   * return value, we can also specify whether a landmark is close to the sensing limit,
   * where there can be mismatch between the probability of detection returned by this
   * function, and the actual probability of detection in real life. This flag is read 
   * by the RB-PHD-Filter and is important in its implementation to prevent landmarks
   * near the sensing limit from dissapearing due to the probability of detection mismatch.
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isClostToSensingLimit true if landmark is close to the sensing limit,
   * where there can be a mismatch between the returned value and the actual probability
   * of detection in reality.
   * \return probability of detection
   */
  virtual double probabilityOfDetection( PoseType &pose,
					 LandmarkType &landmark,
					 bool &isCloseToSensingLimit){
    isCloseToSensingLimit = false;
    return 1;
  }


  /**
   * Abstract function for Determining the clutter intensity in the measurement space
   * \note This should be reimplemented in a derived class
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  virtual double clutterIntensity( MeasurementType &z,
				   int nZ ){
    return 0;
  }

  /**
   * Abstract function for Determining the clutter intensity integral
   * \note This should be reimplemented in a derived class
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity integral
   */
  virtual double clutterIntensityIntegral( int nZ ){
    return 0;
  }  

protected:

  typename MeasurementType::Mat R_; /**< Additive zero-mean Gaussian noise covariance */
  typename MeasurementType::Mat L_; /** Lower triangular part of Cholesky decomposition on R_ */


  
  /** Generate a random number from a normal distribution */
  double randn(){
    return gen_();
  }
  boost::mt19937 rng_;
  boost::normal_distribution<double> nd_;
  boost::variate_generator< boost::mt19937, boost::normal_distribution<double> > gen_;

};


////////// 2d Range-Bearing Measurement Model //////////

/** 
 * \class RangeBearingModel
 * \brief Range and bearing measurement model for 2d point landmarks.
 * The noise is considered to be Gaussian in the range and bearing space
 * \author Felipe Inostroza 
 */
                                                               
class RangeBearingModel : public MeasurementModel <Pose2d, Landmark2d, Measurement2d>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for the 2d RangeBearingModel */
  struct Config{
    double probabilityOfDetection_; /** probability of detection */
    double uniformClutterIntensity_; /** clutter per area */
    double rangeLimMax_; /**< sensing range limit, beyond which Pd = 0 */
    double rangeLimMin_; /**< sensing range limit, beyond which Pd = 0 */
    double rangeLimBuffer_; /**< A buffer for Pd ambiguity */
  }config;

 /** Default constructor */
  RangeBearingModel();

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model
  * \param covZ measurement covariance
  */
  RangeBearingModel(Eigen::Matrix2d &covZ);

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model, 
  * range and bearing are assumed to be uncorrelated
  * \param Sr Range variance
  * \param Sb Bearing variance
  */
  RangeBearingModel(double Sr, double Sb);

 /** Default destructor */
  ~RangeBearingModel();

  /** 
   * Get a measurement
   * \param[in] pose robot pose from which the measurement is made
   * \param[in] landmark The measured landmark
   * \param[out] measurement The measurement
   * \param[out] jacobian If not NULL, the pointed-to matrix will be overwritten 
   * by the Jacobian of the measurement model at the point where the prediction is made
   * \return true if measurement is generated
   */
  bool measure( Pose2d &pose, Landmark2d &landmark, 
		Measurement2d &measurement, Eigen::Matrix2d *jacobian = NULL);

  /** 
   * Inverse measurement
   * \param[in] pose robot pose 
   * \param[in] measurement measurement
   * \param[out] landmark position
   */
  void inverseMeasure(Pose2d &pose, Measurement2d &measurement, Landmark2d &landmark);

  /** 
   * Determine the probability of detection
   * \note This is where we can indirectly specify sensing limits and other sensor characteristics
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isCloseToSensingLimit true if within range limit buffer zone
   * \return probability of detection
   */
  double probabilityOfDetection( Pose2d &pose,
				 Landmark2d &landmark,
				 bool &isCloseToSensingLimit);

  /**
   * Determine the clutter intensity in measurement space
   * \note uniform clutter intensity is assumed
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  double clutterIntensity( Measurement2d &z,
			   int nZ );

  /**
   * Determine the clutter intensity integral in measurement space
   * \brief This is calculated based on the probablity of false alarm,
   * defined as p( NULL | measurement exists)
   * \param[in] nZ the cardinality of Z
   * \return clutter intensity
   */
  double clutterIntensityIntegral( int nZ = 0);

};



////////// 1d Measurement Model //////////

/** 
 * \class  MeasurementModel1d
 * \brief Range measurement model for 1d point landmarks with Gaussian noise.
 * \author Keith Leung
 */
                                                               
class MeasurementModel1d : public MeasurementModel <Pose1d, Landmark1d, Measurement1d>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for the 2d RangeBearingModel */
  struct Config{
    double probabilityOfDetection_; /** probability of detection */
    double uniformClutterIntensity_; /** clutter per length */
    double rangeLimMax_; /**< sensing range limit, beyond which Pd = 0 */
    double rangeLimMin_; /**< sensing range limit, beyond which Pd = 0 */
    double rangeLimBuffer_; /**< A buffer for Pd ambiguity */
  }config;

 /** Default constructor */
  MeasurementModel1d();

 /**
  * Constructor that sets the uncertainty of the measurement model
  * \param Sr measurement variance
  */
  MeasurementModel1d(Eigen::Matrix<double, 1, 1> &Sr);

 /**
  * Constructor that sets the uncertainty of the measurement model
  * \param Sr Range variance
  */
  MeasurementModel1d(double Sr);

 /** Default destructor */
  ~MeasurementModel1d();

  /** 
   * Get a measurement
   * \param[in] pose robot pose from which the measurement is made
   * \param[in] landmark The measured landmark
   * \param[out] measurement The measurement
   * \param[out] jacobian If not NULL, the pointed-to matrix will be overwritten 
   * by the Jacobian of the measurement model at the point where the prediction is made
   * \return true if measurement is generated
   */
  bool measure( Pose1d &pose, Landmark1d &landmark, 
		Measurement1d &measurement, Eigen::Matrix<double, 1, 1> *jacobian = NULL);

  /** 
   * Inverse measurement
   * \param[in] pose robot position 
   * \param[in] measurement measurement
   * \param[out] landmark position
   */
  void inverseMeasure(Pose1d &pose, Measurement1d &measurement, Landmark1d &landmark);

  /** 
   * Determine the probability of detection
   * \note This is where we can indirectly specify sensing limits and other sensor characteristics
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] isCloseToSensingLimit true if in range limit buffer zone
   * \return probability of detection
   */
  double probabilityOfDetection( Pose1d &pose,
				 Landmark1d &landmark,
				 bool &isCloseToSensingLimit);

  /**
   * Determine the clutter intensity in measurement space
   * \note uniform clutter intensity is assumed
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  double clutterIntensity( Measurement1d &z,
			   int nZ );

  /**
   * Determine the clutter intensity integral in measurement space
   * \brief This is calculated based on the probablity of false alarm,
   * defined as p( NULL | measurement exists)
   * \param[in] nZ the cardinality of Z
   * \return clutter intensity
   */
  double clutterIntensityIntegral( int nZ = 0);

};



#endif
