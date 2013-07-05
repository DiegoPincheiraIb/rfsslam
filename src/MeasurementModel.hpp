// Measurement model
// Felipe Inostroza 2013

#ifndef MEASUREMENTMODEL_HPP
#define MEASUREMENTMODEL_HPP

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "Measurement.hpp"
#include "Landmark.hpp"
#include "Pose.hpp"
#include "Tools.hpp"

/** 
 * \class MeasurementModel
 * \brief An abstract class for defining the measurement model
 * \author Felipe Inostroza, Keith Leung
 */
template<class PoseType, class LandmarkType, class MeasurementType>
class MeasurementModel
{
public:

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
    if( R_ != Eigen::Matrix2d::Zero() ){
      Eigen::LLT<Eigen::Matrix2d> cholesky(R_);
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
   */
  virtual void measure( PoseType &pose, LandmarkType &landmark, 
			MeasurementType &meaurement, 
			Eigen::Matrix<double , 
				      MeasurementType::Vec::RowsAtCompileTime ,
				      LandmarkType::Vec::RowsAtCompileTime > 
			*jacobian = NULL ) = 0;
  
  /**
   * Sample a measurement with noise
   * \param[in] pose robot pose 
   * \param[in] landmark The measured landmark
   * \param[out] measurement Sampled measurement
   * \param[in] useAdditiveWhiteGaussianNoise include the zero-mean white 
   * Gaussian noise set for this model
   * \param[in] usePoseWhiteGaussianNoise include the noise set for the pose 
   * and interpret as zero-mean-white-Gaussian noise
   * \param[in] useLandmarkWhiteGaussianNoise include the noise set for 
   * landmark and interpret as zero-mean-white-Gaussian noise
   */
  void sample( PoseType &pose, LandmarkType &landmark, 
	       MeasurementType &measurement,
      	       bool useAdditiveWhiteGaussianNoise = true,		       
	       bool usePoseWhiteGaussianNoise = false,
	       bool useLandmarkWhiteGaussianNoise = false){
	
    typename PoseType::Vec z;
    typename PoseType::Vec noise;
    
    if(usePoseWhiteGaussianNoise){

      typename PoseType::Vec x;
      typename PoseType::Mat Sx, Sx_L;
      pose.get( x, Sx );
      Eigen::LLT<typename PoseType::Mat> cholesky( Sx );
      Sx_L = cholesky.matrixL();
    
      int n = Sx_L.cols();
      typename PoseType::Vec randomVecNormal, randomVecGaussian;
      for(int i = 0; i < n; i++){
	randomVecNormal(i) = randn();
      }
      randomVecGaussian = Sx_L * randomVecNormal;
      x = x + randomVecGaussian;
      pose.set( x, Sx );
    }

    if(useLandmarkWhiteGaussianNoise){

      typename LandmarkType::Vec m;
      typename LandmarkType::Mat Sm, Sm_L;
      landmark.get( m, Sm );
      Eigen::LLT<typename LandmarkType::Mat> cholesky( Sm );
      Sm_L = cholesky.matrixL();
    
      int n = Sm_L.cols();
      typename LandmarkType::Vec randomVecNormal, randomVecGaussian;
      for(int i = 0; i < n; i++){
	randomVecNormal(i) = randn();
      }
      randomVecGaussian = Sm_L * randomVecNormal;
      m = m + randomVecGaussian;
      landmark.set( m, Sm );

    }

    this->measure( pose, landmark, measurement);
    measurement.State<Eigen::Vector2d>::get(z);

    if(useAdditiveWhiteGaussianNoise){
      for(int i = 0; i < TPose::Vec::RowsAtCompileTime; i++){
	noise(i) = randn();
      }
      z += L_ * noise;
    }

    measurement.State<Eigen::Vector2d>::set(z);
    
  };

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
   * \param[in] z measurement point at which clutter intensity will be determined
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
    double probabilityOfDetection_;
    double probabilityOfFalseAlarm_; /**<  interpret as p( NULL | measurement exists) */
    double rangeLim_; /**< sensing range limit, beyond which Pd = 0 */
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
   */
  void measure( Pose2d &pose, Landmark2d &landmark, 
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
   * \param[out] Pd_upper
   * \param[out] Pd_lower
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
  double clutterIntensityIntegral( int nZ );

protected:
  
  double sensingArea_;

};



//////////  Linear Measurement Model //////////

/** 
 * \class LinearModel
 * \brief  This is a linear measurement model that measures directly the state 
 * of a landmark, this is made for testing purposes only. 
 * 
 * \author Felipe Inostroza 
 */
template <class LandmarkType , class MeasurementType> 
class LinearModel: public MeasurementModel <Pose2d, LandmarkType, MeasurementType>{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  /** \brief Configuration for the Model */
  struct Config{
    double probabilityOfDetection_;
    double probabilityOfFalseAlarm_; /**<  interpret as p( NULL | measurement exists) */
  }config;

 /** Default constructor 
  *  Sets the H matrix to the identity matrix (The matrix is cut to have the proper dimensions)
  */
  LinearModel(){

  config.probabilityOfDetection_ = 0.95;
  config.probabilityOfFalseAlarm_ = 0.01;
 
};

 /**
  * Constructor that sets the uncertainty (covariance) of the measurement model
  * \param covZ measurement covariance
  */
  LinearModel(typename MeasurementType::Mat &R , Eigen::Matrix<double, MeasurementType::Vec::RowsAtCompileTime, LandmarkType::Vec::RowsAtCompileTime> &H){
  
  config.probabilityOfDetection_ = 0.95;
  config.probabilityOfFalseAlarm_ = 0.01;
  H_=H;
  H_inv_=(H.transpose()*H).inverse()*H.transpose();
  setNoise(R);
};


 /** Default destructor */
  ~LinearModel(){};

  /** 
   * Get a measurement
   * \param[in] pose robot pose from which the measurement is made
   * \param[in] landmark The measured landmark
   * \param[out] measurement The measurement
   * \param[out] jacobian If not NULL, the pointed-to matrix will be overwritten 
   * by the Jacobian of the measurement model at the point where the prediction is made
   */
  void measure( Pose2d &pose, LandmarkType &landmark, 
		MeasurementType &measurement, Eigen::Matrix<double, 
		MeasurementType::Vec::RowsAtCompileTime, 
		LandmarkType::Vec::RowsAtCompileTime> *jacobian = NULL){
  
  typename LandmarkType::Mat S_land;
  typename LandmarkType::Vec land;
  typename MeasurementType::Vec z;
  typename MeasurementType::Mat Sz;

  landmark.get(land,S_land);
  z=H_*land;
  Sz=H_*S_land*H_.transpose();
  measurement.set(z, Sz);
  if(jacobian != NULL)
    *jacobian = H_;
};

  /** 
   * Inverse measurement
   * \param[in] pose robot pose 
   * \param[in] measurement measurement
   * \param[out] landmark position
   */
  void inverseMeasure(Pose2d &pose, MeasurementType &measurement, LandmarkType &landmark){
    typename MeasurementType::Vec z;
    typename MeasurementType::Mat Sz;  
    typename LandmarkType::Mat S_land;
    typename LandmarkType::Vec land;  

    measurement.get(z,Sz);
    land = H_inv_*z;
    S_land = H_inv_*Sz*H_inv_.transpose();
    landmark.set(land , S_land);
    
    
};

  /** 
   * Determine the probability of detection
   * \note This is where we can indirectly specify sensing limits and other sensor characteristics
   * \param[in] pose robot pose
   * \param[in] landmark landmark position
   * \param[out] Pd_upper
   * \param[out] Pd_lower
   * \return probability of detection
   */
  double probabilityOfDetection( Pose2d &pose,
				 LandmarkType &landmark,
				 bool &isCloseToSensingLimit){
    return config.probabilityOfDetection_;
  };

  /**
   * Determine the clutter intensity in measurement space
   * \note uniform clutter intensity is assumed
   * \param[in] z measurement point at which clutter intensity will be determined
   * \param[in] nZ the cardinality of Z, of which z is a member.
   * \return clutter intensity
   */
  double clutterIntensity( MeasurementType &z,
			   int nZ ){
    return config.probabilityOfFalseAlarm_;
  };

  /**
   * Determine the clutter intensity integral in measurement space
   * \brief This is calculated based on the probablity of false alarm,
   * defined as p( NULL | measurement exists)
   * \param[in] nZ the cardinality of Z
   * \return clutter intensity
   */
  double clutterIntensityIntegral( int nZ ){
    return config.probabilityOfFalseAlarm_*nZ;
  };

protected:
  
  Eigen::Matrix<double, MeasurementType::Vec::RowsAtCompileTime, LandmarkType::Vec::RowsAtCompileTime> H_;
  Eigen::Matrix<double, MeasurementType::Vec::RowsAtCompileTime, LandmarkType::Vec::RowsAtCompileTime> H_inv_;


};


#endif
