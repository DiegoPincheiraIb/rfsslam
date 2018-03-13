#ifndef MEASUREMENTMODELDESCRIPTOR_HPP
#define MEASUREMENTMODELDESCRIPTOR_HPP

#include "MeasurementModel.hpp"
#include "MeasurementWithDescriptor.hpp"
#include "LandmarkWithDescriptor.hpp"
#include "Pose.hpp"

namespace rfs{

/**
 * \class MeasurementModelWithDescriptor
 * \brief A class to add a descriptor to an existing measurement model
 *
 * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f]
 * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
 * \tparam PoseType sensor pose type
 * \author Felipe Inostroza
 */
template<class MeasurementModelType, class DescriptorType>
class MeasurementModelWithDescriptor: public MeasurementModelType
{
public:

  //inherit constructors
  using MeasurementModelType::MeasurementModelType;


  typedef typename MeasurementModelType::TPose TPose;
  typedef typename MeasurementModelType::TLandmark TLandmarkNoDesc;
  typedef typename MeasurementModelType::TMeasurement TMeasurementNoDesc;

  typedef LandmarkWithDescriptor<TLandmarkNoDesc , DescriptorType> TLandmark;
  typedef MeasurementWithDescriptor< TMeasurementNoDesc , DescriptorType> TMeasurement;


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW




  /**
   * Abstract function for predicting a measurement from a robot pose and a landmark position
   * \f[ \mathbf{z} = \mathbf{h}(\mathbf{x}, \mathbf{m} ) + \mathbf{e}, \mathbf{e} \sim (\mathbf{0}, \mathbf{R}) \f]
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position, \f$\mathbf{e}\f$ is the zero-mean Gaussian noise.
   * \note This must be implemented in a derived class
   * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
   * \param[out] measurement \f$\mathbf{x}\f$, the measurement
   * \param[out] jacobian_wrt_lmk if not NULL, the pointed-to matrix is overwritten
   * by the Jacobian of the measurement model w.r.t. the landmark state evaluated at
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \param[out] jacobian_wrt_pose if not NULL, the pointed-to matrix is overwritten
   * by the Jacobian of the measurement model w.r.t. the robot state evaluated at
   * \f$\mathbf{x}\f$ and \f$\mathbf{m}\f$
   * \return true if a valid measurement is produced
   */
  bool measure( const TPose &pose,
                        const TLandmark &landmark,
                        TMeasurement &measurement,
                        ::Eigen::Matrix<double,
                         MeasurementModelType::TMeasurement::Vec::RowsAtCompileTime ,
                         MeasurementModelType::TLandmark::Vec::RowsAtCompileTime> *jacobian_wrt_lmk = NULL,
                        ::Eigen::Matrix<double,
                         MeasurementModelType::TMeasurement::Vec::RowsAtCompileTime,
                         MeasurementModelType::TPose::Vec::RowsAtCompileTime> *jacobian_wrt_pose = NULL
                        ) {
    bool success = MeasurementModelType::measure(pose ,landmark , measurement , jacobian_wrt_lmk , jacobian_wrt_pose);
    measurement.desc = landmark.desc;
    return success;
  };

  /**
   * Sample a measurement with noise parameters from the model, robot pose, and landmark position.
   * \param[in] pose \f$\mathbf{x}\f$, robot pose from which the measurement is made
   * \param[in] landmark \f$\mathbf{m}\f$, the measured landmark
   * \param[out] measurement sampled measurement (which does not contain uncertainty information, i.e., the covariance information is not set for this RandomVec)
   * \param[in] useAdditiveWhiteGaussianNoise if true, sample from the model's noise covariance, \f$\mathbf{R}\f$
   * \param[in] usePoseWhiteGaussianNoise if true, sample the uncertain robot pose using its covariance,
   * \f$\boldsymbol{\Sigma}_{\mathbf{x}}\f$, and interpret as zero-mean-white-Gaussian noise
   * \param[in] useLandmarkWhiteGaussianNoise if true, sample the uncertain landmark position using its covariance,
   * \f$\boldsymbol{\Sigma}_{\mathbf{m}}\f$, and interpret as zero-mean-white-Gaussian noise
   * \return true if sucessfully sampled a measurement
   */
  bool sample( TPose &pose, TLandmark &landmark,
               TMeasurement &measurement,
               bool useAdditiveWhiteGaussianNoise = true,
               bool usePoseWhiteGaussianNoise = false,
               bool useLandmarkWhiteGaussianNoise = false){

    TPose pose_sample; /**< Sampled pose*/
    TLandmark landmark_sample; /**< Sampled landmark*/

    if(usePoseWhiteGaussianNoise){
      pose.sample(pose_sample);
    }else{
      pose_sample = pose;
    }

    if(useLandmarkWhiteGaussianNoise){
      landmark.sample(landmark_sample);
    }else{
      landmark_sample = landmark;
    }

    bool success = this->measure( pose_sample, landmark_sample, measurement);

    if(success && useAdditiveWhiteGaussianNoise){
      measurement.setCov(MeasurementModelType::R_);
      measurement.sample();
    }

    return success;

  }

  /**
   * Abstract function for the inverse measurement model
   * \f[ \mathbf{m} = \mathbf{h}^{-1}(\mathbf{x}, \mathbf{z} )\f]
   * where \f$\mathbf{z}\f$ is a measurement, \f$\mathbf{x}\f$ is the robot pose, \f$\mathbf{m}\f$ is a landmark position
   * \note This must be implemented in a derived class, and both the mean and the covariance of \f$\mathbf{m}\f$ should
   * be calculated. This is used in the RBPHDFilter for generating birth Gaussians
   * \param[in] pose \f$\mathbf{x}\f$, robot pose (the uncertainty is not used here because the pose_sample_
   * RBPHDFilter represents robot pose estimates with particles)
   * \param[in] measurement  \f$\mathbf{z}\f$ measurement, for which the uncertainty is \f$\mathbf{R}\f$
   * \param[out] landmark  \f$\mathbf{m}\f$, predicted landmark position with uncertainty
   */
  void inverseMeasure( const TPose &pose,
                               const TMeasurement &measurement,
                               TLandmark &landmark ){
    MeasurementModelType::inverseMeasure(pose, measurement  , landmark );
    landmark.desc = measurement.desc ;
  };

  double clutterIntensity( const TMeasurement &z,
                                   int nZ ){
    return MeasurementModelType::clutterIntensity(z , nZ) * z.desc.falseAlarmLikelihood();
  }









};

}
#endif
