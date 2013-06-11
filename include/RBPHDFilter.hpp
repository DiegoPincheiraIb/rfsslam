// RB-PHD-Filter Class
// Keith Leung 2013

#ifndef RBPHDFILTER_HPP
#define RBPHDFILTER_HPP

#include <Eigen/Core>
#include "GaussianMixture.hpp"
#include "ParticleFilter.hpp"
#include <vector>

/**
 *  \class RBPHDFilter
 *  \brief Rao-Blackwellized Probability Hypothesis Density Filter class
 *  
 *  This class implements the Rao-Bloackwellized Probability Hypothesis Density
 *  filter. 
 *
 *  \author Keith Leung, Felipe Inostroza
 *  \version 0.1
 *
 *  \todo function for map update (but requires completion of Kalman Filter class)
 *  \todo function for importance weighting
 *  
 */

template< class ProcessModel, class MeasurementModel, class KalmanFilter>
class RBPHDFilter : public ParticleFilter<ProcessModel, MeasurementModel>
{
public:

  typedef typename ProcessModel::tState StateType;
  typedef typename ProcessModel::tInput InputType;
  typedef typename MeasurementModel::tLandmark LandmarkType;
  typedef typename MeasurementModel::tMeasurement MeasureType;

  /** Constructor */
  RBPHDFilter(int n, 
	      StateType &initState,
	      ProcessModel* processModelPtr,
	      MeasurementModel* measurementModelPtr)
    : ParticleFilter<ProcessModel, MeasurementModel>(n, initState, processModelPtr, measurementModelPtr){
    
    for(int i = 0; i < n; i++){
      maps_.push_back( new GaussianMixture<LandmarkType>() );
      nLandmarksBeforeUpdate_.push_back(0);
      nLandmarksAfterUpdate_.push_back(0);
    }

    default_birthGaussianWeight_ = 0.25;
  }

  /** Destructor */
  ~RBPHDFilter();

  /**
   * Update the map with the measurements in measurements_
   * Existing landmarks with probability of detection > 0 will have their Gaussian
   * mixture weight reduced to account for missed detection.
   * For every landmark-measurement pair with probability of detection > 0,
   * a new landmark will be created. 
   */
  void updateMap();

  /** 
   * Importance weighting. Overrides the abstract function in ParticleFilter
   */
  void importanceWeighting();

  /** 
   * Add birth Gaussians for each particle's map using unused_measurements_
   */ 
  void addBirthGaussians();

private:

  std::vector< GaussianMixture<LandmarkType>* > maps_; /**< Particle dependent maps */

  /** indices of unused measurement for each particle for creating birth Gaussians */
  std::vector< std::vector<unsigned int> > unused_measurements_; 

  /** Number of landmarks estimated for each particle before map update */
  std::vector< double > nLandmarksBeforeUpdate_;

  /** Number of landmarks estimated for each particle after map update */
  std::vector< double > nLandmarksAfterUpdate_;
  
  /** The Gaussian mixture weight assigned to new birth Gaussians */
  double default_birthGaussianWeight_; 

};

////////// Implementation //////////

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::~RBPHDFilter(){

  for(int i = 0; i < maps_.size(); i++){
    delete maps_[i];
  }
}

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::updateMap(){

  for(int i = 0; i < this->nParticles_; i++){

  }

  // \todo update unused_measurements
}


template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::importanceWeighting(){

  // \todo implement all three weighting strategies

  for(int i = 0; i < this->nParticles_; i++){

  }

}

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::addBirthGaussians(){

  for(int i = 0; i < this->nParticles_; i++){
    
    while( unused_measurements_[i].size() > 0){
     
      // get measurement
      int unused_idx = unused_measurements_[i].back();
      MeasureType unused_z = this->measurements_[unused_idx];
      unused_measurements_[i].pop_back();

      // use inverse measurement model to get landmark
      StateType robot_pose;
      LandmarkType landmark_pos;
      this->particleSet_[i]->getPose(robot_pose);
      this->pMeasurementModel_->inversePredict( robot_pose, landmark_pos, unused_z );

      // add birth landmark to Gaussian mixture (last param = true to allocate mem)
      maps_[i]->addGaussian( &landmark_pos, default_birthGaussianWeight_, true);
      
    }
    
  }

}

#endif
