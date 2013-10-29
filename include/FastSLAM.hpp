/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FASTSLAM_HPP
#define FASTSLAM_HPP

#include <boost/timer/timer.hpp>
#include <Eigen/Core>
#include "GaussianMixture.hpp"
#include "LinearAssignment.hpp"
#include "KalmanFilter.hpp"
#include "ParticleFilter.hpp"
#include <math.h>
#include <vector>

#include <stdio.h>

/**
 *  \class FastSLAM
 *  \brief Factored Solution to SLAM
 *  
 *  This is an implementation of the FastSLAM v1.0 algorithm. 
 *
 *  @INPROCEEDINGS{Montemerlo02a,
 *  AUTHOR         = {Montemerlo, M. and Thrun, S. and Koller, D. and 
 *                   Wegbreit, B.},
 *  TITLE          = {{FastSLAM}: {A} Factored Solution to the Simultaneous 
 *                   Localization and Mapping Problem},
 *  YEAR           = {2002},
 *  BOOKTITLE      = {Proceedings of the AAAI National Conference on 
 *                   Artificial Intelligence},
 *  PUBLISHER      = {AAAI},
 *  ADDRESS        = {Edmonton, Canada}
 *  }
 *
 *  \tparam RobotProcessModel A robot process model derived from ProcessModel
 *  \tparam LmkProcessModel A landmark process model derived from ProcessModel
 *  \tparam MeasurementModel A sensor model derived from MeasurementModel
 *  \tparam KalmanFilter A Kalman filter that uses LmkProcessModel and MeasurementModel
 *  \author Keith Leung, Felipe Inostroza
 */

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
class FastSLAM : public ParticleFilter<RobotProcessModel, MeasurementModel>
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef typename RobotProcessModel::TState TPose;
  typedef typename RobotProcessModel::TInput TInput;
  typedef typename MeasurementModel::TLandmark TLandmark;
  typedef typename MeasurementModel::TMeasurement TMeasurement;
  typedef typename GaussianMixture<TLandmark>::Gaussian TGaussian;

  /** 
   * \brief Configurations for this RBPHDFilter 
   */
  struct Config{

    /** Gaussian pruning weight threshold, below which Gaussians are eliminated from a Gaussian mixture */
    double gaussianPruningThreshold_;

    /** Minimum timeteps betwen resampling of particles*/
    double minInterSampleTimesteps_;

    /** If true, timing information is written to the console every update*/
    bool reportTimingInfo_;

    /** The probability of existence that is initially assigned to a new landmark */
    double landmarkExistencePrior_;

    /** The log odds threshold for eliminating a landmark from the map*/
    double mapExistencePruneThreshold_;


  } config;

  /** 
   * Constructor 
   * \param n number of particles
   * \param initState initial state of particles
   */
  FastSLAM(int n);

  /** Destructor */
  ~FastSLAM();

  /** 
   * Get the landmark process model
   * \return pointer to the landmark process model
   */
  LmkProcessModel* getLmkProcessModel();

  /**
   * Predict the robot trajectory using the lastest odometry data
   * \param[in] input 
   * \param[in] currentTimestep current timestep;
   */
  void predict( TInput u, int currentTimestep );

  /**
   * Update the map, calculate importance weighting, and perform resampling if necessary
   * \param[in] Z set of measurements to use for the update, placed in a std vector, which
   * gets cleared after the function call. 
   * \param[in] currentTimestep current timestep;
   */
  void update( std::vector<TMeasurement> &Z, int currentTimestep );

  /**
   * Get the size of the Gaussian mixture for a particle
   * \param[in] i particle index
   * \return size if index is valid, else -1
   */
  int getGMSize(int i);

  /**
   * Get the position, covariance, and weight of a Gaussian in particle i's Gaussin mixture
   * \param[in] i particle index
   * \param[in] m Gaussian index
   * \param[out] u mean
   * \param[out] S covariance
   * \return false if the indices specified are invalid 
   */ 
  bool getLandmark(const int i, const int m, 
		   typename TLandmark::Vec &u,
		   typename TLandmark::Mat &S);
  
  /**
   * Get the pointer to the Kalman Filter used for updating the map
   * \return pointer to the Kalman Filter
   */
  KalmanFilter* getKalmanFilter();


private:

  KalmanFilter *kfPtr_; /**< pointer to the Kalman filter */
  LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */

  std::vector< GaussianMixture<TLandmark>* > maps_; /**< Particle dependent maps */

  /** indices of unused measurement for each particle for creating birth Gaussians */
  std::vector< std::vector<unsigned int> > unused_measurements_; 

  int k_currentTimestep_; /**< current time */
  int k_lastResample_; /**< last resample time */
  
  /** 
   * Add birth Gaussians for each particle's map using unused_measurements_
   */ 
  void addBirthGaussians();

  /**
   * Update the map with the measurements in measurements_
   * Existing landmarks with probability of detection > 0 will have their Gaussian
   * mixture weight reduced to account for missed detection.
   * For every landmark-measurement pair with probability of detection > 0,
   * a new landmark will be created. 
   */
  void updateMap();

  /**
   * Resample the particles, along with their individual maps,  according to their 
   * importance weightings.
   */
  void resample();

  /** 
   * Importance weighting. Overrides the abstract function in ParticleFilter
   */
  void importanceWeighting();

  /**
   * Random Finite Set measurement likelihood evaluation
   * \brief The current measurements in measurements_ are used to determine the
   * RFS measurement likelihood given a set of landmarks 
   * \param[in] particleIdx particle for which the likelihood is calcuated
   * \param[in] indices of evaluation points in maps_[particleIdx]
   * \param[in] probability of detection of evaluation point 
   * \return measurement likelihood
   */
  double rfsMeasurementLikelihood( const int particleIdx, 
				   std::vector<unsigned int> &evalPtIdx,
				   std::vector<double> &evalPtPd );

  /**
   * Calculate the sum of all permutations of measurement likelihood from a likelihood
   * table generated from within rfsMeasurementLikelihood
   * \param[in] likelihoodTab likelihood table generated within rfsMeasurementLikelihood
   * \para,[in] A vector of measurement indices (columns) to consider in the likelihoodTab 
   * \return sum of all permutations from the given likelihood table
   */
  double rfsMeasurementLikelihoodPermutations( std::vector< double* > &likelihoodTab, 
					       std::vector< int > &Z_NoClutter);

  /** Checks the Gaussian mixture maps for all particles for errors
   *  \return true if there are no errors 
   */
  bool checkMapIntegrity();

};

////////// Implementation //////////

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::FastSLAM(int n)
  : ParticleFilter<RobotProcessModel, MeasurementModel>(n)
{

  lmkModelPtr_ = new LmkProcessModel;
  kfPtr_ = new KalmanFilter(lmkModelPtr_, this->getMeasurementModel());
  
  for(int i = 0; i < n; i++){
    // printf("Creating map structure for particle %d\n", i);
    maps_.push_back( new GaussianMixture<TLandmark>() );
  }
  
  config.gaussianPruningThreshold_ = 0.2;
  config.minInterSampleTimesteps_ = 5;
  config.reportTimingInfo_ = false;
  config.landmarkExistencePrior_ = 0.5;
  config.mapExistencePruneThreshold_ = -3.0;
  
  k_lastResample_ = -10;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::~FastSLAM(){
  for(int i = 0; i < maps_.size(); i++){
    delete maps_[i];
  }
  delete kfPtr_;
  delete lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
LmkProcessModel* FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getLmkProcessModel(){
  return lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::predict( TInput u,
											      int currentTimestep){

  boost::timer::auto_cpu_timer *timer = NULL;
  if(config.reportTimingInfo_)
    timer = new boost::timer::auto_cpu_timer(6, "Predict time: %ws\n");

  // propagate particles
  k_currentTimestep_ = currentTimestep;
  this->propagate(u);

  // propagate landmarks
  for( int i = 0; i < this->nParticles_; i++ ){
    for( int m = 0; m < maps_[i]->getGaussianCount(); m++){
      TLandmark *plm;
      maps_[i]->getGaussian(m, plm);
      lmkModelPtr_->staticStep(*plm, *plm);
    }
  }

  if(timer != NULL)
    delete timer;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::update( std::vector<TMeasurement> &Z,
												int currentTimestep){

  boost::timer::auto_cpu_timer *timer_mapUpdate = NULL;
  boost::timer::auto_cpu_timer *timer_particleWeighting = NULL;
  boost::timer::auto_cpu_timer *timer_mapMerge = NULL;
  boost::timer::auto_cpu_timer *timer_mapPrune = NULL;
  boost::timer::auto_cpu_timer *timer_particleResample = NULL;

  k_currentTimestep_ = currentTimestep;

  this->setMeasurements( Z ); // Z gets cleared after this call, measurements now stored in this->measurements_

  ////////// Map Update and Particle Weighting//////////
  if(config.reportTimingInfo_){
    timer_mapUpdate = new boost::timer::auto_cpu_timer(6, "Map update time: %ws\n");
  }
  updateMap();
  if(timer_mapUpdate != NULL)
    delete timer_mapUpdate;

  //////////// Particle resampling //////////
  if(config.reportTimingInfo_){
    timer_particleResample = new boost::timer::auto_cpu_timer(6, "Particle resample time: %ws\n");
  }
  resample();
  if(timer_particleResample != NULL)
    delete timer_particleResample;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::updateMap(){

  const unsigned int startIdx = 0;
  const unsigned int stopIdx = this->nParticles_;

  const unsigned int nZ = this->measurements_.size();

  HungarianMethod hm(); // for data association

  for(unsigned int i = startIdx; i < stopIdx; i++){    

    //----------  1. Data Association --------------------

    // Loglikelihood table for Data Association
    const unsigned int nM = maps_[i]->getGaussianCount();
    const unsigned int nMZ = nM;
    if(nZ > nM)
      nMZ = nZ;
    double** likelihoodTable = new double* [ nMZ ];
    for( int m = 0; m < nMZ; m++ ){
      likelihoodTable[m] = new double [ nMZ ];
      for(int z = 0; z < nMZ; z++){
	likelihoodTable[m][z] = -10;
      }
    }

    TPose *pose = new TPose;
    this->particleSet_[i]->getPose(*pose);

    for(unsigned int m = 0; m < nM; m++){

      TLandmark* lm = maps_[i]->getGaussian(m); // landmark position estimate
      TMeasurement measurement_exp; // expected measurement      
      bool isValidExpectedMeasurement = this->pMeasurementModel_->measure( pose , lm , measurement_exp);

      int z = 0;
      for(z = 0; z < nZ; z++){

	if( isValidExpectedMeasurement ){
	  likelihoodTable[m][z] = log(measurement_exp.evalGaussianLikelihood(this->measurements_[z])); 
	}
      }

    }

    int da[nMZ]; // data association result
    double logLikelihoodSum = 0; // used for particle weighting
    hm.run(likelihoodTable, nMZ, da, &logLikelihoodSum);
    

    //----------  2. Kalman Filter map update ----------

    bool zUsed[nZ];
    for(unsigned int z = 0; z < nZ; z++){
      zUsed[z] = false;
    }
    double DA_THRESHOLD = -1;
    
    double nExpectedClutter = this->pMeasurementModel_->clutterIntensityIntegral(nZ);
    double probFalseAlarm = nExpectedClutter / nZ;
    double p_exist_given_Z = 0;

    for(unsigned int m = 0; m < nM; m++){
      
      TLandmark* lm = maps_[i]->getGaussian(m); // landmark position estimate before update

      // Get the probability of detection, for map management caluclations
      bool temp;
      double probDetect = this->pMeasurementModel_->probabilityOfDetection(*pose, *lm, temp); 

      // Update landmark estimate m with the associated measurement
      int z = da[m];
      if( z < nZ && likelihoodTable[m][z] > DA_THRESHOLD){

	zUsed[z] = true; // This flag is for new landmark creation

	// EKF Update
	kfPtr_->correct(*pose, this->measurements_[z], *lm, *lm);

	// Particle weighting
	// logLikelihoodSum += likelihoodTable[m][z];

	// calculate change to existence probability
	p_exist_given_Z = ((1 - probDetect) * probFalseAlarm * config.landmarkExistencePrior_ + probDetect * config.landmarkExistencePrior_)/
	  (probFalseAlarm + (1 - probFalseAlarm) * probDetect * config.landmarkExistencePrior_); // 

      }else{ // landmark estimate m not updated

	// calculate change to existence probability
	p_exist_given_Z = ((1 - probDetect) * config.landmarkExistencePrior_) /
	  ((1 - config.landmarkExistencePrior_) + (1 - probDetect) * config.landmarkExistencePrior_);

      }

      double w = maps_[i]->getWeight(m); // this is the log-odds of existence given previous measurements
      w += log( (p_exist_given_Z) / (1 - p_exist_given_Z) );
      maps_[i]->setWeight(w);
      
    }

    //---------- 3. Map Management (Add and remove landmarks)  ------------

    maps_[i]->prune(config.mapExistencePruneThreshold_); 

    for(unsigned int z = 0; z < nZ; z++){
      if(!zUsed[z]){ // Create new landmarks with inverse measurement model with unused measurements

	TLandmark landmark_pos;
	this->pMeasurementModel_->inverseMeasure( *pose, this->measurements_[z] , landmark_pos );
	double newLandmarkWeight = config.landmarkExistencePrior_;
	maps_[i]->addGaussian( &landmark_pos, newLandmarkWeight, true);

      }
    }

    //---------- 4. Importance Weighting --------------
    // Some of the work has already been done in the map update step

    this->particleSet_[i]->setWeight( exp(logLikelihoodSum) );


    //---------- 5. Cleanup - Free memory ----------

    delete pose;
    for( int n = 0; n < nMZ; n++ ){
      delete[] likelihoodTable[n];
    }
    delete[] likelihoodTable;

  } // particle i loop end

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::resample(){

  bool resampleOccured = false;
  if( k_currentTimestep_ - k_lastResample_ >= config.minInterSampleTimesteps_){
    resampleOccured = this->resample();
  }
  if( resampleOccured ){
    k_lastResample_ = k_currentTimestep_;
  }else{
    this->normalizeWeights();
  }

  // At this point, the resampled particles do not have updated maps
  if( resampleOccured ){   // reassign maps as well according to resampling of particles

    std::vector< GaussianMixture<TLandmark>* > maps_temp( maps_.size(), NULL );
    std::vector< int > useCount ( maps_.size(), 0);

    // Note which GMs get used and how many times
    for(int i = 0; i < this->nParticles_; i++){
      int j = this->particleSet_[i]->getParentId();
      useCount[j]++;
    }

    // for maps that get used, make a (pointer) copy before doing any overwriting
    // Also rid of the maps that die along with particles
    for(int i = 0; i < this->nParticles_; i++){
      if( useCount[i] > 0){
	maps_temp[i] = maps_[i];
      }else{
	delete maps_[i];
      }
    }
    
    // Copy Gaussians
    for(int i = 0; i < this->nParticles_; i++){
      int j = this->particleSet_[i]->getParentId();

      if( useCount[j] == 1){ // map j is only copied over to map i and not to any other particle's map
	
	maps_[i] = maps_temp[j];
	maps_temp[j] = NULL;
	
      }else{ // GM_j is used by more than 1 particle, need to allocate memory for copying map
	
	maps_[i] = new GaussianMixture<TLandmark>;
	maps_temp[j]->copyTo( maps_[i] );
      }
      useCount[j]--;
    }

  }  

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
int FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getGMSize(int i){

  if( i >= 0 && i < maps_.size() )
    return ( maps_[i]->getGaussianCount() );
  else
    return -1;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
getLandmark(const int i, const int m, 
	    typename TLandmark::Vec &u,
	    typename TLandmark::Mat &S)
{

  int sz = getGMSize(i);
  if( sz == -1 || (m < 0) || (m >= sz) )
    {
      return false;
    }
    TLandmark *plm;
    maps_[i]->getGaussian(m, plm);
    plm->get(u, S);
    return true;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
KalmanFilter* FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getKalmanFilter(){
  return kfPtr_;
}

#endif
