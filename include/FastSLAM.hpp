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
class FastSLAM : public ParticleFilter<RobotProcessModel, MeasurementModel, 
				       GaussianMixture< typename MeasurementModel::TLandmark > >
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

    /** Minimum timeteps betwen resampling of particles*/
    double minInterSampleTimesteps_;

    /** If true, timing information is written to the console every update*/
    bool reportTimingInfo_;

    /** The probability of existence that is initially assigned to a new landmark */
    double landmarkExistencePrior_;

    /** The log odds threshold for eliminating a landmark from the map*/
    double mapExistencePruneThreshold_;

    /** Minimum log measurement likelihood for numerical stability */
    double minLogMeasurementLikelihood_;


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
   * \param[out] w log-odds of existence
   * \return false if the indices specified are invalid 
   */ 
  bool getLandmark(const int i, const int m, 
		   typename TLandmark::Vec &u,
		   typename TLandmark::Mat &S,
		   double &w);
  
  /**
   * Get the pointer to the Kalman Filter used for updating the map
   * \return pointer to the Kalman Filter
   */
  KalmanFilter* getKalmanFilter();

  /** Function for initiating particles during startup 
   *  \param[in] i particle index
   *  \param[in] p particle pose
   */
  void setParticlePose(int i, TPose &p);


private:

  KalmanFilter *kfPtr_; /**< pointer to the Kalman filter */
  LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */

  // std::vector< GaussianMixture<TLandmark>* > maps_; /**< Particle dependent maps */

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
   * /return whether the update was successful
   */
  bool updateMap();

  /**
   * Resample the particles, along with their individual maps,  according to their 
   * importance weightings.
   */
  void resampleWithMapCopy();

  /** 
   * Importance weighting. Overrides the abstract function in ParticleFilter
   * \note For this FastSLAM algorithm, we perform weighting as part of 
   * mapUpdate() to be a little more efficient. Therefore, this function is 
   * not called at all. However, we still need to overwrite the virtual function.
   */
  void importanceWeighting(){}

};

////////// Implementation //////////

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::FastSLAM(int n)
: ParticleFilter<RobotProcessModel, MeasurementModel, 
		 GaussianMixture< typename MeasurementModel::TLandmark > >(n)
{

  lmkModelPtr_ = new LmkProcessModel;
  kfPtr_ = new KalmanFilter(lmkModelPtr_, this->getMeasurementModel());
  
  for(int i = 0; i < n; i++){
    // printf("Creating map structure for particle %d\n", i);
    //maps_.push_back( new GaussianMixture<TLandmark>() );
    this->particleSet_[i]->setData( new GaussianMixture<TLandmark>() );
  }
  
  config.minInterSampleTimesteps_ = 5;
  config.reportTimingInfo_ = false;
  config.landmarkExistencePrior_ = 0.5;
  config.mapExistencePruneThreshold_ = -3.0;
  config.minLogMeasurementLikelihood_ = -10.0;
  
  k_lastResample_ = -10;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::~FastSLAM(){
  for(int i = 0; i < this->nParticles_; i++){
    //delete maps_[i];
    this->particleSet_[i]->deleteData();
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
    //for( int m = 0; m < maps_[i]->getGaussianCount(); m++){
    for( int m = 0; m < this->particleSet_[i]->getData()->getGaussianCount(); m++){
      TLandmark *plm;
      //maps_[i]->getGaussian(m, plm);
      this->particleSet_[i]->getData()->getGaussian(m, plm);
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
  if(!updateMap())
    printf("k = %d     Update Failed!\n", k_currentTimestep_);
  if(timer_mapUpdate != NULL)
    delete timer_mapUpdate;

  //////////// Particle resampling //////////
  if(config.reportTimingInfo_){
    timer_particleResample = new boost::timer::auto_cpu_timer(6, "Particle resample time: %ws\n");
  }
  resampleWithMapCopy();
  if(timer_particleResample != NULL)
    delete timer_particleResample;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::updateMap(){

  const unsigned int startIdx = 0;
  const unsigned int stopIdx = this->nParticles_;

  const unsigned int nZ = this->measurements_.size();

  HungarianMethod hm; // for data association

  for(unsigned int i = startIdx; i < stopIdx; i++){    

    //----------  1. Data Association --------------------
    
    // Look for landmarks within sensor range
    const TPose *pose = this->particleSet_[i]->getPose();
    //unsigned int nM = maps_[i]->getGaussianCount();
    unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
    std::vector<int> idx_inRange;
    std::vector<double> pd_inRange;
    std::vector<TLandmark*> lm_inRange;
    for( int m = 0; m < nM; m++ ){
      //TLandmark* lm = maps_[i]->getGaussian(m);
      TLandmark* lm = this->particleSet_[i]->getData()->getGaussian(m);
      bool closeToLimit = false;
      double pd = this->pMeasurementModel_->probabilityOfDetection(*pose, *lm, closeToLimit); 
      if( pd != 0 || closeToLimit ){
	idx_inRange.push_back(m);
	lm_inRange.push_back(lm);
	pd_inRange.push_back(pd);
      }
    }
    nM = lm_inRange.size();

    // Loglikelihood table for Data Association
    unsigned int nMZ = nM;
    if(nZ > nM){
      nMZ = nZ;
    }
    double** likelihoodTable = new double* [ nMZ ];
    for( int m = 0; m < nMZ; m++ ){
      likelihoodTable[m] = new double [ nMZ ];
      for(int z = 0; z < nMZ; z++){
	likelihoodTable[m][z] = config.minLogMeasurementLikelihood_;
      }
    }

    // Fill in table for landmarks within range
    for(unsigned int m = 0; m < nM; m++){

      TLandmark* lm = lm_inRange[m]; // landmark position estimate
      TMeasurement measurement_exp; // expected measurement      
      bool isValidExpectedMeasurement = this->pMeasurementModel_->measure( *pose , *lm , measurement_exp);
      for(int z = 0; z < nZ; z++){
	if( isValidExpectedMeasurement ){
	  likelihoodTable[m][z] = fmax(config.minLogMeasurementLikelihood_, 
				       log(measurement_exp.evalGaussianLikelihood(this->measurements_[z])));
	}
      }
    }

    // Use Hungaian Method for data association
    int da[nMZ]; // data association result
    double logLikelihoodSum = 0; // used for particle weighting
    if(!hm.run(likelihoodTable, nMZ, da, &logLikelihoodSum)){
      printf("Update failed\n");
      hm.run(likelihoodTable, nMZ, da, &logLikelihoodSum, true, true); // true flag at the end is for debug
      return false;
    }

    //----------  2. Kalman Filter map update ----------

    bool zUsed[nZ];
    for(unsigned int z = 0; z < nZ; z++){
      zUsed[z] = false;
    }
    
    double nExpectedClutter = this->pMeasurementModel_->clutterIntensityIntegral(nZ);
    double probFalseAlarm = nExpectedClutter / nZ;
    double p_exist_given_Z = 0;
    double logParticleWeight = 0;

    for(unsigned int m = 0; m < nM; m++){
      
      TLandmark* lm = lm_inRange[m]; // landmark position estimate before update

      // Update landmark estimate m with the associated measurement
      int z = da[m];
      bool isUpdatePerformed = false;
      if(z < nZ){
	isUpdatePerformed = kfPtr_->correct(*pose, this->measurements_[z], *lm, *lm);
      }
      
      // calculate change to existence probability      
      if(isUpdatePerformed){

	zUsed[z] = true; // This flag is for new landmark creation	
	logParticleWeight += likelihoodTable[m][z];
	p_exist_given_Z = ((1 - pd_inRange[m]) * probFalseAlarm * config.landmarkExistencePrior_ + pd_inRange[m] * config.landmarkExistencePrior_) /
	  (probFalseAlarm + (1 - probFalseAlarm) * pd_inRange[m] * config.landmarkExistencePrior_); 

      }else{ // landmark estimate m not updated

	p_exist_given_Z = ((1 - pd_inRange[m]) * config.landmarkExistencePrior_) /
	  ((1 - config.landmarkExistencePrior_) + (1 - pd_inRange[m]) * config.landmarkExistencePrior_);

      }

      //double w = maps_[i]->getWeight( idx_inRange[m] ); // this is the log-odds of existence given previous measurements
      double w = this->particleSet_[i]->getData()->getWeight( idx_inRange[m] );
      w += log( (p_exist_given_Z) / (1 - p_exist_given_Z) ); 
      //maps_[i]->setWeight(idx_inRange[m], w);
      this->particleSet_[i]->getData()->setWeight(idx_inRange[m], w);
    }

    //---------- 3. Map Management (Add and remove landmarks)  ------------
    //int nRemoved = maps_[i]->prune(config.mapExistencePruneThreshold_); 
    int nRemoved = this->particleSet_[i]->getData()->prune(config.mapExistencePruneThreshold_); 

    for(unsigned int z = 0; z < nZ; z++){
      if(!zUsed[z]){ // Create new landmarks with inverse measurement model with unused measurements

	TLandmark landmark_pos;
	this->pMeasurementModel_->inverseMeasure( *pose, this->measurements_[z] , landmark_pos );
	double newLandmarkWeight = log(config.landmarkExistencePrior_ / (1 - config.landmarkExistencePrior_));
	//maps_[i]->addGaussian( &landmark_pos, newLandmarkWeight, true);
	this->particleSet_[i]->getData()->addGaussian( &landmark_pos, newLandmarkWeight, true);
      }
    }

    //---------- 4. Importance Weighting --------------
    // Some of the work has already been done in the map update step
    this->particleSet_[i]->setWeight( exp(logParticleWeight) );


    //---------- 5. Cleanup - Free memory ---------
    for( int n = 0; n < nMZ; n++ ){
      delete[] likelihoodTable[n];
    }
    delete[] likelihoodTable;

  } // particle i loop end
  return true;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::resampleWithMapCopy(){

  bool resampleOccured = false;

  if( k_currentTimestep_ - k_lastResample_ >= config.minInterSampleTimesteps_){
    resampleOccured = this->resample();
  }
  if( resampleOccured ){
    k_lastResample_ = k_currentTimestep_;
  }else{
    this->normalizeWeights();
  }

  /*
  // At this point, the resampled particles do not have updated maps
  if( resampleOccured ){   // reassign maps as well according to resampling of particles

    std::vector< GaussianMixture<TLandmark>* > maps_temp( maps_.size(), NULL );
    std::vector< int > useCount ( maps_.size(), 0); // maps_.size() = size of particle set before resamling

    // Note which GMs get used and how many times
    for(int i = 0; i < this->nParticles_; i++){
      int j = this->particleSet_[i]->getParentId();
      useCount[j]++;
    }

    // for maps that get used, make a (pointer) copy before doing any overwriting
    // Also rid of the maps that die along with particles
    for(int i = 0; i < maps_.size(); i++){
      if( useCount[i] > 0){
	maps_temp[i] = maps_[i];
      }else{
	delete maps_[i];
      }
    }
    // At this point all unused maps have been deleted
    
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
  */

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
int FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getGMSize(int i){

  if( i >= 0 && i < this->nParticles_ )
    //return ( maps_[i]->getGaussianCount() );
    return ( this->particleSet_[i]->getData()->getGaussianCount() );
  else
    return -1;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
getLandmark(const int i, const int m, 
	    typename TLandmark::Vec &u,
	    typename TLandmark::Mat &S,
	    double &w)
{

  int sz = getGMSize(i);
  if( sz == -1 || (m < 0) || (m >= sz) )
    {
      return false;
    }
  TLandmark *plm;
  //maps_[i]->getGaussian(m, plm, w);
  this->particleSet_[i]->getData()->getGaussian(m, plm, w);
  plm->get(u, S);
  return true;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
setParticlePose(int i, TPose &p){
  
  this->particleSet_[i]->setPose(p);

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
KalmanFilter* FastSLAM< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getKalmanFilter(){
  return kfPtr_;
}

#endif
