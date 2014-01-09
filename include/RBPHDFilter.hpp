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

#ifndef RBPHDFILTER_HPP
#define RBPHDFILTER_HPP

#include <boost/timer/timer.hpp>
#include <Eigen/Core>
#include "GaussianMixture.hpp"
#include "ParticleFilter.hpp"
#include "KalmanFilter.hpp"
#include <math.h>
#include <vector>

#include <stdio.h>

/**
 *  \class RBPHDFilter
 *  \brief Rao-Blackwellized Probability Hypothesis Density Filter class
 *  
 *  This class implements the Rao-Bloackwellized Probability Hypothesis Density
 *  filter. The constructor of this class will internally instantiate the 
 *  process model for both the robot and landmarks, the measurement model, 
 *  and the Kalman filter. Users have access to these through pointers that 
 *  can be obtained by calling the appropraite get function.
 *
 *  \tparam RobotProcessModel A robot process model derived from ProcessModel
 *  \tparam LmkProcessModel A landmark process model derived from ProcessModel
 *  \tparam MeasurementModel A sensor model derived from MeasurementModel
 *  \tparam KalmanFilter A Kalman filter that uses LmkProcessModel and MeasurementModel
 *  \author Keith Leung, Felipe Inostroza
 */

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
class RBPHDFilter : public ParticleFilter<RobotProcessModel, MeasurementModel,
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
    
    /**  New birth Gaussians are set with this weight */
    double birthGaussianWeight_;   

    /**  New Gaussians are only created during map update if the innovation mahalanobis distance 
	 is less than this threshold */
    double newGaussianCreateInnovMDThreshold_;

    /**  number of map states to use for evaluating particle weight
	 0 => empty-set strategy,
	 1 => single-feature strategy,
	 >1 => multi-feature strategy
    */
    int importanceWeightingEvalPointCount_;

    /** The mahalanobis distance threshold used to determine if a possible meaurement-landmark
     *  pairing is significant to worth considering 
     */
    double importanceWeightingMeasurementLikelihoodMDThreshold_;

    /** Gaussian merging Mahalanobis distance threshold */
    double gaussianMergingThreshold_;

    /** Gaussian merging covariance inflation factor */
    double gaussianMergingCovarianceInflationFactor_;

    /** Gaussian pruning weight threshold, below which Gaussians are eliminated from a Gaussian mixture */
    double gaussianPruningThreshold_;

    /** Minimum number of updates betwen resampling of particles*/
    int minUpdatesBeforeResample_;

    /** If true, timing information is written to the console every update*/
    bool reportTimingInfo_;

    /** Use the particle weighting strategty from Single-cluster PHD Filtering by Lee, et. al. */
    bool useClusterProcess_;

  } config;

  /** 
   * Constructor 
   * \param n number of particles
   */
  RBPHDFilter(int n);

  /** Destructor */
  ~RBPHDFilter();

  /** 
   * Get the landmark process model
   * \return pointer to the landmark process model
   */
  LmkProcessModel* getLmkProcessModel();

  /**
   * Predict the robot trajectory using the lastest odometry data
   * \param[in] u input 
   * \param[in] currentTimestep current timestep (no longer used);
   * \param[in] useModelNoise use the additive noise for the process model
   * \param[in] useInputNoise use the noise fn the input
   */
  void predict( TInput u, TimeStamp const &dT,
		bool useModelNoise = true,
		bool useInputNoise = false);

  /**
   * Update the map, calculate importance weighting, sample if necessary, and
   * create new birth Gaussians.
   * \param[in] Z set of measurements to use for the update, placed in a std vector, which
   * gets cleared after the function call.
   */
  void update( std::vector<TMeasurement> &Z);

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
   * \param[out] w weight
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

  double **weightingTable_; /**< Weighting table used during map update */
  int weightingTableNRows_; /**< Number of rows in the weighting table */
  int weightingTableNCols_; /**< Number of cols in the weighting table */

  KalmanFilter *kfPtr_; /**< pointer to the Kalman filter */
  LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */

  /** indices of unused measurement for each particle for creating birth Gaussians */
  std::vector< std::vector<unsigned int> > unused_measurements_; 

  unsigned int nUpdatesSinceResample; /**< Number of updates performed since the last resmaple */
  
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
   * Importance weighting. Overrides the abstract function in ParticleFilter
   */
  void importanceWeighting();

  /**
   * Random Finite Set measurement likelihood evaluation
   * \brief The current measurements in measurements_ are used to determine the
   * RFS measurement likelihood given a set of landmarks 
   * \param[in] particleIdx particle for which the likelihood is calcuated
   * \param[in] evalPtIdx indices of evaluation points in maps_[particleIdx]
   * \param[in] evalPtPd probability of detection of evaluation point 
   * \return measurement likelihood
   */
  double rfsMeasurementLikelihood( const int particleIdx, 
				   std::vector<unsigned int> &evalPtIdx,
				   std::vector<double> &evalPtPd );

  /**
   * Calculate the sum of all permutations of measurement likelihood from a likelihood
   * table generated from within rfsMeasurementLikelihood
   * \param[in] likelihoodTab likelihood table generated within rfsMeasurementLikelihood
   * \param[in] Z_NoClutter A vector of measurement indices (columns) to consider in the likelihoodTab 
   * \return sum of all permutations from the given likelihood table
   */
  double rfsMeasurementLikelihoodPermutations( std::vector< double* > &likelihoodTab, 
					       std::vector< int > &Z_NoClutter);

  /** Checks the Gaussian mixture maps for all particles for errors
   *  \return true if there are no errors 
   */
  bool checkMapIntegrity();

  /** Ensure that sufficient memory is allocated for the weighting table
   *  \param[in] nRows desired number of rows
   *  \param[in] nCols desired number of columns
   */
  void checkWeightingTableSize(int nRows, int nCols);

};

////////// Implementation //////////

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::RBPHDFilter(int n)
: ParticleFilter<RobotProcessModel, MeasurementModel, 
		 GaussianMixture< typename MeasurementModel::TLandmark > >(n)
{

  lmkModelPtr_ = new LmkProcessModel;
  kfPtr_ = new KalmanFilter(lmkModelPtr_, this->getMeasurementModel());
  
  for(int i = 0; i < n; i++){
    printf("Creating map structure for particle %d\n", i);
    this->particleSet_[i]->setData( new GaussianMixture<TLandmark>() );
    unused_measurements_.push_back( std::vector<unsigned int>() );
  }
  
  config.birthGaussianWeight_ = 0.25; 
  config.gaussianMergingThreshold_ = 0.5;
  config.gaussianMergingCovarianceInflationFactor_ = 1.5;
  config.gaussianPruningThreshold_ = 0.2;
  config.importanceWeightingEvalPointCount_ = 8;
  config.importanceWeightingMeasurementLikelihoodMDThreshold_ = 3.0;
  config.newGaussianCreateInnovMDThreshold_ = 0.2;
  config.minUpdatesBeforeResample_ = 1;
  config.reportTimingInfo_ = false;
  
  nUpdatesSinceResample = 0;

  weightingTableNRows_ = 100; 
  weightingTableNCols_ = 100;
  weightingTable_ = new double*[weightingTableNRows_];
  for(int n = 0; n < weightingTableNRows_; n++){
    weightingTable_[n] = new double[weightingTableNRows_];
  }
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::~RBPHDFilter(){

  for(int i = 0; i < this->nParticles_; i++){
    this->particleSet_[i]->deleteData();
  }
  delete kfPtr_;
  delete lmkModelPtr_;

  for(int n = 0; n < weightingTableNRows_; n++ ){
    delete[] weightingTable_[n];
  }
  delete[] weightingTable_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
LmkProcessModel* RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getLmkProcessModel(){
  return lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::predict( TInput u,
												 TimeStamp const &dT,
												 bool useModelNoise,
												 bool useInputNoise){

  boost::timer::auto_cpu_timer *timer = NULL;
  if(config.reportTimingInfo_)
    timer = new boost::timer::auto_cpu_timer(6, "Predict time: %ws\n");


  // Add birth Gaussians using pose before prediction
  addBirthGaussians();

  // propagate particles
  this->propagate(u, dT);

  // propagate landmarks
  for( int i = 0; i < this->nParticles_; i++ ){
    for( int m = 0; m < this->particleSet_[i]->getData()->getGaussianCount(); m++){
      TLandmark *plm;
      this->particleSet_[i]->getData()->getGaussian(m, plm);
      lmkModelPtr_->staticStep(*plm, *plm, dT);
    }
  }

  if(timer != NULL)
    delete timer;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::update( std::vector<TMeasurement> &Z){

  boost::timer::auto_cpu_timer *timer_mapUpdate = NULL;
  boost::timer::auto_cpu_timer *timer_particleWeighting = NULL;
  boost::timer::auto_cpu_timer *timer_mapMerge = NULL;
  boost::timer::auto_cpu_timer *timer_mapPrune = NULL;
  boost::timer::auto_cpu_timer *timer_particleResample = NULL;

  nUpdatesSinceResample++;

  this->setMeasurements( Z ); // Z gets cleared after this call, measurements now stored in this->measurements_

  ////////// Map Update //////////
  if(config.reportTimingInfo_){
    timer_mapUpdate = new boost::timer::auto_cpu_timer(6, "Map update time: %ws\n");
  }
  updateMap();
  if(timer_mapUpdate != NULL)
    delete timer_mapUpdate;

  ////////// Particle Weighintg //////////
  if(config.reportTimingInfo_){
    timer_particleWeighting = new boost::timer::auto_cpu_timer(6, "Particle weighting time: %ws\n");
  }
  if(!config.useClusterProcess_){
    importanceWeighting();
  }
  if(timer_particleWeighting != NULL)
    delete timer_particleWeighting;

  //////////// Merge and prune //////////
  int maxMapSize = -1;
  int i_maxMapSize = -1;
  if(config.reportTimingInfo_){
    timer_mapMerge = new boost::timer::auto_cpu_timer(6, "Map merge time: %ws\n");
  }
  for( int i = 0; i < this->nParticles_; i++){ 

    this->particleSet_[i]->getData()->merge( config.gaussianMergingThreshold_, 
					     config.gaussianMergingCovarianceInflationFactor_);   
  } 
  if(timer_mapMerge != NULL)
    delete timer_mapMerge;
  
  if(config.reportTimingInfo_){
    timer_mapPrune = new boost::timer::auto_cpu_timer(6, "Map prune time: %ws\n");
  }
  for( int i = 0; i < this->nParticles_; i++){ 

    this->particleSet_[i]->getData()->prune( config.gaussianPruningThreshold_ );    

  }
  if(timer_mapPrune != NULL)
    delete timer_mapPrune;

  //////////// Particle resampling //////////
  if(config.reportTimingInfo_){
    timer_particleResample = new boost::timer::auto_cpu_timer(6, "Particle resample time: %ws\n");
  }
  bool resampleOccured = false;
  if( nUpdatesSinceResample >= config.minUpdatesBeforeResample_){
    resampleOccured = this->resample();
  }

  if( resampleOccured ){
    nUpdatesSinceResample = 0;
  }else{
    this->normalizeWeights();
  }
 
  if(timer_particleResample != NULL)
    delete timer_particleResample;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::updateMap(){

  const unsigned int startIdx = 0;
  const unsigned int stopIdx = this->nParticles_;
  const unsigned int nZ = this->measurements_.size();

  for(unsigned int i = startIdx; i < stopIdx; i++){    

    //---------- 1. setup / book-keeping ----------
   
    const unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
    unused_measurements_[i].clear();    
    if(nM == 0){ // No existing landmark case -> flag all measurements as unused and go to next particles
      for(int z = 0; z < nZ; z++){
	unused_measurements_[i].push_back( z );
      }
      continue; // goto next particle
    }
    double Pd[nM];
    int landmarkCloseToSensingLimit[nM];

    // For cluster process particle weighting
    double w_km_sum = std::numeric_limits<double>::denorm_min();
    double likelihoodProd = 1;
    if(config.useClusterProcess_){
      for(int m = 0; m < nM; m++){
	w_km_sum += this->particleSet_[i]->getData()->getWeight(m);
      }
    }

    // nM x nZ table for Gaussian weighting
    checkWeightingTableSize(nM, nZ);
    TLandmark*** newLandmarkPointer = new TLandmark** [ nM ];
    for( int n = 0; n < nM; n++ ){
      newLandmarkPointer[n] = new TLandmark* [ nZ ];
    }
    for(int m = 0; m < nM; m++){
      for(int z = 0; z < nZ; z++){
	newLandmarkPointer[m][z] = NULL;
	weightingTable_[m][z] = 0;
      }
    }

    //----------  2. Kalman Filter map update ----------

    const TPose *pose = this->particleSet_[i]->getPose();
    TLandmark* lmNew = NULL;

    for(unsigned int m = 0; m < nM; m++){

      TLandmark* lm = this->particleSet_[i]->getData()->getGaussian(m);
      bool isCloseToSensingLimit;
      Pd[m] = this->pMeasurementModel_->probabilityOfDetection( *pose, *lm, 
								isCloseToSensingLimit); 
      if(isCloseToSensingLimit){
	landmarkCloseToSensingLimit[m] = 1;
	Pd[m] = 1;
      }else{
	landmarkCloseToSensingLimit[m] = 0;
      }
      double w_km = this->particleSet_[i]->getData()->getWeight(m);
      double Pd_times_w_km = Pd[m] * w_km;

      if(Pd[m] != 0){
	for(int z = 0; z < nZ; z++){

	  if(lmNew == NULL)
	    lmNew = new TLandmark;

	  newLandmarkPointer[m][z] = NULL;
	  weightingTable_[m][z] = 0;
	  double innovationLikelihood = 0;
	  double innovationMahalanobisDist2 = 0;
	  double threshold = config.newGaussianCreateInnovMDThreshold_ * config.newGaussianCreateInnovMDThreshold_;
	
	  // RUN KF, create new landmark for likely updates but do not add to map_[i] yet
	  // because we cannot determine actual weight until the entire weighting table is
	  // filled in
	  bool updateMade = kfPtr_->correct(*pose, this->measurements_[z], *lm, *lmNew, 
					    &innovationLikelihood, &innovationMahalanobisDist2);

	  if ( !updateMade || innovationMahalanobisDist2 > threshold ){
	    newLandmarkPointer[m][z] = NULL;
	    weightingTable_[m][z] = 0;
	  }else{
	    newLandmarkPointer[m][z] = lmNew;
	    lmNew = NULL;
	    weightingTable_[m][z] = Pd_times_w_km * innovationLikelihood;
	  }	

	} // z forloop end

      }else{ // Pd = 0
	for(int z = 0; z < nZ; z++){
	  newLandmarkPointer[m][z] = NULL;
	  weightingTable_[m][z] = 0;
	}
      }

    } // m forloop end

    if(lmNew != NULL)
      delete lmNew;
    
    // Now calculate the weight of each new Gaussian
    for(int z = 0; z < nZ; z++){

      double clutter = this->pMeasurementModel_->clutterIntensity( this->measurements_[z], nZ );
      double sum = clutter;
      for(unsigned int m = 0; m < nM; m++){
	sum += weightingTable_[m][z];
      }

      if(config.useClusterProcess_){
	likelihoodProd *= sum;
      }

      for(unsigned int m = 0; m < nM; m++){
	weightingTable_[m][z] = weightingTable_[m][z] / sum;
      }
    }
    if(config.useClusterProcess_){
      double prev_particle_i_weight = this->particleSet_[i]->getWeight();
      this->particleSet_[i]->setWeight( exp(w_km_sum) * likelihoodProd);
    }

    // ---------- 3. Add new Gaussians to map  ----------
    // New Gaussians will have indices >= nM 
    for(int m = 0; m < nM; m++){
      for(int z = 0; z < nZ; z++){
	if(newLandmarkPointer[m][z] != NULL && weightingTable_[m][z] > 0){
	  this->particleSet_[i]->getData()->addGaussian( newLandmarkPointer[m][z],  weightingTable_[m][z]);  
	}
      }
    }

    //----------  4. Determine weights for existing Gaussians (missed detection) ----------
    for(int m = 0; m < nM; m++){
      
      double w_km = this->particleSet_[i]->getData()->getWeight(m);
      double w_k = (1 - Pd[m]) * w_km;

      // For landmarks close to sensing limit
      if (landmarkCloseToSensingLimit[m] == 1 && w_km > config.birthGaussianWeight_){
	double weight_sum_m = 0;
	for(int z = 0; z < nZ; z++){
	  weight_sum_m += weightingTable_[m][z];
	}
	double delta_w = Pd[m] * w_km - weight_sum_m;
	if( delta_w > 0 ){
	  w_k += delta_w; // This is just a heuristic that works well
	  if(w_k > 1)
	    w_k = 1;
	}
      }

      this->particleSet_[i]->getData()->setWeight(m, w_k);
    }

    //----------  5. Identify unused measurements for adding birth Gaussians later ----------
    unused_measurements_[i].clear();
    for(int z = 0; z < nZ; z++){
      int useCount = 0;
      for(int m = 0; m < nM; m++){
	if (weightingTable_[m][z] != 0){
	  useCount++;
	}
      }
      if (useCount == 0)
	unused_measurements_[i].push_back( z );
    }

    //----------  6. Cleanup - Free memory ----------

    for( int n = 0; n < nM; n++ ){
      delete[] newLandmarkPointer[n];
    }
    delete[] newLandmarkPointer;

  }

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::importanceWeighting(){

  for(int i = 0; i < this->nParticles_; i++){

    //printf("Importance weighting for particle %d\n", i);
    TPose x;
    this->particleSet_[i]->getPose( x );

    // 1. select evaluation points from highest-weighted Gaussians after update, that are within sensor FOV
    const unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
    int nEvalPoints = config.importanceWeightingEvalPointCount_ > nM ? nM : config.importanceWeightingEvalPointCount_ ;
    std::vector<unsigned int> evalPointIdx;
    std::vector<double> evalPointPd;
    evalPointIdx.reserve(nEvalPoints);
    evalPointPd.reserve(nEvalPoints);
    if( nEvalPoints == 0 ){
      this->particleSet_[i]->setWeight( std::numeric_limits<double>::denorm_min() );
      continue;
    }
    this->particleSet_[i]->getData()->sortByWeight(); // sort by weight so that we can pick off the top nEvalPoints Gaussians
    for(int m = 0; m < nM; m++){
      TLandmark* plm_temp;
      double w, w_prev;
      bool closeToSensingLim;
      this->particleSet_[i]->getData()->getGaussian(m, plm_temp, w, w_prev);
      double Pd = this->pMeasurementModel_->probabilityOfDetection(x,*plm_temp, closeToSensingLim);
      if( Pd > 0 ){
	evalPointIdx.push_back(m);
	evalPointPd.push_back(Pd);
      }
      if(evalPointIdx.size() >= nEvalPoints)
	break;
    }
    nEvalPoints = evalPointIdx.size();

    // 2. evaluate sum of Gaussian weights
    double gaussianWeightSumBeforeUpdate = 0;
    double gaussianWeightSumAfterUpdate = 0;
    for(int m = 0; m < nM; m++){
      TLandmark* plm_temp;
      double w, w_prev;
      this->particleSet_[i]->getData()->getGaussian(m, plm_temp, w, w_prev); // for newly created Gaussians, w_prev = 0
      gaussianWeightSumBeforeUpdate += w_prev;
      gaussianWeightSumAfterUpdate += w;
      //printf("w_[%d] = %f\n", i, w);
    }
    // Check for NaN
    if( gaussianWeightSumBeforeUpdate != gaussianWeightSumBeforeUpdate ||
	gaussianWeightSumAfterUpdate != gaussianWeightSumAfterUpdate ){
      //printf("Particle %d map size before update = %f\n", i, gaussianWeightSumBeforeUpdate);
      //printf("Particle %d map size after update = %f\n", i, gaussianWeightSumAfterUpdate);
    }
    //printf("Particle %d map size after - before update = %f\n", i, gaussianWeightSumAfterUpdate - gaussianWeightSumBeforeUpdate);
    
    // 3. evaluate intensity function at eval points and take their product
    double intensityProd_beforeUpdate = 1;
    double intensityProd_afterUpdate = 1;
    for(int e = 0; e < nEvalPoints; e++){

      int p = evalPointIdx[e];
      TLandmark* lm_evalPt;
      double w_temp;
      this->particleSet_[i]->getData()->getGaussian(p, lm_evalPt, w_temp);

      double intensity_at_evalPt_beforeUpdate = std::numeric_limits<double>::denorm_min();
      double intensity_at_evalPt_afterUpdate = std::numeric_limits<double>::denorm_min();

      for(int m = 0; m < nM; m++){
	TLandmark* plm;
	double w, w_prev;
	this->particleSet_[i]->getData()->getGaussian(m, plm, w, w_prev);
	// New Gaussians from update will have w_prev = 0
	// Out Gaussians (missed-detection) will not have been updated, but weights will have changed
	double likelihood = plm->evalGaussianLikelihood( *lm_evalPt );
	intensity_at_evalPt_beforeUpdate += w_prev * likelihood; // w_prev for newly created Gaussians are 0
	intensity_at_evalPt_afterUpdate += w * likelihood;
	// NaN check
	if( likelihood != likelihood || 
	    intensity_at_evalPt_beforeUpdate != intensity_at_evalPt_beforeUpdate || 
	    intensity_at_evalPt_afterUpdate != intensity_at_evalPt_afterUpdate){
	  printf("Particle %d map intensity error for eval point %d\n", i, m);
	  printf("intensity before update = %f\n", intensity_at_evalPt_beforeUpdate);
	  printf("intensity after update = %f\n", intensity_at_evalPt_afterUpdate);
	}
      }
      intensityProd_beforeUpdate *= intensity_at_evalPt_beforeUpdate;
      intensityProd_afterUpdate *= intensity_at_evalPt_afterUpdate;
      // NaN Check
      if( intensityProd_beforeUpdate != intensityProd_beforeUpdate ||
	  intensityProd_afterUpdate != intensityProd_afterUpdate ){
	printf("Particle %d map intensity product error\n", i);
	printf("intensity product before update = %f\n", intensityProd_beforeUpdate);
	printf("intensity product after update = %f\n", intensityProd_afterUpdate);
      } 
    }
    //printf("Particle %d intensity product before / after update = %f\n", i, intensityProd_beforeUpdate / intensityProd_afterUpdate);

    // 4. calculate measurement likelihood at eval points
    // note that rfsMeasurementLikelihood uses maps_[i] which is already sorted by weight
    double measurementLikelihood = rfsMeasurementLikelihood( i, evalPointIdx, evalPointPd );
    //printf("Particle %d measurement likelihood = %f\n", i, measurementLikelihood);

    // 5. calculate overall weight
    double overall_weight = measurementLikelihood * intensityProd_beforeUpdate / intensityProd_afterUpdate *
      exp( gaussianWeightSumAfterUpdate - gaussianWeightSumBeforeUpdate); 
    
    double prev_weight = this->particleSet_[i]->getWeight();
    this->particleSet_[i]->setWeight( overall_weight * prev_weight );
    //printf("Particle %d overall weight = %f\n\n", i, overall_weight);

  }

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
double RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
rfsMeasurementLikelihood( const int particleIdx, 
			  std::vector<unsigned int> &evalPtIdx,
			  std::vector<double> &evalPtPd ){
  // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight; 

  const int i = particleIdx;
  const int nM = evalPtIdx.size();
  const int nZ = this->measurements_.size();

  // Fill in likelihood table
  TPose x;
  this->particleSet_[i]->getPose( x );
  std::vector< double* > likelihoodTab;
  likelihoodTab.reserve(nM);
  for( int m = 0; m < nM; m++ ){
    
    double* row = new double[nZ];
    
    TLandmark* evalPt; 
    this->particleSet_[i]->getData()->getGaussian( evalPtIdx[m], evalPt );
    double Pd = evalPtPd[m];
    
    double row_sum = 0;
    for( int z = 0; z < nZ; z++ ){

      TMeasurement expected_z;
      this->pMeasurementModel_->measure( x, *evalPt, expected_z);
      //TMeasurement actual_z = this->measurements_[z];
      
      double md2;
      double threshold = config.importanceWeightingMeasurementLikelihoodMDThreshold_;
      threshold *= threshold;
      double likelihood = this->measurements_[z].evalGaussianLikelihood( expected_z, &md2);
      if( md2 <= threshold ){
	row[z] = likelihood * Pd;
	row_sum += row[z];
      }else{
	row[z] = 0;
      }
      
    }
    
    // Check to see if likelihood to all measurements is 0, if so remove this eval point
    // as it will not contribute anything to the output likelihood
    if( row_sum == 0 ){
      delete[] row;
    }else{
      likelihoodTab.push_back(row);
    }
    
  }
  const int likelihoodTabSizeWithoutClutter = likelihoodTab.size();

  // Check measurements (columns) with 0 likelihood to all eval points
  // if so, that measurement is considered clutter
  int nClutter = 0;
  double clutterLikelihood = 1; // we will multiply the likelihood sum of all d.a. permuations with this at the end
  std::vector<int> z_noClutter; // we only want the non-clutter measurements when we permutate over all data assocation pairs later
  z_noClutter.reserve(nZ);

  for( int z = 0; z < nZ; z++ ){
    double isClutter = true;
    for( int m = 0; m < likelihoodTab.size(); m++ ){
      if( likelihoodTab[m][z] > 0 ){
	isClutter = false;
	break;
      }
    }
    if( isClutter ){
      nClutter++;
      // TMeasurement actual_z = this->measurements_[z];
      clutterLikelihood *= this->pMeasurementModel_->clutterIntensity(this->measurements_[z], nZ);;
    }else{
      z_noClutter.push_back(z);
    }
  }

  // If the number of measurements is greater than the number of
  // eval points, then some measurements have to be assigned as clutter.
  // We will add extra rows for these assignments in likelihoodTab
  double *clutterRow = NULL;
  int nR = z_noClutter.size() - likelihoodTab.size();
  if (nR > 0){
    clutterRow = new double[nZ];
    for( int z = 0; z < nZ; z++ ){
      //TMeasurement actual_z = this->measurements_[z];
      clutterRow[z] = this->pMeasurementModel_->clutterIntensity(this->measurements_[z], nZ);
    }
    for( int r = 0; r < nR; r++ ){
      likelihoodTab.push_back(clutterRow);
    }
  }

  // Go through all permutations of eval point - measurement pairs
  // to calculate the likelihood
  double likelihood = 0;
  while (likelihood == 0){

    if( likelihoodTab.size() == 0 ){
      likelihood = 1;
      break;
    }

    likelihood = rfsMeasurementLikelihoodPermutations( likelihoodTab, z_noClutter);
    if( likelihood != likelihood ){
      printf("RFS Measurement likelihood = %f", likelihood);
    } 

    if( likelihood == 0 ){

      // Add another row to of clutter to likelihoodTab
      // printf("Adding clutter row for RFS measurement likelihood calculation\n");
      if( clutterRow == NULL ){
	clutterRow = new double[nZ];
	for( int z = 0; z < nZ; z++ ){
	  //TMeasurement actual_z = this->measurements_[z];
	  clutterRow[z] = this->pMeasurementModel_->clutterIntensity(this->measurements_[z], nZ);
	}
      }
      likelihoodTab.push_back( clutterRow );
      
    }

  }

  // Deallocate likelihood table
  for( int m = 0; m < likelihoodTabSizeWithoutClutter; m++ ){
    delete[] likelihoodTab[m];
  }
  if( clutterRow != NULL )
    delete[] clutterRow;

  if (nClutter > 0){
    likelihood /= this->pMeasurementModel_->clutterIntensityIntegral( nZ );
    likelihood *= clutterLikelihood;
  }

  return likelihood;
}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
double RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
rfsMeasurementLikelihoodPermutations( std::vector< double* > &likelihoodTab, 
				      std::vector< int > &Z_NoClutter){
  // Note that nM is always >= nZ
  // We will find all eval point permutations of (0, 1, 2, ... , nM - 1)
  // and use the first nZ of each permutation to calculate the likelihood

  // A function required for sorting
  struct sort{
    static bool descend(int i, int j){ return (i > j); }
  };

  const int nM = likelihoodTab.size();
  const int nZ = Z_NoClutter.size();
  double allPermutationLikelihood = 0;
  bool lastPermutationSequence = false;
  std::vector<int> currentPermutation(nM);
  for(int m = 0; m < nM; m++){
    currentPermutation[m] = m;
  }

  while( !lastPermutationSequence ){

    // find the likelihood of the current permutation
    
    double currentPermutationLikelihood  = 1;
    for(int z = 0; z < nZ; z++){
      int m = currentPermutation[ z ];
      currentPermutationLikelihood *= likelihoodTab[m][ Z_NoClutter[z] ];
      
      // Fast-forward permutation if we know that following sequences will also
      // have 0 likelihood
      if( currentPermutationLikelihood == 0 && z < nZ - 1){
	std::sort(currentPermutation.begin() + z + 1, currentPermutation.end(), sort::descend);
	break;
      }
    }

    allPermutationLikelihood += currentPermutationLikelihood;

    // Fast-forward if nM > nZ (i.e., the last nM - nZ elements in the permutation sequence does not matter)
    if( nM > nZ ){
      std::sort(currentPermutation.begin() + nZ, currentPermutation.end(), sort::descend);
    }

    // Generate the next permutation sequence
    for(int m = nM - 2; m >= -1; m--){

      if( m == -1){
	lastPermutationSequence = true;
	break;
      }
      
      // Find the highest index m such that currentPermutation[m] < currentPermutation[m+1]
      if(currentPermutation[m] < currentPermutation[ m + 1 ]){

	// Find highest index i such that currentPermutation[i] > currentPermutation[m] 
	// then swap the elements
	for(int i = nM - 1; i >= 0; i--){
	  if( currentPermutation[i] > currentPermutation[m] ){
	    int temp = currentPermutation[i];
	    currentPermutation[i] = currentPermutation[m];
	    currentPermutation[m] = temp;
	    break;
	  }
	}

	// reverse order of elements after currentPermutation[m]
	int nElementsToSwap = nM - (m + 1);
	int elementsToSwapMidPt = nElementsToSwap / 2;
	int idx1 = m + 1;
	int idx2 = nM - 1;
	for(int i = 1; i <= elementsToSwapMidPt; i++){
	  int temp = currentPermutation[idx1];
	  currentPermutation[idx1] = currentPermutation[idx2];
	  currentPermutation[idx2] = temp;
	  idx1++;
	  idx2--;
	}

	break;
      }

    }

    // now we should have the next permutation sequence
  }

  return allPermutationLikelihood;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::addBirthGaussians(){

  for(int i = 0; i < this->nParticles_; i++){
    
    while( unused_measurements_[i].size() > 0){
     
      // get measurement
      int unused_idx = unused_measurements_[i].back();
      TMeasurement unused_z = this->measurements_[unused_idx];
      unused_measurements_[i].pop_back();

      // use inverse measurement model to get landmark
      TPose robot_pose;
      TLandmark landmark_pos;
      this->particleSet_[i]->getPose(robot_pose);
      this->pMeasurementModel_->inverseMeasure( robot_pose,  unused_z, landmark_pos );
      
      // add birth landmark to Gaussian mixture (last param = true to allocate mem)
      this->particleSet_[i]->getData()->addGaussian( &landmark_pos, config.birthGaussianWeight_, true);
      
    }
    
  }

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::checkMapIntegrity(){

  for( int i = 0; i < this->nParticles_; i++ ){

    unsigned int nM = this->particleSet_[i]->getData()->getGaussianCount();
    for(int m = 0; m < nM; m++){

      TLandmark *lm;
      typename TLandmark::Vec lm_x;
      typename TLandmark::Mat lm_S;
      double w;

      this->particleSet_[i]->getData()->getGaussian(m, lm, w);
      if( lm != NULL ){
	lm->get(lm_x, lm_S);

	bool vecError = false;
	for( int c = 0; c < lm_x.rows(); c++){
	  if(lm_x(c) != lm_x(c)){
	    vecError = true;
	    break;
	  }
	}
	if(vecError){
	  printf("particle %d, landmark index %d, vector error\n", i, m);
	  std::cout << lm_x << std::endl;
	  return false;
	}

	bool matError = false;
	for( int r = 0; r < lm_S.rows(); r++){
	  for( int c = 0; c < lm_S.cols(); c++){
	    if(lm_S(r,c) != lm_S(r,c)){
	      matError = true;
	      break;
	    }
	  }
	}
	lm_x.setOnes();
	double posDefCheck = lm_x.transpose() * lm_S * lm_x;
	if( posDefCheck != posDefCheck || posDefCheck <= 0){
	  matError = true;
	}
	if(matError){
	  printf("particle %d, landmark index %d, covariance error\n", i, m);
	  std::cout << lm_S << std::endl;
	  return false;
	}

	if( w != w ){
	  printf("particle %d, landmark index %d, w = %f\n", i, m, w);
	  return false;
	}

      }

    }

  }

  return true;

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
int RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getGMSize(int i){

  if( i >= 0 && i < this->nParticles_ )
    return ( this->particleSet_[i]->getData()->getGaussianCount() );
  else
    return -1;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
bool RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
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
    this->particleSet_[i]->getData()->getGaussian(m, plm, w);
    plm->get(u, S);
    return true;
}

#endif


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
setParticlePose(int i, TPose &p){
  
  this->particleSet_[i]->setPose(p);

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
KalmanFilter* RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getKalmanFilter(){
  return kfPtr_;
}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::checkWeightingTableSize(int nRows, int nCols){

  nRows *= 1.2;
  nCols *= 1.2;
  if( weightingTableNRows_ < nRows ){ // Row and Col need to increase
    if( weightingTableNCols_ < nCols ){
      for(int m = 0; m < weightingTableNRows_; m++ ){
	delete[] weightingTable_[m];
      }
      delete[] weightingTable_;
      weightingTable_ = new double* [nRows];
      for(int m = 0; m < nRows ; m++ ){
	weightingTable_[m] = new double [nCols];
      }
      weightingTableNRows_ = nRows;
      weightingTableNCols_ = nCols; 
    }else{ // Only increase row
      double** weightingTableOld = weightingTable_;
      weightingTable_ = new double* [nRows];
      for(int m = 0; m < nRows ; m++ ){
	if( m < weightingTableNRows_){
	  weightingTable_[m] = weightingTableOld[m];
	}else{
	  weightingTable_[m] = new double [weightingTableNCols_];
	}
      }
      weightingTableNRows_ = nRows;
    }
      }else if(weightingTableNCols_ < nCols){ // Only increase Col
    for(int m = 0; m < weightingTableNRows_; m++ ){
      delete[] weightingTable_[m];
      weightingTable_[m] = new double[nCols];
    }
    weightingTableNCols_ = nCols;
  }

}
