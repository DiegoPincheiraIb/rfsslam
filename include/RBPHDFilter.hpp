// RB-PHD-Filter Class
// Keith Leung 2013

#ifndef RBPHDFILTER_HPP
#define RBPHDFILTER_HPP

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
 *  \version 0.1 
 */

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
class RBPHDFilter : public ParticleFilter<RobotProcessModel, MeasurementModel>
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef typename RobotProcessModel::TState TPose;
  typedef typename RobotProcessModel::TInput TInput;
  typedef typename MeasurementModel::TLandmark TLandmark;
  typedef typename MeasurementModel::TMeasurement TMeasurement;
  typedef typename GaussianMixture<TLandmark>::Gaussian TGaussian;

  /** 
   * \brief Configurations for running the RB-PHD-Filter 
   */
  struct Config{
    
    /**  new birth Gaussians are set with this weight */
    double birthGaussianWeight_;   

    /**  new Gaussians are only created during map update if the innovation mahalanobis distance 
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

    /** Gaussian pruning weight threshold */
    double gaussianPruningThreshold_;

    /** Minimum timeteps betwen resampling */
    double minInterSampleTimesteps_;


  } config;

  /** 
   * Constructor 
   * \param n number of particles
   * \param initState initial state of particles
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
   * \param[in] input 
   * \param[in] currentTimestep current timestep;
   */
  void predict( TInput u, int currentTimestep );

  /**
   * Update the map, calculate importance weighting, sample if necessary, and
   * create new birth Gaussians.
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
   * \param[out] w weight
   * \return false if the indices specified are invalid 
   */ 
  bool getLandmark(const int i, const int m, 
		   typename TLandmark::Vec &u,
		   typename TLandmark::Mat &S,
		   double &w);
  
  void setParticlePose(int i, TPose &p);

private:

  KalmanFilter *kfPtr_; /**< pointer to the Kalman filter */
  LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */

  std::vector< GaussianMixture<TLandmark>* > maps_; /**< Particle dependent maps */

  /** indices of unused measurement for each particle for creating birth Gaussians */
  std::vector< std::vector<unsigned int> > unused_measurements_; 

  int k_currentTimestep_; /** current time */
  int k_lastResample_; /** last resample time *
  
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
   * \param[in] maxNumOfEvalPoints the maximum number of map evaluation points to use
   * \return measurement likelihood
   */
  double rfsMeasurementLikelihood( const int particleIdx, const int maxNumOfEvalPoints );

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
RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::RBPHDFilter(int n)
  : ParticleFilter<RobotProcessModel, MeasurementModel>(n)
{

  lmkModelPtr_ = new LmkProcessModel;
  kfPtr_ = new KalmanFilter(lmkModelPtr_, this->getMeasurementModel());
  
  for(int i = 0; i < n; i++){
    printf("Creating map structure for particle %d\n", i);
    maps_.push_back( new GaussianMixture<TLandmark>() );
    maps_[0]->getGaussianCount();
    unused_measurements_.push_back( std::vector<unsigned int>() );
  }
  
  config.birthGaussianWeight_ = 0.25; 
  config.gaussianMergingThreshold_ = 0.5;
  config.gaussianMergingCovarianceInflationFactor_ = 1.5;
  config.gaussianPruningThreshold_ = 0.2;
  config.importanceWeightingEvalPointCount_ = 8;
  config.importanceWeightingMeasurementLikelihoodMDThreshold_ = 3.0;
  config.newGaussianCreateInnovMDThreshold_ = 0.2;
  config.minInterSampleTimesteps_ = 5;
  
  k_lastResample_ = -10;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::~RBPHDFilter(){

  for(int i = 0; i < maps_.size(); i++){
    delete maps_[i];
  }
  delete kfPtr_;
  delete lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
LmkProcessModel* RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::getLmkProcessModel(){
  return lmkModelPtr_;
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::predict( TInput u,
												 int currentTimestep){

  // Add birth Gaussians using pose before prediction
  addBirthGaussians();

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
}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::update( std::vector<TMeasurement> &Z,
												int currentTimestep){
  k_currentTimestep_ = currentTimestep;

  this->setMeasurements( Z );

  if(!checkMapIntegrity())
    std::cin.get();

  updateMap();

  if(!checkMapIntegrity())
    std::cin.get();

  importanceWeighting();

  // Merge and prune
  for( int i = 0; i < maps_.size(); i++){ // maps_size is same as number of particles

    maps_[i]->merge( config.gaussianMergingThreshold_, 
		     config.gaussianMergingCovarianceInflationFactor_ );

    maps_[i]->prune( config.gaussianPruningThreshold_ );
  }

  if(!checkMapIntegrity())
    std::cin.get();

  bool resampleOccured = false;
  if( k_currentTimestep_ - k_lastResample_ >= config.minInterSampleTimesteps_){
    resampleOccured = this->resample();
  }

  if( resampleOccured ){
    k_lastResample_ = k_currentTimestep_;
  }else{
    this->normalizeWeights();
  }

  // reassign maps as well according to resampling of particles
  if( resampleOccured){

    std::vector< GaussianMixture<TLandmark>* > maps_temp( maps_.size(), NULL );
    std::vector< int > useCount ( maps_.size(), 0);

    // Note which GMs get used and how many times
    for(int i = 0; i < this->nParticles_; i++){
      int j = this->particleSet_[i]->getParentId();
      useCount[j]++;
    }

    // for maps that get used, make a (pointer copy) before doing any overwriting
    // Also rid of the maps that die along with particles
    for(int i = 0; i < this->nParticles_; i++){
      if( useCount[i] > 0){
	maps_temp[i] = maps_[i];
      }else{
	delete maps_[i];
      }
    }
    
    // Copy GMs
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

  if(!checkMapIntegrity())
    std::cin.get();

}

template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::updateMap(){


  const unsigned int nZ = this->measurements_.size();

  for(int i = 0; i < this->nParticles_; i++){

    //---------- 1. setup / book-keeping ----------
    const unsigned int nM = maps_[i]->getGaussianCount();

    if(nM == 0){ // No existing landmark case -> flag all measurements as unused and go to next particle
      unused_measurements_[i].clear();
      for(int z = 0; z < nZ; z++){
	unused_measurements_[i].push_back( z );
      }
      continue;
    }

    double Pd[nM];
    int landmarkCloseToSensingLimit[nM];

    // Mn x nZ table for Gaussian weighting
    double** weightingTable = new double* [ nM ];
    TLandmark*** newLandmarkPointer = new TLandmark** [ nM ];
    for( int n = 0; n < nM; n++ ){
      weightingTable[n] = new double [ nZ ];
      newLandmarkPointer[n] = new TLandmark* [ nZ ];
    }



    //----------  2. Kalman Filter map update ----------

    TPose *pose = new TPose;
    this->particleSet_[i]->getPose(*pose);

    for(unsigned int m = 0; m < nM; m++){

      TLandmark* lm = maps_[i]->getGaussian(m);
      bool isCloseToSensingLimit;
      Pd[m] = this->pMeasurementModel_->probabilityOfDetection( *pose, *lm, 
								isCloseToSensingLimit); 
      landmarkCloseToSensingLimit[m] = ( isCloseToSensingLimit ) ? 1 : 0;
      double w_km = maps_[i]->getWeight(m);
      double Pd_times_w_km = Pd[m] * w_km;

      if(Pd[m] != 0){

	for(int z = 0; z < nZ; z++){

	  TLandmark* lmNew = new TLandmark;
	  newLandmarkPointer[m][z] = NULL;
	  weightingTable[m][z] = 0;
	  double innovationLikelihood = 0;
	  double innovationMahalanobisDist2 = 0;
	  double threshold = config.newGaussianCreateInnovMDThreshold_ * config.newGaussianCreateInnovMDThreshold_;
	
	  // RUN KF, create new landmark for likely updates but do not add to map_[i] yet
	  // because we cannot determine actual weight until the entire weighting table is
	  // filled in
	
	  kfPtr_->correct(*pose, this->measurements_[z], *lm, *lmNew, 
			  &innovationLikelihood, &innovationMahalanobisDist2);

	  if ( innovationMahalanobisDist2 > threshold ){
	    innovationLikelihood = 0;
	    delete lmNew;
	  }else{
	    newLandmarkPointer[m][z] = lmNew;
	  }	
	  weightingTable[m][z] = Pd_times_w_km * innovationLikelihood;
	}

      }else{ // Pd = 0

	for(int z = 0; z < nZ; z++){

	  newLandmarkPointer[m][z] = NULL;
	  weightingTable[m][z] = 0;

	}
      }
    }
    

    // Now calculate the weight of each new Gaussian
    for(int z = 0; z < nZ; z++){
      double clutter = this->pMeasurementModel_->clutterIntensity( this->measurements_[z], nZ );
      double sum = clutter;
      for(unsigned int m = 0; m < nM; m++){
	sum += weightingTable[m][z];
      }
      for(unsigned int m = 0; m < nM; m++){
	weightingTable[m][z] /= sum;
      }
    }


    // ---------- 3. Add new Gaussians to map  ----------
    // New Gaussians will have indices >= nM 
    for(int m = 0; m < nM; m++){
      for(int z = 0; z < nZ; z++){
	if(newLandmarkPointer[m][z] != NULL && weightingTable[m][z] > 0){
	  maps_[i]->addGaussian( newLandmarkPointer[m][z], weightingTable[m][z]);  
	}
      }
    }

    //----------  4. Determine weights for existing Gaussians (missed detection) ----------
    for(int m = 0; m < nM; m++){
      
      double w_km = maps_[i]->getWeight(m);
      double w_k = (1 - Pd[m]) * w_km;

      // For landmarks close to sensing limit
      if (landmarkCloseToSensingLimit[m] == 1){
	double weight_sum_m = 0;
	for(int z = 0; z < nZ; z++){
	  weight_sum_m += weightingTable[m][z];
	}
	double delta_w = Pd[m] * w_km - weight_sum_m;
	if( delta_w > 0 ){
	  w_k += delta_w;
	}
      }

      maps_[i]->setWeight(m, w_k);
    }

    //----------  5. Identify unused measurements for adding birth Gaussians later ----------
    unused_measurements_[i].clear();
    for(int z = 0; z < nZ; z++){
      int useCount = 0;
      for(int m = 0; m < nM; m++){
	if (weightingTable[m][z] != 0){
	  useCount++;
	}
      }
      if (useCount == 0)
	unused_measurements_[i].push_back( z );
    }

    //----------  6. Cleanup - Free memory ----------
    delete pose;

    for( int n = 0; n < nM; n++ ){
      delete[] weightingTable[n];
      delete[] newLandmarkPointer[n];
    }
    delete[] weightingTable;
    delete[] newLandmarkPointer;

  }

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::importanceWeighting(){

  for(int i = 0; i < this->nParticles_; i++){

    //printf("Importance weighting for particle %d\n", i);

    // 1. select evaluation points from highest-weighted Gaussians after update
    const unsigned int nM = maps_[i]->getGaussianCount();
    const int nEvalPoints = config.importanceWeightingEvalPointCount_ > nM ? nM : config.importanceWeightingEvalPointCount_ ;
    if( nEvalPoints == 0 ){
      this->particleSet_[i]->setWeight( std::numeric_limits<double>::denorm_min() );
      continue;
    }
    maps_[i]->sortByWeight(); // sort by weight so that we can pick off the top nEvalPoints Gaussians

    // 2. evaluate sum of Gaussian weights
    double gaussianWeightSumBeforeUpdate = 0;
    double gaussianWeightSumAfterUpdate = 0;
    for(int m = 0; m < nM; m++){
      TLandmark* plm_temp;
      double w, w_prev;
      maps_[i]->getGaussian(m, plm_temp, w, w_prev); // for newly created Gaussians, w_prev = 0
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
    for(int p = 0; p < nEvalPoints; p++){

      TLandmark* lm_evalPt;
      double w_temp;
      maps_[i]->getGaussian(p, lm_evalPt, w_temp);

      double intensity_at_evalPt_beforeUpdate = std::numeric_limits<double>::denorm_min();
      double intensity_at_evalPt_afterUpdate = std::numeric_limits<double>::denorm_min();

      for(int m = 0; m < nM; m++){
	TLandmark* plm;
	double w, w_prev;
	maps_[i]->getGaussian(m, plm, w, w_prev);
	// New Gaussians from update will have w_prev = 0
	// Out Gaussians (missed-detection) will not have been updated, but weights will have changed
	double likelihood = RandomVecMathTools< TLandmark >::evalGaussianLikelihood( *plm, *lm_evalPt);
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
    double measurementLikelihood = rfsMeasurementLikelihood( i, nEvalPoints );
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
      maps_[i]->addGaussian( &landmark_pos, config.birthGaussianWeight_, true);
      
    }
    
  }

}


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
double RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
rfsMeasurementLikelihood( const int particleIdx, const int maxNumOfEvalPoints ){

  // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight; 

  const int i = particleIdx;
  const int nM = maxNumOfEvalPoints;
  const int nZ = this->measurements_.size();

  // Fill in likelihood table
  TPose x;
  this->particleSet_[i]->getPose( x );
  std::vector< double* > likelihoodTab;
  for( int m = 0; m < nM; m++ ){
    
    double* row = new double[nZ];
    
    TLandmark* evalPt; 
    maps_[i]->getGaussian( m, evalPt );
    bool temp;
    double Pd = this->pMeasurementModel_->probabilityOfDetection( x, *evalPt, temp);
    
    double row_sum = 0;
    for( int z = 0; z < nZ; z++ ){

      TMeasurement expected_z;
      this->pMeasurementModel_->measure( x, *evalPt, expected_z);
      TMeasurement actual_z = this->measurements_[z];
      
      double md2;
      double threshold = config.importanceWeightingMeasurementLikelihoodMDThreshold_;
      threshold *= threshold;
      double likelihood = RandomVecMathTools<TMeasurement>::
	evalGaussianLikelihood( actual_z, expected_z, &md2 );
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
      TMeasurement actual_z = this->measurements_[z];
      clutterLikelihood *= this->pMeasurementModel_->clutterIntensity(actual_z, nZ);;
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
      TMeasurement actual_z = this->measurements_[z];
      clutterRow[z] = this->pMeasurementModel_->clutterIntensity(actual_z, nZ);
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
	  TMeasurement actual_z = this->measurements_[z];
	  clutterRow[z] = this->pMeasurementModel_->clutterIntensity(actual_z, nZ);
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
bool RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::checkMapIntegrity(){

  for( int i = 0; i < maps_.size(); i++ ){

    unsigned int nM = maps_[i]->getGaussianCount();
    for(int m = 0; m < nM; m++){

      TLandmark *lm;
      typename TLandmark::Vec lm_x;
      typename TLandmark::Mat lm_S;
      double w;

      maps_[i]->getGaussian(m, lm, w);
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

  if( i >= 0 && i < maps_.size() )
    return ( this->maps_[i]->getGaussianCount() );
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
    maps_[i]->getGaussian(m, plm, w);
    plm->get(u, S);
    return true;
}

#endif


template< class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter >::
setParticlePose(int i, TPose &p){
  
  this->particleSet_[i]->setPose(p);

}
