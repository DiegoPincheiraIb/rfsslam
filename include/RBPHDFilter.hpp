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

  typedef typename ProcessModel::TPose TPose;
  typedef typename ProcessModel::TInput TInput;
  typedef typename MeasurementModel::tLandmark TLandmark;
  typedef typename MeasurementModel::tMeasurement TMeasure;
  typedef typename GaussianMixture<TLandmark>::Gaussian TGaussian;

  /** 
   * \brief Configurations for running the RB-PHD-Filter 
   */
  struct Config{
    
    /**  new birth Gaussians are set with this weight */
    double birthGaussianWeight_;   

    /**  number of map states to use for evaluating particle weight
	 0 => empty-set strategy,
	 1 => single-feature strategy,
	 >1 => multi-feature strategy
    */
    int nImportanceWeightingEvalPoints_;

    /**  new Gaussians are only created during map update if the innovation likelihood 
	 is greater than this threshold */
    double newGaussianCreateLikelihoodThreshold_; 

  } config;

  /** 
   * Constructor 
   * \param n number of particles
   * \param initState initial state of particles
   * \param processModelPtr pointer to the process model
   * \param measurementModelPtr pointer to the measurement model
   */
  RBPHDFilter(int n, 
	      TPose &initState,
	      ProcessModel* processModelPtr,
	      MeasurementModel* measurementModelPtr)
    : ParticleFilter<ProcessModel, MeasurementModel>(n, initState, processModelPtr, measurementModelPtr){
    
    for(int i = 0; i < n; i++){
      maps_.push_back( new GaussianMixture<TLandmark>() );
    }

    config.birthGaussianWeight_ = 0.25; 
    config.nImportanceWeightingEvalPoints_ = 8;
    config.newGaussianCreateLikelihoodThreshold_ = 0.2; 
  }

  /** Destructor */
  ~RBPHDFilter();

  /**
   * Update the map, calculate importance weighting, sample if necessary, and
   * create new birth Gaussians.
   * \param Z set of measurements to use for the update, placed in a stl vector, which
   * gets cleared after the function call. 
   */
  void update( std::vector<TMeasure> &Z );


private:

  std::vector< GaussianMixture<TLandmark>* > maps_; /**< Particle dependent maps */

  /** indices of unused measurement for each particle for creating birth Gaussians */
  std::vector< std::vector<unsigned int> > unused_measurements_; 
  
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

};

////////// Implementation //////////

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::~RBPHDFilter(){

  for(int i = 0; i < maps_.size(); i++){
    delete maps_[i];
  }
}

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::update( std::vector<TMeasure> &Z ){

  this->setMeasurements( Z );
  updateMap();
  importanceWeighting();
  // \todo merge and prune
  this->resample();
  // \todo need to reassign maps as well according to resampling
  addBirthGaussians();

}

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::updateMap(){

  for(int i = 0; i < this->nParticles_; i++){

    //---------- 1. setup / book-keeping ----------
    const unsigned int nM = maps_[i]->getGaussianCount();
    const unsigned int nZ = this->measurements_.size();
    double Pd[nM];
    int landmarkCloseToSensingLimit[nM];

    // nM x nZ table for Gaussian weighting
    double** weightingTable = new double* [ nM ];
    TGaussian** newLandmarkPointer = new TGaussian* [ nM ];
    for( int n = 0; n < nM; n++ ){
      weightingTable[n] = new double [ nZ ];
      newLandmarkPointer[n] = new TGaussian [ nZ ];
    }

    for(int m = 0; m < nM; m++){
      double w = maps_[i]->getWeight(m);
    }

    //----------  2. Kalman Filter map update ----------

    const TPose pose = this->particleSet_[i]->getPose();

    for(int m = 0; m < nM; m++){

      TLandmark* lm = maps_[i]->landmark;
      bool isCloseToSensingLimit;
      Pd[m] = this->pMeasurementModel_->probabilityOfDetection( pose, *lm, isCloseToSensingLimit); 
      landmarkCloseToSensingLimit[m] = ( isCloseToSensingLimit ) ? 1 : 0;
      double w_km = maps_[i]->getWeight(m);
      double Pd_times_w_km = Pd[m] * w_km;

      for(int z = 0; z < nZ; z++){

	newLandmarkPointer[m][z] = NULL;
	weightingTable[m][z] = 0;

	// RUN KF, create new landmark for likely updates but do not add to map_[i] yet
	// because we cannot determine actual weight until the entire weighting table is
	// filled in

	double innovationLikelihood = 0;
	// \todo innovationLikelihood = ???
	if ( innovationLikelihood < config.newGaussianCreateLikelihoodThreshold ){
	  innovationLikelihood = 0;
	}
	TGaussian* lmPtr = new TGaussian;
	newLandmarkPointer[m][z] = lmPtr;
	weightingTable[m][z] = Pd_times_w_km * innovationLikelihood;

      }

    }

    for(int z = 0; z < nZ; z++){
      double clutter = this->pMeasurementModel_->clutterIntensity( this->measurements_[z], nZ );
      double sum = clutter;
      for(int m = 0; m < nM; m++){
	sum += weightingTable[m][z];
      }
      for(int m = 0; m < nM; m++){
	weightingTable[m][z] /= sum;
      }
    }


    // ---------- 3. Add new Gaussians to map  ----------

    for(int m = 0; m < nM; m++){
      for(int z = 0; z < nZ; z++){
	maps_[i]->addGaussian( newLandmarkPointer[m][z], weightingTable[m][z]);  
      }
    }

    //----------  4. Determine weights for existing Gaussians (missed detection) ----------
    for(int m = 0; m < nM; m++){
      
      double w_km = maps_[i]->getWeight(m);
      double w_k = (1 - Pd[m]) * w_km;

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

    //----------  5. Identity unused measurements for adding birth Gaussians later ----------
    unused_measurements_[i].clear();
    for(int z = 0; z < nZ; z++){
      bool unused = true;
      for(int m = 0; m < nM; m++){
	if (weightingTable[m][z] != 0){
	  unused = false;
	  break;
	}
      }
      if (unused)
	unused_measurements_[i].push_back( z );
    }

    //----------  6. Cleanup - Free memory ----------
    for( int n = 0; n < nM; n++ ){
      delete[] weightingTable[n];
      delete[] newLandmarkPointer[n];
    }
    delete[] weightingTable;
    delete[] newLandmarkPointer;

  }

}


template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::importanceWeighting(){

  for(int i = 0; i < this->nParticles_; i++){

    // 1. select evaluation points from highest-weighted Gaussians after update
    const int nEvalPoints = config.nImportanceWeightingEvalPoints_;
    maps_[i]->sortByWeight();
    const unsigned int nM = maps_[i]->getGaussianCount();

    // 2. evaluate sum of Gaussian weights
    double gaussianWeightSumBeforeUpdate = 0;
    double gaussianWeightSumAfterUpdate = 0;
    for(int m = 0; m < nM; m++){
      TLandmark* plm_temp;
      double w, w_prev;
      maps_[i]->getGaussian(m, plm_temp, w, w_prev);
      gaussianWeightSumBeforeUpdate += w_prev;
      gaussianWeightSumAfterUpdate += w;
    }

    // 3. evaluate intensity function at eval points and take their product
    double intensityProd_beforeUpdate = 1;
    double intensityProd_afterUpdate = 1;
    for(int p = 0; p < nEvalPoints; p++){

      TLandmark* lm_evalPt;
      double w_temp;
      maps_[i]->getGaussian(p, lm_evalPt, w_temp);

      double intensity_at_evalPt_beforeUpdate = 1;
      double intensity_at_evalPt_afterUpdate = 1;
	
      for(int m = 0; m < nM; m++){
	TLandmark* plm;
	double w, w_prev;
	maps_[i]->getGaussian(m, plm, w, w_prev);
	intensity_at_evalPt_beforeUpdate += w_prev * plm->evalGaussianLikelihood( *lm_evalPt );
	intensity_at_evalPt_afterUpdate += w * plm->evalGaussianLikelihood( *lm_evalPt );
      }

      intensityProd_beforeUpdate *= intensity_at_evalPt_beforeUpdate;
      intensityProd_afterUpdate *= intensity_at_evalPt_afterUpdate;
    }

    // 4. calculate measurement likelihood at eval points
    
    double measurementLikelihood = 1;
    std::vector<TLandmark> a;
    for(int p = 0; p < nEvalPoints; p++){

      TLandmark* lm_evalPt;
      double w_temp;
      maps_[i]->getGaussian(p, lm_evalPt, w_temp);

      rfsMeasurementLikelihood( config.nImportanceWeightingEvalPoints_ );

    }

    // 5. calculate overall weight
    double overall_weight = measurementLikelihood * intensityProd_beforeUpdate / intensityProd_afterUpdate *
      exp( gaussianWeightSumAfterUpdate - gaussianWeightSumAfterUpdate); 
    this->particleSet_[i]->setWeight( overall_weight );

  }

}

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::addBirthGaussians(){

  for(int i = 0; i < this->nParticles_; i++){
    
    while( unused_measurements_[i].size() > 0){
     
      // get measurement
      int unused_idx = unused_measurements_[i].back();
      TMeasure unused_z = this->measurements_[unused_idx];
      unused_measurements_[i].pop_back();

      // use inverse measurement model to get landmark
      TPose robot_pose;
      TLandmark landmark_pos;
      this->particleSet_[i]->getPose(robot_pose);
      this->pMeasurementModel_->inversePredict( robot_pose, landmark_pos, unused_z );

      // add birth landmark to Gaussian mixture (last param = true to allocate mem)
      maps_[i]->addGaussian( &landmark_pos, config.birthGaussianWeight_, true);
      
    }
    
  }

}



template< class ProcessModel, class MeasurementModel, class KalmanFilter >
double RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::
rfsMeasurementLikelihood( const int particleIdx, const int maxNumOfEvalPoints ){

  // eval points are first nEvalPoints elements of maps_[i]; 

  const int i = particleIdx;
  const int nM = maxNumOfEvalPoints;
  const int nZ = this->measurements_->size();
  int nClutter = 0;

  // Allocate memory for likelihood table
  std::vector< double* > likelihoodTab(nM);
  for( int m = 0; m < nM; m++ ){
    likelihoodTab[m] = new double [nZ];
  }

  // Fill in likelihood table
  TPose x;
  this->particleSet_[i]->getPose( x );

  for( int m = 0; m < nM; m++ ){

    TLandmark evalPt; 
    maps_[i]->getGaussian( m, evalPt );
    bool temp;
    double Pd = this->pMeasurementModel_->probabilityOfDetection( x, evalPt, temp);

     for( int z = 0; z < nZ; z++ ){

       TMeasure expected_z;
       this->pMeasurementModel_->predict( x, evalPt, expected_z);
       TMeasure actual_z = this->measurements_[z];

       double likelihood = actual_z.evalGaussianLikelihood( expected_z );
       
       if( 1 ){ // \todo check against threshold
	 likelihoodTab[m][z] = likelihood * Pd;
       }else{
	 likelihoodTab[m][z] = 0;
       }

     }
  }

  // Check measurements with 0 likelihood to all eval points
  // if so, that measurement is considered clutter
  for( int z = 0; z < nZ; z++ ){
    double z_sum = 0;
    for( int m = 0; m < nM; m++ ){
      z_sum += likelihoodTab[m][z];
    }
    if( z_sum = 0 ){
      nClutter++;
      TMeasure actual_z = this->measurements_[z];
      double c = this->measurementModel_->clutterIntensity(actual_z, nZ);
      for( int m = 0; m < nM; m++ ){
	likelihoodTab[m][z] = c;
      }
    }
  }

  // If the number of measurements is greater than the number of
  // eval points, then some measurements must be clutter.
  // We will add extra rows for these clutter measurements in 
  // likelihoodTab
  int nC = nZ - nM;
  if (nC < 0){
    nC = 0;
  }else{
    nClutter += nC;
  }
  for( int c = 0; c < nC; c++ ){
    likelihoodTab.push_back(new double [nZ]);
    int m = likelihoodTab.size() - 1;
    for( int z = 0; z < nZ; z++ ){
      TMeasure actual_z = this->measurements_[z];
      likelihoodTab[m][z] = this->measurementModel_->clutterIntensity(actual_z, nZ);
    }
  }

  // Go through all permutations of eval point - measurement
  // to calculate the likelihood
  double likelihood = 0;
  while (likelihood == 0){

    // likelihood = ?


    if( likelihood == 0 ){
      // Add another row to of clutter to likelihoodTab
      likelihoodTab.push_back(new double [nZ]);
      int m = likelihoodTab.size() - 1;
      for( int z = 0; z < nZ; z++ ){
	TMeasure actual_z = this->measurements_[z];
	likelihoodTab[m][z] = this->measurementModel_->clutterIntensity(actual_z, nZ);
      }
    }

  }

  // Deallocate likelihood table
  for( int m = 0; m < nM + nC; m++ ){
    delete[] likelihoodTab[m];
  }

  if (nClutter > 0){
    likelihood /= this->measurementModel_->clutterIntensityIntegral( nZ );
  }

  return likelihood;
}

#endif
