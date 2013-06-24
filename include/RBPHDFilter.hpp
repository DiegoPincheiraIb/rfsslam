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
      gaussianWeightSumBeforeUpdate_.push_back(0);
      gaussianWeightSumAfterUpdate_.push_back(0);
      gaussianWeightsBeforeUpdate_.push_back( std::vector<double>() );
    }

    default_birthGaussianWeight_ = 0.25; // \todo make this part of config
    newGaussianCreateLikelihoodThreshold_ = 0.2; // \todo make this part of config
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

  /** Number of landmarks estimated for each particle before map update */
  std::vector< double > gaussianWeightSumBeforeUpdate_;

  /** Number of landmarks estimated for each particle after map update */
  std::vector< double > gaussianWeightSumAfterUpdate_;

  /** Weight of Gaussians before map update. This is required for the important weighting step */
  std::vector< std::vector<double> > gaussianWeightsBeforeUpdate_;
  
  /** The Gaussian mixture weight assigned to new birth Gaussians */
  double default_birthGaussianWeight_; 
  
  /** The likelihood threshold for creating new Gaussians during map update */
  double newGaussianCreateLikelihoodThreshold_;

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
  this->resample();
  addBirthGaussians();

}

template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::updateMap(){

  for(int i = 0; i < this->nParticles_; i++){

    //---------- 1. setup / book-keeping ----------
    const unsigned int nM = maps_[i]->getGaussianCount();
    const unsigned int nZ = this->measurements_.size();
    double Pd[nM];

    // nM x nZ table of measurement likelihood, given a map prior
    double** likelihoodTable = new double* [ nM ];
    TLandmark** newLandmarkPointer = new TGaussian* [ nM ];
    for( int n = 0; n < nM; n++ ){
      likelihoodTable[n] = new double [ nZ ];
      newLandmarkPointer[n] = new TGaussian [ nZ ];
    }

    gaussianWeightsBeforeUpdate_[i].clear();
    gaussianWeightsBeforeUpdate_[i].resize(nM);
    gaussianWeightSumBeforeUpdate_[i] = 0;
    for(int m = 0; m < nM; m++){
      double w = maps_[i]->getWeight(m);
      gaussianWeightsBeforeUpdate_[i][m] = w;
      gaussianWeightSumBeforeUpdate_[i] += w;
    }

    //----------  2. Kalman Filter map update ----------
    for(int m = 0; m < nM; m++){

      Pd[m] = 0; // \todo in measurement model
      double Pfa = 0; // \todo in measurement model
      double c = 0; // \todo clutter model

      for(int z = 0; z < nZ; z++){

	newLandmarkPointer[m][z] = NULL;
	likelihoodTable[m][z] = 0;

	// Run Kalman Filter
	// Create new landmark for likely updates but do not add to map_[i] yet
	TGaussian* lmPtr = new TLandmark;
	// lmPtr->landmark = todo;
	newLandmarkPointer[m][z] = lmPtr;
	likelihoodTable[m][z] = 0; // \todo update likelihood table

      }

    }

    //----------  3. Identity unused measurements for adding birth Gaussians later ----------
    unused_measurements_[i].clear();
    for(int z = 0; z < nZ; z++){
      bool unused = true;
      for(int m = 0; m < nM; m++){
	if (likelihoodTable[m][z] != 0){
	  unused = false;
	  break;
	}
      }
      if (unused)
	unused_measurements_[i].push_back( z );
    }


    // ---------- 4. Add new Gaussians to map and determine weights based on the likelihood table ----------
    gaussianWeightSumAfterUpdate_[i] = 0;
    for(int m = 0; m < nM; m++){
      for(int z = 0; z < nZ; z++){

	double weight_m_z = 0; // \todo calculate weight for new Gaussian
	gaussianWeightSumAfterUpdate_[i] += weight_m_z;
	maps_[i]->addGaussian( newLandmarkPointer[m][z], weight_m_z);  

      }
    }

    //----------  5. Determine weights for existing Gaussians (missed detection) ----------
    for(int m = 0; m < nM; m++){
      double w_km = maps_[i]->getWeight(m);
      double w_k = (1 - Pd[m]) * w_km;
      maps_[i]->setWeight(m, w_k);
      gaussianWeightSumAfterUpdate_[i] += w_k;
    }

    //----------  6. Cleanup - Free memory ----------
    for( int n = 0; n < nM; n++ ){
      delete[] likelihoodTable[n];
      delete[] newLandmarkPointer[n];
    }
    delete[] likelihoodTable;
    delete[] newLandmarkPointer;

  }

  // \todo implement measurement ambiguity zone
  
  // \todo multithread this
}


template< class ProcessModel, class MeasurementModel, class KalmanFilter >
void RBPHDFilter< ProcessModel, MeasurementModel, KalmanFilter >::importanceWeighting(){

  // \todo implement all three weighting strategies

  for(int i = 0; i < this->nParticles_; i++){

    // \todo select evaluation points from highest peaks in the map after update
    // upper limit of map set size should be the number of measurements
    // select evaluation points from highest-weighted Gaussians

    const int nEvalPoints = 0; // \todo

    const unsigned int nM = maps_[i]->getGaussianCount();

    // \todo sort and get top nEvalPoints points

    // \todo evaluate intensity function at eval points
    for(int m = 0; m < nM; m++){

    }

    // \todo calculate measurement likelihood at eval points

    // \todo calculate overall weight

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
      maps_[i]->addGaussian( &landmark_pos, default_birthGaussianWeight_, true);
      
    }
    
  }

}

#endif
