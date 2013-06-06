// Particle filter class
// Keith Leung 2013

#ifndef PARTICLE_FILTER_HPP
#define PARTICLE_FILTER_HPP

#include <cmath> // for random number generation
#include <ctime> // for seeding random number generation
#include <limits> // for random number generation
#include "Particle.hpp"
#include "ProcessModel.hpp"
#include <vector>

/** 
 * \class ParticleFilter
 * \brief A class containing functions for implementing the particle filter
 * \tparam StateType container for the state that the filter updates
 * \tparam SystemInputType container for process model input 
 * \author Keith Leung
 *
 * \todo Function for setting measurement model
 * \todo Function for particle weighting based on measurement likelihood
 * \todo Test this class
 */
template< class StateType, class SystemInputType>
class ParticleFilter
{
public:

  typedef Particle<StateType>* pParticle;
  
  ParticleFilter(){};

  /** 
   * Constructor 
   * \param n number of particles
   * \param initState initial state of particles
   * \param processModelPtr pointer to process model
   */ 
  ParticleFilter(int n, 
		 StateType &initState,
		 ProcessModel<StateType, SystemInputType>* processModelPtr);

  /** Default destructor */
  ~ParticleFilter();

  /** 
   * Set the process model to use for particle propagation
   * \param model pointer to process model
   */
  void setProcessModel( ProcessModel<StateType, SystemInputType>* modelPtr );

  /** 
   * Propagate particles using the process model
   * \param input to the process model
   * \double dT time-step of input (not used by all process models)
   */
  void propagate( SystemInputType &input, double const dt = 0);


  /** 
   * Normalize particle weights so that they sum to 1
   */
  void normalizeWeights();

  /**
   * Particling resampling using a low variance sampling
   * \param n number of particles in the resampled set.
   *        The default value of 0 will keep the same number of particles 
   */
  void resample( unsigned int n = 0);
  

protected:

  int nParticles_; /**< number of particles */
  std::vector< pParticle > particleSet_; /**< container for particle pointers */

  ProcessModel<StateType, SystemInputType>* pProcessModel_;

private:
};

////////// Implementation //////////

template< class StateType, class SystemInputType >
ParticleFilter<StateType, SystemInputType>::ParticleFilter(int n, 
					  StateType &initState,
					  ProcessModel<StateType, SystemInputType>* processModelPtr){
  
  // initiate particles
  nParticles_ = n;
  particleSet_.resize(nParticles_);
  double newParticleWeight = 1;
  for( int i = 0 ; i < nParticles_ ; i++ ){
    particleSet_[i] = new Particle<StateType>(n, initState, newParticleWeight);
  }
  
  // set prcoess model
  setProcessModel( processModelPtr );

  // set random seed for particle resampling
  srand((unsigned int)time(NULL));
  
}


template< class StateType, class SystemInputType >
ParticleFilter<StateType, SystemInputType>::~ParticleFilter(){
   for( int i = 0 ; i < nParticles_ ; i++ ){
     delete particleSet_[i];
   }
}

template< class StateType, class SystemInputType >
void ParticleFilter<StateType, SystemInputType>::setProcessModel( ProcessModel<StateType, SystemInputType>* modelPtr ){
  pProcessModel_ = modelPtr;
}


template< class StateType, class SystemInputType >
void ParticleFilter<StateType, SystemInputType>::propagate( SystemInputType &input, 
					   double const dt){
   
  for( int i = 0 ; i < nParticles_ ; i++ ){
    StateType x_km, x_k;
    particleSet_[i]->getPose( x_km );
      pProcessModel_->step( x_k, x_km, input, dt);
    particleSet_[i]->setPose( x_k );
  } 
}

template< class StateType, class SystemInputType >
void ParticleFilter<StateType, SystemInputType>::normalizeWeights(){
  
  double sum = 0;
  for( int i = 0; i < nParticles_; i++ ){
    sum += particleSet_[i]->getWeight();
  }
  for( int i = 0; i < nParticles_; i++ ){
    particleSet_[i]->setWeight( particleSet_[i]->getWeight() / sum );
  }

}

template< class StateType, class SystemInputType >
void ParticleFilter<StateType, SystemInputType>::resample( unsigned int n ){

  if( n == 0 )
    n = nParticles_; // number of particles to sample

  // normalize particle weights so they sum to 1
  normalizeWeights();
  
  // Sampler settings
  double randomNum_0_to_1 = ((double) rand() / (RAND_MAX + 1));
  unsigned int idx = 0;
  const double sample_interval = 1.0 / double(n); 
  const double sampler_offset = sample_interval * randomNum_0_to_1;
  double sample_point = sampler_offset;
  double cumulative_weight = particleSet_[idx]->getWeight();

  // book-keeping
  std::vector<char> flag_particle_sampled (nParticles_, 0); 
  std::vector<unsigned int> sampled_idx (n, 0);
  
  // Sample
  for( int i = 0; i < n; i++ ){
    while( sample_point > cumulative_weight ){
      // particle[idx] not sampled
      idx++;
      cumulative_weight += particleSet_[idx]->getWeight();
    }
    // particle[idx] sampled
    sampled_idx[i] = idx;
    flag_particle_sampled[idx] = 1;
    sample_point += sample_interval;
  }

  // Do the actual data copying
  unsigned int idx_prev = 0;
  unsigned int next_unsampled_idx = 0;
  for( int i = 0; i < n; i++ ){
    
    bool firstTime = true;
    idx = sampled_idx[i]; // particle[idx] was sampled
    if( i > 0){
      if( idx == idx_prev ){
	firstTime = false;
      }
    }
    idx_prev = idx;

    while( flag_particle_sampled[next_unsampled_idx] == 1)
      next_unsampled_idx++;

    // cases:
    // 1. idx < n AND idx appears for first time -> do nothing
    // 2. idx < n AND it is not the first time that idx appears -> make copy
    // 3. idx > n AND idx appears for first time -> make copy
    // 4. idx > n AND it is not first time that idx appears -> make copy

    if( idx < n && firstTime){ // case 1
      particleSet_[idx]->setParentId( particleSet_[idx]->getId() );
    }else{ // case 2, 3, 4
      particleSet_[next_unsampled_idx]->setParentId( particleSet_[idx]->getId() );
      particleSet_[i]->copyStateTo( particleSet_[next_unsampled_idx] );
    }
      
  }

  // Delete all particles with idx >= n
  for( int i = n; i < nParticles_; i++ ){
    delete particleSet_[i];
  }
  nParticles_ = n;

  // Reset weight of all particles
  for( int i = 0; i < n; i++ ){
    particleSet_[i]->setWeight(1);
  }
  
}


#endif
