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

  /** \todo function for setting measurement model */

  /** 
   * Propagate particles using the process model
   * \param input to the process model
   * \double dT time-step of input (not used by all process models)
   */
  void propagate( SystemInputType &input, double const dt = 0);

  /** \todo function for particle weighting - requires measurement model to be completed first */

  /** \todo function for resampling of particles */

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
    n = nParticles_;

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

  // Replaced unsampled particles with sampled particles
  
  
}


#endif
