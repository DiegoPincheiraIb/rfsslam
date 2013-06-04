// Particle filter class
// Keith Leung 2013

#ifndef PARTICLE_FILTER_HPP
#define PARTICLE_FILTER_HPP

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


#endif
