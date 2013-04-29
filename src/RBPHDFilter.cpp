// RB-PHD-Filter Class Implementation
// Keith Leung 2013

#include "RBPHDFilter.hpp"

RBPHDFilter::RBPHDFilter(){

}

RBPHDFilter::~RBPHDFilter(){
  
  // deallocate particles
  for(int i = 0; i < _nParticles; i++){
    delete _particles[i];
  }
}

void RBPHDFilter::init( int nParticles ){
  
  _nParticles = nParticles;
  _particles = std::vector<pParticle>(_nParticles);
  for(int i = 0; i < _nParticles; i++){
    _particles[i] = new Particle();
  }

}

void RBPHDFilter::particles_propagate(){

  for(int i = 0; i < _nParticles; i++){
    /** \todo implement propagation of each particle state usig motion model */
  }

  return;
}
