// RB-PHD-Filter Class
// Keith Leung 2013

#include <Eigen/Core>
#include <vector>

#include "Particle.hpp"

typedef Particle* pParticle;

/**
 *  \class RBPHDFilter
 *  \brief Rao-Blackwellized Probability Hypothesis Density Filter class
 *  
 *  This class implements the Rao-Bloackwellized Probability Hypothesis Density
 *  filter. 
 *
 *  \author Keith Leung
 *  \version 0.1
 */

class RBPHDFilter
{
public:

  /** Default consturctor */
  RBPHDFilter();

  /** Destructor */
  ~RBPHDFilter();

  /** 
   *  Initialize the filter
   *  \param nParticles The number of particles
   */
  void init(int nParticles);

private:

  int _nParticles; /**< number of particles in the filter */
  std::vector<pParticle> _particles; /**< vector of particle pointers */

  /** 
   *  Propagate particles using the motion model
   */
  void particles_propagate();

};
