// Particle Class used in the RB-PHD-Filter
// Keith Leung 2013

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "GaussianMixture.hpp"

/** 
 *  \class particle
 *  \brief A class for a particle in the rao-blackwellized particle filter
 *
 *  \author Keith Leung
 */
class Particle
{

public:

  /** Default constructor */
  Particle();

  /** Destructor */
  ~Particle();

private:

  /** \todo make the robot state a template */
  double _x; /**< robot state **/
  GaussianMixture _map; /**< feature -- a mixture of Gaussians */

};


#endif
