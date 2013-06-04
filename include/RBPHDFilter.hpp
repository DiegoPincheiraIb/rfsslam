// RB-PHD-Filter Class
// Keith Leung 2013

#include <Eigen/Core>
#include <vector>
#include "ParticleFilter.hpp"

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

// template< class PoseType, class SystemInputType >
class RBPHDFilter //: public ParticleFilter<PoseType, SystemInputType>
{
public:

  // typedef Particle<PoseType>* pParticle;

  /** Default consturctor */
  RBPHDFilter();

  /** Destructor */
  ~RBPHDFilter();

private:

};
