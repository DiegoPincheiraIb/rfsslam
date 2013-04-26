// RB-PHD-Filter Class
// Keith Leung 2013

#include <eigen>
#include <vector>

/**
 *  @class RBPHDFilter
 *  @brief Rao-Blackwellized Probability Hypothesis Density Filter class
 *  
 *  This class implements the Rao-Bloackwellized Probability Hypothesis Density
 *  filter. 
 *
 *  @author Keith Leung
 *  @version 0.1
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
   *  @param[in] nParticles The number of particles
   */
  void init(int nParticles);

private:

  int _nParticles;

};
