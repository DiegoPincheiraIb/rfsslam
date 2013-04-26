#ifndef GAUSSIAN_MIXTURE_HPP
#define GAUSSIAN_MIXTURE_HPP

/** 
 * \class GaussianMixture
 * \brief A class for mixture of Gaussians
 *
 * This class represents a mixture of weighted Gaussians using a binary tree.
 * It is designed to work in a Rao-Blackwellized particle filter
 *
 * \author Keith Leung
 */

class GaussianMixture
{
public:
  
  /** Default constructor */
  GaussianMixture();
  
  /** Desctructor */
  ~GaussianMixture();

private:

  int _n; /**< number of Gaussians */
};

#endif
