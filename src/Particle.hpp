// Particle Class used in the RB-PHD-Filter
// Keith Leung 2013

#ifndef PARTICLE_HPP
#define PARTICLE_HPP

/** 
 *  \class Particle
 *  \brief A class for a particle for the particle filter
 *  \tparam PoseType container for the state
 *  \author Keith Leung
 */
template< class PoseType >
class Particle
{

public:

  /** Default constructor */
  Particle();

  /* Constructor 
   * \param id particle id
   * \param x_k_i particle pose
   * \param w particle weight
   */
  Particle( unsigned int id, PoseType &x_k_i, double w = 0 );

  /** Destructor */
  ~Particle();

  /**
   *  Set particle pose
   *  \param x_k_i pose
   */
  void setPose( PoseType &x_k_i );

  /**
   *  Get particle pose
   *  \param x_k_i pose [overwritten]
   */
  void getPose( PoseType &x_k_i );

  /**
   *  Get particle pose
   *  \return pointer to pose
   */
  const PoseType* getPose();

  /**
   *  Set particle importance weight
   *  \param w weight
   */
  void setWeight( double w );

  /** 
   *  Get particle importance weight
   *  \return particle weight
   */
  double getWeight();

protected:

  unsigned int id_; /**< particle id number */
  PoseType x_k_i_; /**< robot pose at timestep k **/
  double w_; /**< Particle weight */

};

// Implementation

template< class PoseType >
Particle<PoseType>::Particle(){
  id_ = 0;
  w_ = 0;
}

template< class PoseType >
Particle<PoseType>::Particle( unsigned int id, PoseType &x_k_i, double w ){
  id_ = id;
  x_k_i_ = x_k_i;
  w_ = w;
}

template< class PoseType >
Particle<PoseType>::~Particle(){}

template< class PoseType >
void Particle<PoseType>::setPose( PoseType &x_k_i ){
  x_k_i = x_k_i_;
}

template< class PoseType >
void Particle<PoseType>::getPose( PoseType &x_k_i ){
  x_k_i_ = x_k_i;
}

template< class PoseType >
const PoseType* Particle<PoseType>::getPose(){
  return &x_k_i_;
}

template< class PoseType >
void Particle<PoseType>::setWeight( double w ){
  w_ = w;
}

template< class PoseType >
double Particle<PoseType>::getWeight(){
  return w_;
}


#endif
