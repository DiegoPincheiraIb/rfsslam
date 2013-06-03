// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include <Eigen/Core>
#include "State.hpp"

/** 
 * \class Landmark
 * \brief An abstract class for defining landmark state
 * \author Keith Leung
 */
template<class StateType, class UncertaintyType>
class Landmark : public StateWithUncertainty<StateType, UncertaintyType>
{
public:

  /** Default constructor */
  Landmark(){};

  /** Default destructor */
  ~Landmark(){};

};




class Landmark1D : public Landmark<double, double>
{
public: 
  Landmark1D();
  ~Landmark1D();
};

#endif
