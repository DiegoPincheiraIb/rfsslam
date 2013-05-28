// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include <Eigen/Core>

/** 
 * \class Landmark
 * \brief An abstract class for defining landmark state
 * \author Keith Leung
 */
class Landmark
{
public:

  /** Default constructor */
  Landmark(){};

  /** Default destructor */
  ~Landmark(){};

  /** Abstract function for getting landmark state */
  virtual void getState(){};

  /** Abstract function for setting landmark state */  
  virtual void setState(){};

};

/** 
 * \class Landmark2d
 * \brief 2d Landmark
 * \author Keith Leung
 */
class Landmark2d : public Landmark
{
public:

  typedef Eigen::Vector2d vec;

  /** Default constructor */
  Landmark2d();

  /** Constructor */
  Landmark2d(vec x);

  /** Default destructor */
  ~Landmark2d();

  /** 
   * Function for setting landmark state
   */
  void setState(vec x);

  /** 
   * Function for getting landmark state
   * \param x state passed in by reference to be overwritten
   */
  void getState(vec& x);

private:

  vec x_; /**< state */

};


/** 
 * \class Landmark1d
 * \brief 1d Landmark
 * \author Keith Leung
 */
class Landmark1d : public Landmark
{
public:

  /** Default constructor */
  Landmark1d();

  /** Constructor */
  Landmark1d(double x);

  /** Default destructor */
  ~Landmark1d();

  /** 
   * Function for setting landmark state
   */
  void setState(double x);

  /** 
   * Function for getting landmark state
   * \param x state passed in by reference to be overwritten
   */
  void getState(double& x);

private:

  double x_; /**< state */

};

#endif
