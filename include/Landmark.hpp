// Landmark classes for defining map feature state
// Keith Leung 2013

#ifndef LANDMARK_HPP
#define LANDMARK_HPP

#include "RandomVec.hpp"

/** 
 * \class Landmark
 * \brief An abstract class for defining landmark state
 * \author Keith Leung
 */

template<class VecType, class MatType, class DescriptorType = int>
class Landmark : public RandomVec<VecType, MatType>
{
public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /** Default constructor */
  Landmark(){}

  /** 
   * Constructor
   */
  Landmark(VecType &x, MatType &Sx){
    this->set(x, Sx);
  }

  /** 
   * Constructor
   */
  Landmark(VecType &x, MatType &Sx, DescriptorType &d){
    this->set(x, Sx);
    desc_ = d;
  }

  /** Default destructor */
  ~Landmark(){};

  /** Set descriptor for landmark 
   *  \param[in] d descriptor
   */
  void setDescriptor(DescriptorType &d){
    desc_ = d;
  }

  /** Get descriptor for landmark 
   *  \param[in] d descriptor
   */
  void getDescriptor(DescriptorType &d){
    d = desc_;
  }

private:
  
  DescriptorType desc_;

};

typedef Landmark< Eigen::Matrix<double, 1, 1>, Eigen::Matrix<double, 1, 1> >
Landmark1d;

typedef Landmark<Eigen::Vector2d, Eigen::Matrix2d> Landmark2d;

#endif
