#ifndef MEASUREMENT_AMPLITUDE_HPP
#define MEASUREMENT_AMPLITUDE_HPP

#include "Measurement.hpp"
namespace rfs{

/**
 * This Class implements a Measurement that includes amplitude information and calculates
 * its likelihood using a one dimensional Normal distribution.
 */
class Measurement2d_amplitude : public Measurement2d{

public:
  /** Default constructor */
  Measurement2d_amplitude() : Measurement2d(){
    amplitude_=1;
    mean_=10;
    cov_=2;
  }
  /**
   * 
   */
  Measurement2d_amplitude(const Measurement2d::Vec x, const Measurement2d::Mat Sx, const TimeStamp t = TimeStamp() ) :
    Measurement2d(x,Sx,t){
    amplitude_=1;
    mean_=10;
    cov_=2;
  }
  /** 
   * Constructor 
   * \param[in] t time
   */ 
  Measurement2d_amplitude(const Measurement2d::Vec x, const TimeStamp t = TimeStamp()) :
    Measurement2d(x,t){
    amplitude_=1;
    mean_=10;
    cov_=2; 
  }
  /**
   * Copy constructor
   * \param[in] other the Measurement_amplitude being copied from
   */
  Measurement2d_amplitude(const Measurement2d_amplitude& other):
    Measurement2d(other){
    amplitude_=other.amplitude_;
  }
  /**
   * Assignment operator
   * \param[in] rhs the right-hand-side from which data is copied
   */  
  
  Measurement2d_amplitude& operator=(const Measurement2d_amplitude& rhs){
    Measurement2d::operator=(rhs);
    amplitude_=rhs.amplitude_;
    return  *this;
  }
  /**
   * Setter for the amplitude value.
   * \param[in] a Amplitude value to calculate the measurement likelihood.
   */
  void setAmplitude(double a){
    amplitude_ = a;
  }
  /**
   * Getter for the amplitude value.
   * \param[in] a Amplitude value to calculate the measurement likelihood.
   * \return the amplitude value
   */
  double getAmplitude(){
    return amplitude_;
  }
  
  /**
   * Modified likelihood including the likelihood of the measurement amplitude value.
   * \param[in] x_eval the measurement whose mean and amplitude will be used to calculate the likelihood
   * \param[out] mDist2 Squared Mahalanobis distance (does not include amplitude information)
   */
  double evalGaussianLikelihood(const Measurement2d_amplitude &x_eval,
				 double* mDist2 = NULL){

  return  Measurement2d::evalGaussianLikelihood(x_eval,mDist2)*exp(-0.5*pow(x_eval.amplitude_-mean_,2)/cov_)/sqrt(2*PI*cov_);
  }

  
  
protected:
  double mean_;
  double cov_;
  double amplitude_;

  };

}

#endif

