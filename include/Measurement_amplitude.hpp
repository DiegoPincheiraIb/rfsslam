#ifndef MEASUREMENT_AMPLITUDE_HPP
#define MEASUREMENT_AMPLITUDE_HPP

#include "Measurement.hpp"
namespace rfs{


class Measurement2d_amplitude : public Measurement2d{

public:
  /** Default constructor */
  Measurement2d_amplitude() : Measurement2d(){
    amplitude_=1;
  }
  
  Measurement2d_amplitude(const typename Measurement2d::Vec x, const typename Measurement2d::Mat Sx, const TimeStamp t = TimeStamp() ) :
    Measurement2d(x,Sx,t){
    amplitude_=1;
  }
  Measurement2d_amplitude(const typename Measurement2d::Vec x, const TimeStamp t = TimeStamp()) :
    Measurement2d(x,t){
    amplitude_=1;  
  }
  
  Measurement2d_amplitude(const Measurement2d_amplitude& other):
    Measurement2d(other){
    amplitude_=other.amplitude_;
  }
  Measurement2d_amplitude& operator=(const Measurement2d_amplitude& rhs){
    Measurement2d::operator=(rhs);
    amplitude_=rhs.amplitude_;
    return  *this;
  }
  
  void setAmplitude(double a){
    amplitude_=a;
  }
  double getAmplitude(){
    return amplitude_;
  }
  
  double evalGaussianLikelihood(const Measurement2d_amplitude &x_eval,
				 double* mDist2 = NULL){

  return  Measurement2d::evalGaussianLikelihood(x_eval,mDist2)*exp(-0.5*pow(x_eval.amplitude_-10,2)/2.0)/sqrt(2*PI*2.0);
  }

  
  
protected:
  double amplitude_;

};

}

#endif

