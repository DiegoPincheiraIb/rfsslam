/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2015, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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

    /** \todo Documentation */
    Measurement2d_amplitude(const Measurement2d::Vec &x,
			    const Measurement2d::Mat &Sx,
			    const TimeStamp t = TimeStamp() ) :
      Measurement2d(x,Sx,t){
      amplitude_=1;
    }

    /** \todo Documentation */
    Measurement2d_amplitude(const Measurement2d::Vec &x,
			    const TimeStamp t = TimeStamp()) :
      Measurement2d(x,t){
      amplitude_=1;  
    }

    /* copy constructors and operator= not necessary
    Measurement2d_amplitude(const Measurement2d_amplitude& other):
      Measurement2d(other){
      amplitude_=other.amplitude_;
    }
    
    Measurement2d_amplitude& operator=(const Measurement2d_amplitude& rhs){
      Measurement2d::operator=(rhs);
      amplitude_=rhs.amplitude_;
      return  *this;
    }
    */

    /** \todo Documentation */
    void setAmplitude(double a){
      amplitude_=a;
    }

    /** \todo Documentation */
    double getAmplitude(){
      return amplitude_;
    }
    

    /** \todo Documentation */
    double evalGaussianLikelihood(const Measurement2d_amplitude &x_eval,
				  double* mDist2 = NULL){

      return  Measurement2d::evalGaussianLikelihood(x_eval,mDist2)*exp(-0.5*pow(x_eval.amplitude_-10,2)/2.0)/sqrt(2*PI*2.0);
    }
  
  protected:
    
    double amplitude_;

  };

}

#endif

