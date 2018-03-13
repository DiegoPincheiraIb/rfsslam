/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2014, Keith Leung, Felipe Inostroza
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

#ifndef RFSPSOPARTICLE_HPP
#define RFSPSOPARTICLE_HPP

#include "TimeStamp.hpp"


namespace rfs {
  /**
   *  \class RFSPSOParticle
   *  \brief Status for a batch optimization
   *
   *  This class stores a solution to SLAM, intended for use in batch optimizers.
   *  It can calculate its own Measurement likelihood.
   *  \todo Add a function to calculate the gradient for gradient decent based algorithms.
   *
   *
   *  \author  Felipe Inostroza
   */
  template<class RobotProcessModel, class MeasurementModel>
    class RFSPSOParticle
  {
  public:

    std::vector<typename RobotProcessModel::TInput> inputs , bestInputs;
    std::vector<typename RobotProcessModel::TState> trajectory , bestTrajectory;
    std::vector<typename MeasurementModel::TLandmark> landmarks , bestLandmarks;

    std::vector<typename RobotProcessModel::TInput> inputs_velocity, bestInputs_velocity;
    std::vector<typename RobotProcessModel::TState> trajectory_velocity , bestTrajectory_velocity;
    std::vector<typename MeasurementModel::TLandmark> landmarks_velocity, bestLandmarks_velocity;

    double currentLikelihood = -std::numeric_limits<double>::infinity()
        , bestLikelihood = -std::numeric_limits<double>::infinity();



  };




}



#endif
