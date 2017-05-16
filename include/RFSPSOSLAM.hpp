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
#ifndef RFSPSOSLAM_HPP
#define RFSPSOSLAM_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Timer.hpp"
#include <Eigen/Core>
#include "MurtyAlgorithm.hpp"
#include "PermutationLexicographic.hpp"
#include "RFSPSOParticle.hpp"
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include "OSPA.hpp"
#include "RandomVecMathTools.hpp"

namespace rfs {

  /**
   *  \class RFSPSOoptimizer
   *  \brief Random Finite Set Particle Swarm optimization solver for  feature based SLAM
   *
   *
   *
   *  \tparam RobotProcessModel A robot process model derived from ProcessModel
   *  \tparam MeasurementModel A sensor model derived from MeasurementModel
   *  \author  Felipe Inostroza
   */
  template<class RobotProcessModel, class MeasurementModel>
    class RFSPSOSLAM {

    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      typedef typename RobotProcessModel::TState TPose;
      typedef typename RobotProcessModel::TInput TInput;
      typedef typename MeasurementModel::TLandmark TLandmark;
      typedef typename MeasurementModel::TMeasurement TMeasurement;
      typedef RFSPSOParticle<RobotProcessModel, MeasurementModel> TParticle;
      typedef std::vector<TParticle> TParticleSet;
      /**
       * \brief Configurations for this RFSBatchPSO optimizer
       */
      struct Config {

        int nParticles_; /**< number of particles for the PSO algorithm */

        double ospa_c_; /**< ospa-like c value for calculating set differences and speed*/

        /** The threshold used to determine if a possible meaurement-landmark
         *  pairing is significant to worth considering
         */
        double MeasurementLikelihoodThreshold_;

        double mapFromMeasurementProb_; /**< probability that each measurement will initialize a landmark on map initialization*/

        double w; /**<  Particle Inertia */

        double phi_p; /**<  particle minimum influence   */

        double card_phi_p; /**<  particle minimum influence on cardinality  */

        double phi_g; /**< global minimum influence  */

        double card_phi_g; /**< global minimum influence on cardinality */

        int K; /**< average number of neighbors */

        bool use_global; /**< use global  topology*/

      } config;

      /**
       * Constructor
       */
      RFSPSOSLAM ();

      /** Destructor */
      ~RFSPSOSLAM ();

      /**
       * Add a single measurement
       * @param z The measurement to add, make sure the timestamp is valid.
       */
      void
      addMeasurement (TMeasurement z);

      /**
       * Add a set of measurements
       * @param Z The measurements to add, make sure the timestamp is valid.
       */
      void
      addMeasurement (std::vector<TMeasurement> Z);

      /**
       * Add a single odometry input
       * @param u the odometry input
       */
      void
      addInput (TInput u);

      /**
       * set all the odometry inputs
       * @param U vector containing all odometry inputs
       */
      void
      setInputs (std::vector<TInput> U);

      /**
       * Generate random trajectories based on the motion model of the robot
       */
      void
      initTrajectories ();

      /**
       * Generate random maps based on the already initialized trajectories
       */
      void
      initMaps ();

      /**
       * initialize the particles
       */
      void
      init ();

      /**
       * Get the PSO particles
       * @return pointer to the particle set
       */
      std::vector<TParticle>*
      getParticles ();

      /**
       * Get the best  PSO particle
       * @return pointer to the particle
       */
      TParticle*
      getBestParticle ();

      /**
       * Calculates the measurement likelihood of particle particleIdx at time k
       * @param particleIdx The particle, trajectory and map
       * @param k the time for which to calculate the likelihood
       * @return the measurement likelihood
       */

      double
      rfsMeasurementLikelihood (const int particleIdx, const int k);

      /**
       * Calculates the measurement likelihood of particle particleIdx  including all available times
       * @param particleIdx The particle, trajectory and map
       * @return the measurement likelihood
       */

      double
      rfsMeasurementLikelihood (const int particleIdx);

      /**
       * Evaluate the current likelihood of all the particles
       */
      void
      evaluateLikelihoods ();

      /**
       * Optimize each particle's state using gradient-based algorithms
       */
      void
      optimizeParticles ();

      /**
       * Update the particles velocities according to the PSO algorithm and move the particles accordingly
       */
      void
      moveParticles ();
      /**
       * Reset the topology randomly according to standard PSO see chapter 1 of the "Handbook of Swarm Intelligence" Springer 2011
       *
       */
      void
      updateTopology();

      /***
       * Update the particle map velocity using the ospa inspired motion
       * @param particleIdx The particle index , trajectory and map
       * @param bestParticleIdx The particle index of the best particle in neighborhood
       */
      void
      updateLandmarkVelocity (const int particleIdx, const int bestParticleIdx);

      MeasurementModel *mModelPtr_;
      RobotProcessModel *robotProcessModelPtr_;
    private:

      int nThreads_; /**< Number of threads  */
      int iteration_;
      bool hasImproved_;
      std::vector<TInput> inputs_; /**< vector containing all odometry inputs */
      std::vector<std::vector<TMeasurement> > Z_; /**< vector containing all feature measurements */
      std::vector<TimeStamp> time_;
      std::vector<TParticle> particles_;
      TParticle bestParticle_;

      std::vector< std::vector<int > > topology_;

    };

  //////////////////////////////// Implementation ////////////////////////

  template<class RobotProcessModel, class MeasurementModel>
    typename RFSPSOSLAM<RobotProcessModel, MeasurementModel>::TParticleSet*
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::getParticles () {

      return &particles_;
    }
  template<class RobotProcessModel, class MeasurementModel>
    typename RFSPSOSLAM<RobotProcessModel, MeasurementModel>::TParticle*
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::getBestParticle () {

      return &bestParticle_;
    }
  template<class RobotProcessModel, class MeasurementModel>
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::RFSPSOSLAM () {
      nThreads_ = 1;
      iteration_ = 0;
      hasImproved_ = false;

#ifdef _OPENMP
      nThreads_ = omp_get_max_threads();
#endif
      mModelPtr_ = new MeasurementModel();
      robotProcessModelPtr_ = new RobotProcessModel();

    }

  template<class RobotProcessModel, class MeasurementModel>
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::~RFSPSOSLAM () {

      delete mModelPtr_;
      delete robotProcessModelPtr_;
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::addInput (TInput u) {

      inputs_.push_back(u);
      time_.push_back(u.getTime());
      std::sort(inputs_.begin(), inputs_.begin());
      std::sort(time_.begin(), time_.begin());
      Z_.resize(inputs_.size() + 1);

    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::setInputs (std::vector<TInput> U) {

      inputs_ = U;
      std::sort(inputs_.begin(), inputs_.begin());

      time_.resize(inputs_.size() + 1);
      for (int i = 0; i < inputs_.size(); i++) {
        time_[i + 1] = inputs_[i].getTime();
      }
      Z_.resize(inputs_.size() + 1);
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::addMeasurement (TMeasurement z) {

      TimeStamp zTime = z.getTime();

      auto it = std::lower_bound(time_.begin(), time_.end(), zTime);
      if (*it != zTime) {
        std::cerr << "Measurement time does not match with any of the odometries\n zTime: " << zTime.getTimeAsDouble() << "\n";
        std::exit(1);
      }
      int k = it - time_.begin();
      Z_[k].push_back(z);

    }
  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::updateTopology(){

      //std::cout << "updating topology\n";
      if (config.use_global) {
        if (topology_.size() != config.nParticles_) {
          topology_.resize(config.nParticles_);
          for (int i = 0; i < config.nParticles_; i++) {
            topology_[i].resize(config.nParticles_);
            for (int j = 0; j < config.nParticles_; j++) {
            topology_[i][j] = j;
            }
          }
        }

        return;
      }

    if (topology_.size() !=  config.nParticles_){
      topology_.resize(config.nParticles_);
    }
    for (int i = 0; i < config.nParticles_; i++) {
        topology_[i].resize(1);
        topology_[i][0] = i;
    }

    for (int i = 0; i < config.nParticles_; i++) {
      for (int j =0; j < config.K; j++){
        topology_[floor(drand48()*config.nParticles_)].push_back(i);
      }
    }
  }
  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::moveParticles () {
      if (! hasImproved_){
        updateTopology();
      }


      for (int i = 0; i < config.nParticles_; i++) {

        // find the best particle in the neighborhood defined by current topology
        int best_i = i;
        double best_w = particles_[i].bestLikelihood;
        for (int n = 0; n < topology_[i].size() ;  n++){
          if (particles_[topology_[i][n]].bestLikelihood > best_w){
            best_i = topology_[i][n];
            best_w = particles_[topology_[i][n]].bestLikelihood;
          }
        }


        for (int k = 0; k < particles_[i].inputs.size(); k++) {
          typename TInput::Vec v = particles_[i].inputs_velocity[k].get(), p = particles_[i].inputs[k].get();

          v *= config.w;
          for (int dim = 0; dim < particles_[i].inputs_velocity[k].getNDim(); dim++) {
            v[dim] += config.phi_p * drand48() * (particles_[i].bestInputs[k].get()[dim] - particles_[i].inputs[k].get()[dim]);
            v[dim] += config.phi_g * drand48() * (particles_[best_i].bestInputs[k].get()[dim] - particles_[i].inputs[k].get()[dim]);
          }
          particles_[i].inputs_velocity[k].set(v);

          p += v;
          particles_[i].inputs[k].set(p);
          robotProcessModelPtr_->step(particles_[i].trajectory[k + 1], particles_[i].trajectory[k], particles_[i].inputs[k], time_[k + 1] - time_[k]);
        }

          updateLandmarkVelocity(i , best_i);

        for (int m = 0; m < particles_[i].landmarks.size(); m++) {
          typename TLandmark::Vec lm = particles_[i].landmarks[m].get();
          lm += particles_[i].landmarks_velocity[m].get();
          particles_[i].landmarks[m].set(lm);
        }

      }
      iteration_++;
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::optimizeParticles (){


  }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::updateLandmarkVelocity (const int particleIdx, const int bestParticleIdx) {

      // first use the particle optimum

      unsigned int nM1 = particles_[particleIdx].landmarks.size();
      unsigned int nM2 = particles_[particleIdx].bestLandmarks.size();
      unsigned int nM3 = particles_[bestParticleIdx].bestLandmarks.size();

      // Associate with OSPA and euclidian distance
      rfs::OSPA<TLandmark> ospa_particle(particles_[particleIdx].landmarks, particles_[particleIdx].bestLandmarks, config.ospa_c_, 1,
                                         &(RandomVecMathTools<TLandmark>::euclidianDist));
      rfs::OSPA<TLandmark> ospa_global(particles_[particleIdx].landmarks, particles_[bestParticleIdx].bestLandmarks, config.ospa_c_, 1,
                                       &(RandomVecMathTools<TLandmark>::euclidianDist));

      boost::shared_array<int> SolnArr_particle = ospa_particle.getOptAssignment();
      boost::shared_array<int> SolnArr_global = ospa_global.getOptAssignment();
      std::vector<bool> toadd_particle(nM2, true);
      std::vector<bool> toadd_global(nM3, true);
      std::vector<bool> todelete(nM1);
/*
      std::cout << "\n\ncurrentmap: " ;
      for (int i = 0; i < nM1; i++) {
        std::cout << particles_[particleIdx].landmarks[i][0] << "   ";
      }
      std::cout << "\n";
      std::cout << "best   map: " ;
           for (int i = 0; i < nM3; i++) {
             std::cout << bestParticle_.landmarks[i][0] << "   ";
           }
           std::cout << "\n";

           std::cout << "best  dist: " ;
                for (int i = 0; i < std::min(nM3,nM1); i++) {
                  std::cout << RandomVecMathTools<TLandmark>::euclidianDist(bestParticle_.landmarks[i],particles_[particleIdx].landmarks[i]) << "   ";
                }
                std::cout << "\n";

      std::cout << "nM1:  " << nM1 <<"  nM2: "<< nM2 << "  nM3:  " << nM3<<"  assoc: " ;*/
      for (int i = 0; i < nM1; i++) {
        typename TLandmark::Vec v = particles_[particleIdx].landmarks_velocity[i].get();
        v *= config.w;
        int j_p = SolnArr_particle[i];

        if (j_p >= nM2) {
          if (drand48() < config.card_phi_p)
            todelete[i] = true;
        }
        else {
          if (drand48() < 0.01)
            todelete[i] = true;
          toadd_particle[j_p] = false;
          for (int dim = 0; dim < particles_[particleIdx].landmarks_velocity[i].getNDim(); dim++) {
            v[dim] += config.phi_p * drand48() * (particles_[particleIdx].bestLandmarks[j_p].get()[dim] - particles_[particleIdx].landmarks[i].get()[dim]);
          }
        }
        int j_g = SolnArr_global[i];
        //std::cout << "  " << j_g;
        if (j_g >= nM3) {
          if (drand48() < config.card_phi_g)
            todelete[i] = true;
        }
        else {
          toadd_global[j_g] = false;
          for (int dim = 0; dim < particles_[particleIdx].landmarks_velocity[i].getNDim(); dim++) {
            v[dim] += config.phi_g * drand48() * (particles_[bestParticleIdx].bestLandmarks[j_g].get()[dim] - particles_[particleIdx].landmarks[i].get()[dim]);
          }
        }

        particles_[particleIdx].landmarks_velocity[i].set(v);
      }
      //std::cout << "\n";


      //
      for (int i = nM1 - 1; i >= 0; i--) {
        if (todelete[i]) {
          particles_[particleIdx].landmarks[i] = particles_[particleIdx].landmarks[particles_[particleIdx].landmarks.size() - 1];
          particles_[particleIdx].landmarks_velocity[i] = particles_[particleIdx].landmarks_velocity[particles_[particleIdx].landmarks_velocity.size() - 1];
          particles_[particleIdx].landmarks.pop_back();
          particles_[particleIdx].landmarks_velocity.pop_back();
        }
      }
      for (int i = 0; i < nM2; i++) {
        if (toadd_particle[i]) {
          if (drand48() < config.card_phi_p) {
            particles_[particleIdx].landmarks.push_back(particles_[particleIdx].bestLandmarks[i]);
            particles_[particleIdx].landmarks_velocity.push_back(particles_[particleIdx].bestLandmarks_velocity[i]);
          }
        }
      }
      for (int i = 0; i < nM3; i++) {
        if (toadd_global[i]) {
          if (drand48() < config.card_phi_g) {
            particles_[particleIdx].landmarks.push_back(particles_[bestParticleIdx].bestLandmarks[i]);
            particles_[particleIdx].landmarks_velocity.push_back(particles_[bestParticleIdx].bestLandmarks_velocity[i]);
          }
        }
      }
/*
      std::cout << "new    map: " ;
            for (int i = 0; i < particles_[particleIdx].landmarks.size(); i++) {
              std::cout << particles_[particleIdx].landmarks[i][0] << "   ";
            }
            std::cout << "\n";
*/


      //std::cout << "newSize: " << particles_[particleIdx].landmarks.size() << "\n";
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::addMeasurement (std::vector<TMeasurement> Z) {

      for (int i = 0; i < Z.size(); i++) {
        this->addMeasurement(Z[i]);
      }
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::init () {
      particles_.resize(config.nParticles_);

      // initialize the number of neighbors distribution
      double c=0;

      for(int n = 0 ; n < config.K; n++){

      }

      initTrajectories();
      initMaps();
      updateTopology();
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::initTrajectories () {

      for (int i = 0; i < config.nParticles_; i++) {

        particles_[i].trajectory.resize(inputs_.size() + 1);
        particles_[i].inputs.resize(inputs_.size());
        particles_[i].inputs_velocity.resize(inputs_.size());
        for (int k = 0; k < inputs_.size(); k++) {
          TimeStamp dT = inputs_[k].getTime() - particles_[i].trajectory[k].getTime();
          robotProcessModelPtr_->sample(particles_[i].trajectory[k + 1], particles_[i].trajectory[k], inputs_[k], dT, false, true, &particles_[i].inputs[k]);
        }
        particles_[i].bestTrajectory = particles_[i].trajectory;
        particles_[i].bestInputs = particles_[i].inputs;

        particles_[i].bestInputs_velocity.resize(inputs_.size());

      }
    }
  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::initMaps () {

      for (int i = 0; i < config.nParticles_; i++) {

        for (int k = 0; k < particles_[i].trajectory.size(); k++) {
          for (int nz = 0; nz < Z_[k].size(); nz++) {

            if (drand48() < config.mapFromMeasurementProb_) {

              TLandmark lm;
              this->mModelPtr_->inverseMeasure(particles_[i].trajectory[k], Z_[k][nz], lm);
              particles_[i].landmarks.push_back(lm);
            }
          }
        }

        particles_[i].landmarks_velocity.resize(particles_[i].landmarks.size());
        particles_[i].bestLandmarks_velocity.resize(particles_[i].landmarks.size());
        particles_[i].bestLandmarks = particles_[i].landmarks;

      }

    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::evaluateLikelihoods () {
      int besti = 0;
#pragma omp parallel for
      for (int i = 0; i < config.nParticles_; i++) {
        particles_[i].currentLikelihood = rfsMeasurementLikelihood(i);
        //std::cout << "likeli: " << particles_[i].currentLikelihood << "\n";

        if (particles_[i].currentLikelihood > particles_[i].bestLikelihood) {
          particles_[i].bestLikelihood = particles_[i].currentLikelihood;
          particles_[i].bestTrajectory = particles_[i].trajectory;
          particles_[i].bestInputs = particles_[i].inputs;
          particles_[i].bestInputs_velocity = particles_[i].inputs_velocity;
          particles_[i].bestLandmarks = particles_[i].landmarks;
          particles_[i].bestLandmarks_velocity = particles_[i].landmarks_velocity;
        }
#pragma omp critical
        {
          if (particles_[i].currentLikelihood > particles_[besti].currentLikelihood) {
            besti = i;
          }
        }
      }
      hasImproved_ = false;
      if (particles_[besti].currentLikelihood > bestParticle_.currentLikelihood) {
        bestParticle_ = particles_[besti];
        hasImproved_ = true;
        std::cout << "new best particle :   " << bestParticle_.currentLikelihood << "\n";
      }

    }

  template<class RobotProcessModel, class MeasurementModel>
    double
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::rfsMeasurementLikelihood (const int particleIdx) {
      double l =  log(rfsMeasurementLikelihood(particleIdx, 0));
      TimeStamp dT;
      for (int k = 1; k < particles_[particleIdx].trajectory.size(); k++) {
        l += log(rfsMeasurementLikelihood(particleIdx, k));
        dT = time_[k] - time_[k - 1];
        l += log(robotProcessModelPtr_->likelihood(particles_[particleIdx].trajectory[k], particles_[particleIdx].trajectory[k-1], inputs_[k - 1], dT));
      }
      //std::cout << "likelihood   " << l << "\n";
      return l;
    }

  template<class RobotProcessModel, class MeasurementModel>
    double
    RFSPSOSLAM<RobotProcessModel, MeasurementModel>::rfsMeasurementLikelihood (const int particleIdx, const int k) {

      const int i = particleIdx;
      TPose* pose = &particles_[particleIdx].trajectory[k];
      const int nZ = this->Z_[k].size();
      const unsigned int mapSize = particles_[particleIdx].landmarks.size();
      // Find map points within field of view and their probability of detection
      std::vector<unsigned int> lmInFovIdx;
      std::vector<double> lmInFovPd;
      std::vector<int> landmarkCloseToSensingLimit;

      for (unsigned int m = 0; m < mapSize; m++) {

        bool isCloseToSensingLimit = false;

        TLandmark* lm = &this->particles_[particleIdx].landmarks[m];

        bool isClose;
        double Pd = this->mModelPtr_->probabilityOfDetection(*pose, *lm, isCloseToSensingLimit);

        if (Pd > 0) {
          lmInFovIdx.push_back(m);
          lmInFovPd.push_back(Pd);

          landmarkCloseToSensingLimit.push_back(isCloseToSensingLimit);
        }

      }
      const unsigned int nM = lmInFovIdx.size();

      // If map is empty everything must be a false alarm

      double clutter[nZ];
      for (int n = 0; n < nZ; n++) {
        clutter[n] = this->mModelPtr_->clutterIntensity(this->Z_[k][n], nZ);
      }

      if (nM == 0) {
        double l = 1;
        for (int n = 0; n < nZ; n++) {
          l *= clutter[n];
        }
        return l;
      }

      TLandmark* evalPt;
      TLandmark evalPt_copy;
      TMeasurement expected_z;

      double md2; // Mahalanobis distance squared

      // Create and fill in likelihood table (nM x nZ)
      double** L;
      CostMatrixGeneral likelihoodMatrix(L, nM, nZ);

      for (int m = 0; m < nM; m++) {

        evalPt = &this->particles_[i].landmarks[lmInFovIdx[m]]; // get location of m
        evalPt_copy = *evalPt; // so that we don't change the actual data //
        evalPt_copy.setCov(MeasurementModel::TLandmark::Mat::Zero()); //
        this->mModelPtr_->measure(*pose, evalPt_copy, expected_z); // get expected measurement for m
        double Pd = lmInFovPd[m]; // get the prob of detection of m

        for (int n = 0; n < nZ; n++) {

          // calculate measurement likelihood with detection statistics
          L[m][n] = expected_z.evalGaussianLikelihood(this->Z_[k][n], &md2) * Pd; // new line
          if (L[m][n] < config.MeasurementLikelihoodThreshold_) {
            L[m][n] = 0;
          }
        }
      }

      // Partition the Likelihood Table and turn into a log-likelihood table
      int nP = likelihoodMatrix.partition();
      double l = 1;
      double const BIG_NEG_NUM = -1000; // to represent log(0)

      // Go through each partition and determine the likelihood
      for (int p = 0; p < nP; p++) {

        double partition_likelihood = 0;

        unsigned int nCols, nRows;
        double** Cp;
        unsigned int* rowIdx;
        unsigned int* colIdx;

        bool isZeroPartition = !likelihoodMatrix.getPartitionSize(p, nRows, nCols);
        bool useMurtyAlgorithm = true;
        if (nRows + nCols <= 8 || isZeroPartition)
          useMurtyAlgorithm = false;

        isZeroPartition = !likelihoodMatrix.getPartition(p, Cp, nRows, nCols, rowIdx, colIdx, useMurtyAlgorithm);

        if (isZeroPartition) { // all landmarks in this partition are mis-detected. All measurements are outliers

          partition_likelihood = 1;
          for (int r = 0; r < nRows; r++) {
            partition_likelihood *= lmInFovPd[rowIdx[r]];
          }

          for (int c = 0; c < nCols; c++) {
            partition_likelihood *= clutter[colIdx[c]];
          }

        }
        else {
          // turn the matrix into a log likelihood matrix with detection statistics,
          // and fill in the extended part of the partition

          for (int r = 0; r < nRows; r++) {
            for (int c = 0; c < nCols; c++) {
              if (Cp[r][c] == 0)
                Cp[r][c] = BIG_NEG_NUM;
              else {
                Cp[r][c] = log(Cp[r][c]);
                if (Cp[r][c] < BIG_NEG_NUM)
                  Cp[r][c] = BIG_NEG_NUM;
              }
            }
          }

          if (useMurtyAlgorithm) { // use Murty's algorithm

            // mis-detections
            for (int r = 0; r < nRows; r++) {
              for (int c = nCols; c < nRows + nCols; c++) {
                if (r == c - nCols)
                  Cp[r][c] = log(1 - lmInFovPd[rowIdx[r]]);
                else
                  Cp[r][c] = BIG_NEG_NUM;
              }
            }

            // clutter
            for (int r = nRows; r < nRows + nCols; r++) {
              for (int c = 0; c < nCols; c++) {
                if (r - nRows == c)
                  Cp[r][c] = log(clutter[colIdx[c]]);
                else
                  Cp[r][c] = BIG_NEG_NUM;
              }
            }

            // the lower right corner
            for (int r = nRows; r < nRows + nCols; r++) {
              for (int c = nCols; c < nRows + nCols; c++) {
                Cp[r][c] = 0;
              }
            }

            Murty murtyAlgo(Cp, nRows + nCols);
            Murty::Assignment a;
            partition_likelihood = 0;
            double permutation_log_likelihood = 0;
            murtyAlgo.setRealAssignmentBlock(nRows, nCols);
            for (int k = 0; k < 200; k++) {
              int rank = murtyAlgo.findNextBest(a, permutation_log_likelihood);
              if (rank == -1 || permutation_log_likelihood < BIG_NEG_NUM)
                break;
              partition_likelihood += exp(permutation_log_likelihood);
            }

          }
          else { // use lexicographic ordering

            partition_likelihood = 0;
            double permutation_log_likelihood = 0;

            uint o[nRows + nCols];

            PermutationLexicographic pl(nRows, nCols, true);
            unsigned int nPerm = pl.next(o);
            while (nPerm != 0) {
              permutation_log_likelihood = 0;
              for (int a = 0; a < nRows; a++) {
                if (o[a] < nCols) { // detection
                  permutation_log_likelihood += Cp[a][o[a]];
                }
                else { // mis-detection
                  permutation_log_likelihood += log(1 - lmInFovPd[rowIdx[a]]);
                }
              }
              for (int a = nRows; a < nRows + nCols; a++) { // outliers
                if (o[a] < nCols) {
                  permutation_log_likelihood += log(clutter[colIdx[o[a]]]);
                }
              }
              partition_likelihood += exp(permutation_log_likelihood);
              nPerm = pl.next(o);
            }

          } // End lexicographic ordering

        } // End non zero partition

        l *= partition_likelihood;

      } // End partitions

      return (l / this->mModelPtr_->clutterIntensityIntegral(nZ));
    }

}

#endif
