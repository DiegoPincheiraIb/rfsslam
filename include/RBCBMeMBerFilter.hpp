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
#ifndef RBCBMEMBERFILTER_HPP
#define RBCBMEMBERFILTER_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Timer.hpp"
#include <Eigen/Core>
#include "MurtyAlgorithm.hpp"
#include "PermutationLexicographic.hpp"
#include "GMBernoulliComponent.hpp"
#include "ParticleFilter.hpp"
#include "KalmanFilter.hpp"
#include <math.h>
#include <vector>
#include <stdio.h>

namespace rfs {
/**
 *  \class RBCBMeMBerFilter
 *  \brief Rao-Blackwellized Cardinality Balanced Multi Target Multi Bernoulli Filter class
 *  
 *  This class implements the Rao-Bloackwellized Cardinality Balanced Multi Target Multi Bernoulli
 *  filter. The constructor of this class will internally instantiate the 
 *  process model for both the robot and landmarks, the measurement model, 
 *  and the Kalman filter. Users have access to these through pointers that 
 *  can be obtained by calling the appropraite get function.
 *
 *  \tparam RobotProcessModel A robot process model derived from ProcessModel
 *  \tparam LmkProcessModel A landmark process model derived from ProcessModel
 *  \tparam MeasurementModel A sensor model derived from MeasurementModel
 *  \tparam KalmanFilter A Kalman filter that uses LmkProcessModel and MeasurementModel
 *  \author Keith Leung, Felipe Inostroza
 */
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  class RBCBMeMBerFilter : public ParticleFilter<RobotProcessModel, MeasurementModel,
      GMMultiBernoulli<typename MeasurementModel::TLandmark> >
  {
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef typename RobotProcessModel::TState TPose;
    typedef typename RobotProcessModel::TInput TInput;
    typedef typename MeasurementModel::TLandmark TLandmark;
    typedef typename MeasurementModel::TMeasurement TMeasurement;
    typedef typename GaussianMixture<TLandmark>::Gaussian TGaussian;

    /** 
     * \brief Configurations for this RBCBMeMBerFilter
     */
    struct Config
    {

      /**  New birth Tracks are set this probability of existence */
      double birthTrackPE_;

      /**  New Gaussians are only created during map update if the innovation mahalanobis distance
       is less than this threshold */
      double newGaussianCreateInnovMDThreshold_;

      /**  number of map states to use for evaluating particle weight
       0 => empty-set strategy,
       1 => single-feature strategy,
       >1 => multi-feature strategy
       */
      int importanceWeightingEvalPointCount_;

      /** The weight above which a Gaussian's mean is considered as a evaluation point for particle importance factor */
      double importanceWeightingEvalPointTrackPE_;

      /** The Mahalanobis distance threshold used to determine if a possible meaurement-landmark
       *  pairing is significant to worth considering
       */
      double importanceWeightingMeasurementLikelihoodMDThreshold_;

      /** Gaussian merging Mahalanobis distance threshold */
      double gaussianMergingThreshold_;

      /** Gaussian merging covariance inflation factor */
      double gaussianMergingCovarianceInflationFactor_;

      /** Gaussian pruning weight threshold, below which Gaussians are eliminated from a Gaussian mixture */
      double gaussianPruningThreshold_;

      /** Track pruning probability of existence threshold, below which tracks are eliminated from a MultiBernoulli distribution */
      double trackPruningThreshold_;

      /** Minimum number of updates between resampling of particles*/
      int minUpdatesBeforeResample_;

      /** Use the particle weighting strategy from Single-cluster PHD Filtering by Lee, et. al. */
      bool useClusterProcess_;

    } config;

    /**
     * \brief Elapsed timing information
     */
    struct TimingInfo
    {
      long long predict_wall;
      long long predict_cpu;
      long long mapUpdate_wall;
      long long mapUpdate_cpu;
      long long mapUpdate_kf_wall;
      long long mapUpdate_kf_cpu;
      long long particleWeighting_wall;
      long long particleWeighting_cpu;
      long long mapMerge_wall;
      long long mapMerge_cpu;
      long long mapPrune_wall;
      long long mapPrune_cpu;
      long long particleResample_wall;
      long long particleResample_cpu;
    } timingInfo_;

    /**
     * Constructor
     * \param n number of particles
     */
    RBCBMeMBerFilter(int n);

    /** Destructor */
    ~RBCBMeMBerFilter();
    /**
     * Get the landmark process model
     * \return pointer to the landmark process model
     */
    LmkProcessModel* getLmkProcessModel();

    /**
     * Predict the robot trajectory using the lastest odometry data
     * \param[in] u input
     * \param[in] Timestep timestep, which the motion model may or may not use;
     * \param[in] useModelNoise use the additive noise for the process model
     * \param[in] useInputNoise use the noise in the input
     */
    void predict(TInput u, TimeStamp const &dT, bool useModelNoise = true, bool useInputNoise = false);
    /**
     * Update the map, calculate importance weighting, sample if necessary, and
     * create new birth Tracks.
     * \param[in] Z set of measurements to use for the update, placed in a std vector, which
     * gets cleared after the function call.
     */
    void update(std::vector<TMeasurement> &Z);

    /**
     * Get the pointer to the Kalman Filter used for updating the map
     * \return pointer to the Kalman Filter
     */
    KalmanFilter* getKalmanFilter();

    /** Function for initiating particles during startup
     *  \param[in] i particle index
     *  \param[in] p particle pose
     */
    void setParticlePose(int i, TPose &p);
    /**
     * Get elapsed timing information for various steps of the filter
     */
    TimingInfo* getTimingInfo();

    /**
     * Get the number of tracks of a particle
     * \param[in] i particle index
     * \return size if index is valid, else -1
     */
    int getTrackNum(int i);
    /**
     * Get the track m  of particle i
     * \param[in] i particle index
     * \param[in] m Track index
     * \param[out] track The track
     * \return false if the indices specified are invalid
     */
    bool getTrack(const int i, const int m, GMBernoulliComponent<TLandmark> &track);

  private:

    int nThreads_; /**< Number of threads for running map update for all particles */

    std::vector<KalmanFilter> kfs_; /**< Kalman filters (one for each thread)*/
    LmkProcessModel *lmkModelPtr_; /**< pointer to landmark process model */
    double threshold_mahalanobisDistance2_mapUpdate_; /**< Mahalanobis distance squared threshold, above which Kalman Filter update will not occur*/
    /** indices of unused measurement for each particle for creating birth Gaussians */
    std::vector<std::vector<unsigned int> > unused_measurements_;

    unsigned int nUpdatesSinceResample; /**< Number of updates performed since the last resmaple */

    Timer timer_predict_; /**< Timer for prediction step */
    Timer timer_mapUpdate_; /**< Timer for map update */
    Timer timer_mapUpdate_kf_; /**<Timer for the Kalman filter part of map update */
    Timer timer_particleWeighting_; /**< Timer for particle weighting */
    Timer timer_mapMerge_; /**< Timer for map merging */
    Timer timer_mapPrune_; /**< Timer for map pruning */
    Timer timer_particleResample_; /**<Timer for particle resampling */

    /**
     * Add birth tracks for each particle's map using unused_measurements_
     */
    void addBirthTracks();

    /**
     * Update the map with the measurements in measurements_
     * Existing landmark tracks with probability of detection > 0 will have their probability of existence
     * reduced to account for missed detection.
     * For every landmark-measurement pair with probability of detection > 0,
     * a new landmark track will be created.
     * \param[in] particleIdx index of the particle for which the map will be updated
     */
    void updateMap(const uint particleIdx);

    /**
     * Importance weighting. Overrides the abstract function in ParticleFilter
     * \param[in] idx particle index
     */
    void importanceWeighting(const uint idx);

    /**
     * Random Finite Set measurement likelihood evaluation
     * \brief The current measurements in measurements_ are used to determine the
     * RFS measurement likelihood given a set of landmarks
     * \param[in] particleIdx particle for which the likelihood is calcuated
     * \param[in] evalPtIdx indices of evaluation points in maps_[particleIdx]
     * \param[in] evalPtPd probability of detection of evaluation point
     * \return measurement likelihood
     */
    double rfsMeasurementLikelihood(const int particleIdx, std::vector<unsigned int> &evalPtIdx,
                                    std::vector<double> &evalPtPd);

    /**
     * Random Finite Set map likelihood evaluation
     * \brief The current map estimate is  used to determine the
     * RFS likelihood of the set set of landmarks provided in evalPtIdx
     * \param[in] particleIdx particle for which the likelihood is calcuated
     * \param[in] evalPtIdx indices of evaluation points in maps_[particleIdx]
     * \return measurement likelihood
     */

    double rfsMapLikelihood(const int particleIdx, std::vector<unsigned int> &evalPtIdx);
    /**
     * Random Finite Set map likelihood evaluation for the unupdated map, this uses the legacy tracks and recalculates the previous weights.
     *\brief The previous map estimate is  used to determine the
     * RFS likelihood of the set set of landmarks provided in evalPtIdx
     * \param[in] particleIdx particle for which the likelihood is calcuated
     * \param[in] evalPtIdx indices of evaluation points in maps_[particleIdx]
     * \return measurement likelihood
     */

    double rfsPrevMapLikelihood(const int particleIdx, std::vector<unsigned int> &evalPtIdx);

    /**
     * Calculate the sum of all permutations of measurement likelihood from a likelihood
     * table generated from within rfsMeasurementLikelihood
     * \param[in] likelihoodTab likelihood table generated within rfsMeasurementLikelihood
     * \param[in] Z_NoClutter A vector of measurement indices (columns) to consider in the likelihoodTab
     * \return sum of all permutations from the given likelihood table
     */
    double rfsMeasurementLikelihoodPermutations(std::vector<double*> &likelihoodTab, std::vector<int> &Z_NoClutter);

  };

////////// Implementation //////////

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  LmkProcessModel* RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getLmkProcessModel() {
    return lmkModelPtr_;
  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::predict(
      TInput u, TimeStamp const &dT, bool useModelNoise, bool useInputNoise) {

    timer_predict_.resume();

    // Add birth tracks using pose before prediction
    addBirthTracks();

    // propagate particles
    this->propagate(u, dT, useModelNoise, useInputNoise);

    // propagate landmarks
    for (int i = 0; i < this->nParticles_; i++) {
      for (int m = 0; m < this->particleSet_[i]->getData()->tracks_.size(); m++) {
        for (int g = 0; g < this->particleSet_[i]->getData()->tracks_[m].getGaussianCount(); g++) {
          TLandmark *plm;
          this->particleSet_[i]->getData()->tracks_[m].getGaussian(g, plm);
          lmkModelPtr_->staticStep(*plm, *plm, dT);
        }
      }
    }

    timer_predict_.stop();
  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::addBirthTracks() {

    for (int i = 0; i < this->nParticles_; i++) {

      while (unused_measurements_[i].size() > 0) {

        // get measurement
        int unused_idx = unused_measurements_[i].back();
        TMeasurement unused_z = this->measurements_[unused_idx];
        unused_measurements_[i].pop_back();

        // use inverse measurement model to get landmark
        TPose robot_pose;
        TLandmark landmark_pos;
        this->particleSet_[i]->getPose(robot_pose);
        this->pMeasurementModel_->inverseMeasure(robot_pose, unused_z, landmark_pos);

        // add birth landmark to Gaussian mixture (last param = true to allocate mem)
        this->particleSet_[i]->getData()->tracks_.resize(this->particleSet_[i]->getData()->tracks_.size() + 1);
        this->particleSet_[i]->getData()->tracks_[this->particleSet_[i]->getData()->tracks_.size() - 1].addGaussian(
            &landmark_pos, 1, true);
        this->particleSet_[i]->getData()->tracks_[this->particleSet_[i]->getData()->tracks_.size() - 1].setP(
            config.birthTrackPE_);

      }

    }

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::RBCBMeMBerFilter(int n) :
      ParticleFilter<RobotProcessModel, MeasurementModel, GMMultiBernoulli<typename MeasurementModel::TLandmark> >(n) {
    nThreads_ = 1;

#ifdef _OPENMP
    nThreads_ = omp_get_max_threads();
#endif

    lmkModelPtr_ = new LmkProcessModel;
    kfs_ = std::vector<KalmanFilter>(nThreads_, KalmanFilter(lmkModelPtr_, this->getMeasurementModel()));

    for (int i = 0; i < n; i++) {
      // printf("Creating map structure for particle %d\n", i);
      this->particleSet_[i]->setData(
          boost::shared_ptr<GMMultiBernoulli<TLandmark> >(new GMMultiBernoulli<TLandmark>()));
      unused_measurements_.push_back(std::vector<unsigned int>());
    }

    config.birthTrackPE_ = 0.5;
    config.gaussianMergingThreshold_ = 0.5;
    config.gaussianMergingCovarianceInflationFactor_ = 1.5;
    config.gaussianPruningThreshold_ = 0.2;
    config.trackPruningThreshold_ = 0.001;
    config.importanceWeightingEvalPointCount_ = 8;
    config.importanceWeightingMeasurementLikelihoodMDThreshold_ = 3.0;
    config.newGaussianCreateInnovMDThreshold_ = 0.2;
    config.minUpdatesBeforeResample_ = 1;

    nUpdatesSinceResample = 0;

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  KalmanFilter* RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getKalmanFilter() {
    return &(kfs_[0]);
  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::setParticlePose(int i,
                                                                                                             TPose &p) {

    this->particleSet_[i]->setPose(p);

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  bool RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getTrack(
      const int i, const int m, GMBernoulliComponent<TLandmark> &track) {

    int sz = getTrackNum(i);
    if (sz == -1 || (m < 0) || (m >= sz)) {
      return false;
    }
    TLandmark *plm;
    this->particleSet_[i]->getData()->tracks_.at(m).copyTo(&track);
    return true;
  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  typename RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::TimingInfo* RBCBMeMBerFilter<
      RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getTimingInfo() {

    timer_predict_.elapsed(timingInfo_.predict_wall, timingInfo_.predict_cpu);
    timer_mapUpdate_.elapsed(timingInfo_.mapUpdate_wall, timingInfo_.mapUpdate_cpu);
    timer_mapUpdate_kf_.elapsed(timingInfo_.mapUpdate_kf_wall, timingInfo_.mapUpdate_kf_cpu);
    timer_particleWeighting_.elapsed(timingInfo_.particleWeighting_wall, timingInfo_.particleWeighting_cpu);
    timer_mapMerge_.elapsed(timingInfo_.mapMerge_wall, timingInfo_.mapMerge_cpu);
    timer_mapPrune_.elapsed(timingInfo_.mapPrune_wall, timingInfo_.mapPrune_cpu);
    timer_particleResample_.elapsed(timingInfo_.particleResample_wall, timingInfo_.particleResample_cpu);

    return &timingInfo_;
  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  int RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getTrackNum(int i) {

    if (i >= 0 && i < this->nParticles_)
      return (this->particleSet_[i]->getData()->tracks_.size());
    else
      return -1;
  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::~RBCBMeMBerFilter() {

    for (int i = 0; i < this->nParticles_; i++) {
      this->particleSet_[i]->deleteData();
    }
    delete lmkModelPtr_;

  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::update(
      std::vector<TMeasurement> &Z) {

    nUpdatesSinceResample++;

    this->setMeasurements(Z); // Z gets cleared after this call, measurements now stored in this->measurements_
    if (this->measurements_.size() == 0)
      return;

    const unsigned int startIdx = 0;
    const unsigned int stopIdx = this->nParticles_;
    const unsigned int nZ = this->measurements_.size();

    // make sure any setting changes to the Kalman Filter are set for all threads
    if (nThreads_ > 1) {
      for (int j = 1; j < nThreads_; j++) {
        kfs_[j] = kfs_[0];
      }
    }

    if (nThreads_ > 1) {
      timer_mapUpdate_.resume();
    }
#pragma omp parallel
    {
      ////////// Map Update //////////
      if (nThreads_ == 1) {
        timer_mapUpdate_.resume();
      }
#pragma omp for
      for (unsigned int i = startIdx; i < stopIdx; i++) {
        updateMap(i);
      }
      if (nThreads_ == 1) {
        timer_mapUpdate_.stop();
      }

      ////////// Particle Weighintg //////////
      if (!config.useClusterProcess_) {
        if (nThreads_ == 1) {
          timer_particleWeighting_.resume();
        }
#pragma omp for
        for (unsigned int i = startIdx; i < stopIdx; i++) {
          importanceWeighting(i);
        }
        if (nThreads_ == 1) {
          timer_particleWeighting_.stop();
        }
      }

      //////////// Merge and prune //////////
      if (nThreads_ == 1) {
        timer_mapMerge_.resume();
      }
#pragma omp for
      for (unsigned int i = startIdx; i < stopIdx; i++) {
        this->particleSet_[i]->getData()->mergeGM(config.gaussianMergingThreshold_,
                                                  config.gaussianMergingCovarianceInflationFactor_);
      }
      if (nThreads_ == 1) {
        timer_mapMerge_.stop();
      }

      if (nThreads_ == 1) {
        timer_mapPrune_.resume();
      }
#pragma omp for
      for (unsigned int i = startIdx; i < stopIdx; i++) {
        this->particleSet_[i]->getData()->pruneTracks(config.trackPruningThreshold_);
        this->particleSet_[i]->getData()->pruneGM(config.gaussianPruningThreshold_);
      }
      if (nThreads_ == 1) {
        timer_mapPrune_.stop();
      }
    }
    if (nThreads_ > 1) {
      timer_mapUpdate_.stop();
    }

    //////////// Particle resampling //////////
    timer_particleResample_.resume();
    bool resampleOccured = false;
    if (nUpdatesSinceResample >= config.minUpdatesBeforeResample_) {
      resampleOccured = this->resample();
    }

    if (resampleOccured) {
      nUpdatesSinceResample = 0;
    }
    else {
      this->normalizeWeights();
    }
    timer_particleResample_.stop();

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::updateMap(
      const uint particleIdx) {

    const unsigned int i = particleIdx;
    const unsigned int nZ = this->measurements_.size();

    int threadnum = 0;
#ifdef _OPENMP
    threadnum = omp_get_thread_num();
#endif

    //---------- 1. setup / book-keeping ----------

    const unsigned int nM = this->particleSet_[i]->getData()->tracks_.size();
    unused_measurements_[i].clear();
    if (nM == 0) { // No existing landmark case -> flag all measurements as unused and go to next particles
      for (int z = 0; z < nZ; z++) {
        unused_measurements_[i].push_back(z);
      }
      return; // goto next particle
    }
    double Pd[nM];
    int landmarkCloseToSensingLimit[nM];

    //----------  2. Kalman Filter map update ----------

    if (nThreads_ == 1)
      timer_mapUpdate_kf_.resume();
    const TPose *pose = this->particleSet_[i]->getPose();
    threshold_mahalanobisDistance2_mapUpdate_ = config.newGaussianCreateInnovMDThreshold_
        * config.newGaussianCreateInnovMDThreshold_;

    double innovationLikelihood;
    double innovationMahalanobisDist2;
    this->particleSet_[i]->getData()->tracks_.resize(nM + nZ);
    TLandmark lmNew;
    // determine the probability of detection
    for (unsigned int m = 0; m < nM; m++) {

      bool isCloseToSensingLimit = false;
      Pd[m] = 0;
      for (unsigned int g = 0; g < this->particleSet_[i]->getData()->tracks_[m].getGaussianCount(); g++) {
        TLandmark* lm = this->particleSet_[i]->getData()->tracks_[m].getGaussian(g);

        bool isClose;
        Pd[m] += this->particleSet_[i]->getData()->tracks_[m].getWeight(g)
            * this->pMeasurementModel_->probabilityOfDetection(*pose, *lm, isClose);
        isCloseToSensingLimit = isCloseToSensingLimit || isClose;
      } // gaussian in track m for loop
      landmarkCloseToSensingLimit[m] = isCloseToSensingLimit;

      // update weights for legacy tracks
      double prevP = this->particleSet_[i]->getData()->tracks_[m].getP();
      this->particleSet_[i]->getData()->tracks_[m].setP(prevP * (1 - Pd[m]) / (1 - prevP * Pd[m]));
      this->particleSet_[i]->getData()->tracks_[m].setPrevP(prevP);

    }
    // generate updated tracks
    for (unsigned int z = 0; z < nZ; z++) {
      double num = 0; // numerator in eq 56 in CBMeMBer paper
      double denom = 0; // denominator in eq 56 in CBMeMBer paper

      for (unsigned int m = 0; m < nM; m++) {

        double p = this->particleSet_[i]->getData()->tracks_[m].getPrevP(); // previous probability of existence of track m (changed earlier in this method)
        double rho = 0; // from eq 56 in CBMeMBer paper

        for (unsigned int g = 0; g < this->particleSet_[i]->getData()->tracks_[m].getGaussianCount(); g++) {
          TLandmark* lm = this->particleSet_[i]->getData()->tracks_[m].getGaussian(g);

          if (kfs_[threadnum].correct(*pose, this->measurements_[z], *lm, lmNew, &innovationLikelihood,
                                      &innovationMahalanobisDist2)) {

            rho += innovationLikelihood * this->particleSet_[i]->getData()->tracks_[m].getWeight(g) * Pd[m];
            this->particleSet_[i]->getData()->tracks_[nM + z].addGaussian(
                &lmNew,
                innovationLikelihood * this->particleSet_[i]->getData()->tracks_[m].getWeight(g) * Pd[m] * p / (1 - p),
                true);
          }

        } // gaussian in track m for loop

        num += rho * p * (1 - p) / pow(1 - p * Pd[m], 2);
        denom += rho * p / (1 - p * Pd[m]);
      } // existing track for loop
      this->particleSet_[i]->getData()->tracks_[nM + z].normalizeWeights();
      denom += this->pMeasurementModel_->clutterIntensity(this->measurements_[z], nZ);
      this->particleSet_[i]->getData()->tracks_[nM + z].setP(num / denom);

    } // measurements for loop

    if (nThreads_ == 1)
      timer_mapUpdate_kf_.stop();

    //----------  5. Identify unused measurements for adding birth Gaussians later ----------
    unused_measurements_[i].clear();
    int num_unused=0;
    for (int z = 0; z < nZ; z++) {

      if (this->particleSet_[i]->getData()->tracks_[nM + z].getP() == 0){
        this->particleSet_[i]->getData()->tracks_[nM + nZ-1-num_unused].copyTo(&(this->particleSet_[i]->getData()->tracks_[nM + z]));
        unused_measurements_[i].push_back(z);
        num_unused++;
        }
    }
    this->particleSet_[i]->getData()->tracks_.resize(nM+nZ-num_unused);

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::importanceWeighting(
      const uint idx) {

    const uint particleIdx = idx;
    TPose x;
    this->particleSet_[particleIdx]->getPose(x);

    // 1. select evaluation points from highest probability of existence tracks after update, that are within sensor FOV
    const unsigned int nM = this->particleSet_[particleIdx]->getData()->tracks_.size();
    int nEvalPoints = config.importanceWeightingEvalPointCount_ > nM ? nM : config.importanceWeightingEvalPointCount_;
    std::vector<unsigned int> evalPointIdx;
    std::vector<double> evalPointPd;
    evalPointIdx.reserve(nEvalPoints);
    evalPointPd.reserve(nEvalPoints);
    if (nEvalPoints == 0) {
      this->particleSet_[particleIdx]->setWeight(std::numeric_limits<double>::denorm_min());
      return;
    }
    this->particleSet_[particleIdx]->getData()->sortByPE(); // sort by probability of existence so that we can pick off the top nEvalPoints Gaussians
    for (int m = 0; m < nM; m++) {
      TLandmark* plm_temp;
      double w, w_prev;
      bool closeToSensingLim;
      if (this->particleSet_[particleIdx]->getData()->tracks_[m].getP() < config.importanceWeightingEvalPointTrackPE_)
        break;
      double Pd = 0;
      for (int i = 0; i < this->particleSet_[particleIdx]->getData()->tracks_[m].getGaussianCount(); i++) {
        this->particleSet_[particleIdx]->getData()->tracks_[m].getGaussian(i, plm_temp, w, w_prev);
        if (plm_temp != NULL) {
          Pd += this->pMeasurementModel_->probabilityOfDetection(x, *plm_temp, closeToSensingLim);

        }
      } // loop over gaussians
      Pd = Pd / this->particleSet_[particleIdx]->getData()->tracks_[m].getGaussianCount();

      if (Pd > 0) {
        evalPointIdx.push_back(m);
        evalPointPd.push_back(Pd);
      }
      if (nEvalPoints != -1 && evalPointIdx.size() >= nEvalPoints)
        break;
    } // loop over tracks
    nEvalPoints = evalPointIdx.size();

    double mapLikelihood = this->rfsMapLikelihood(particleIdx, evalPointIdx);
    double prevMapLikelihood = this->rfsPrevMapLikelihood(particleIdx, evalPointIdx);

    // 4. calculate measurement likelihood at eval points
    // note that rfsMeasurementLikelihood uses maps_[particleIdx] which is already sorted by weight
    double measurementLikelihood = rfsMeasurementLikelihood(particleIdx, evalPointIdx, evalPointPd);
    //printf("Particle %d measurement likelihood = %f\n", i, measurementLikelihood);

    // 5. calculate overall weight
    double overall_weight = measurementLikelihood * mapLikelihood / prevMapLikelihood;

    double prev_weight = this->particleSet_[particleIdx]->getWeight();
    this->particleSet_[particleIdx]->setWeight(overall_weight * prev_weight);
    //printf("Particle %d overall weight = %f\n\n", i, overall_weight);

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  double RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::rfsPrevMapLikelihood(
      const int particleIdx, std::vector<unsigned int>& evalPtIdx) {

    TPose x;
    this->particleSet_[particleIdx]->getPose(x);
    // find all the tracks in the field of view of the robot in the previous map estimate.
    std::vector<unsigned int> tracksInFoVIdx;
    tracksInFoVIdx.reserve(evalPtIdx.size());

    for (int m = 0; m < this->particleSet_[particleIdx]->getData()->tracks_.size(); m++) {
      TLandmark* plm_temp;
      double w, w_prev;
      bool closeToSensingLim;

      if (this->particleSet_[particleIdx]->getData()->tracks_[m].getP() > 0) {

        for (int i = 0; i < this->particleSet_[particleIdx]->getData()->tracks_[m].getGaussianCount(); i++) {
          this->particleSet_.at(particleIdx)->getData()->tracks_[m].getGaussian(i, plm_temp, w, w_prev);
          if (plm_temp != NULL && this->particleSet_[particleIdx]->getData()->tracks_[m].getPrevP() > 0) {
            double Pd = this->pMeasurementModel_->probabilityOfDetection(x, *plm_temp, closeToSensingLim);
            if (Pd > 0) {
              tracksInFoVIdx.push_back(m);
              break;
            }
          }
        } // loop over gaussians
      } // loop over tracks

    }

    // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight;

    const int i = particleIdx;

    const int nM = evalPtIdx.size();
    const int nMM = tracksInFoVIdx.size();

    TLandmark* evalPt;
    TLandmark evalPt_copy;
    double const threshold = config.importanceWeightingMeasurementLikelihoodMDThreshold_
        * config.importanceWeightingMeasurementLikelihoodMDThreshold_;
    double md2; // Mahalanobis distance squared

    // Create and fill in likelihood table (nMM x nM)
    double** L;
    CostMatrixGeneral likelihoodMatrix(L, nMM, nM);

    for (int m = 0; m < nM; m++) {

      this->particleSet_[i]->getData()->tracks_[evalPtIdx[m]].getMaxGaussian(evalPt); // get location of m
      evalPt_copy = *evalPt; // so that we don't change the actual data //
      evalPt_copy.setCov(MeasurementModel::TLandmark::Mat::Zero()); //

      for (int n = 0; n < nMM; n++) {

        // calculate the likelihood with existence probability
        L[n][m] = 0;
        for (unsigned int g = 0; g < this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getGaussianCount(); g++) {
          TLandmark* lm = this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getGaussian(g);
          double weight = this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getWeight(g);
          double like = lm->evalGaussianLikelihood(*evalPt);
          L[n][m] += like*weight;
        } // gaussian in track m for loop
        L[n][m] *= this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getPrevP(); // new line

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
      if (nRows < nCols){
    	  std::cerr << "one of the evalpoints does not exist!!\n";
    	  return 0;
      }

      if (isZeroPartition ) { // No tracks in this partition exist, nCols should be zero

    	  if (nCols != 0){
    		  std::cerr << "one of the evalpoints does not exist!!\n";
    		  std::cerr << "Error: Prev Map likelihood is zero \n";
    		  std::cerr << "ZeroP: " << isZeroPartition << "  nRows: " << nRows << "  nCols:  " << nCols << "\n";
    		  std::cerr << "nP:  " << nP << "  p: "<<p  << "  nMM: " << nMM << "  nM: " << nM <<"\n";
    		  return 0;
    	  }
    	  partition_likelihood = 1;
    	  for (int r = 0; r<nRows; r++){
    		  partition_likelihood *=1 - this->particleSet_[i]->getData()->tracks_[rowIdx[r]].getPrevP();
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
                Cp[r][c] = log(1 - this->particleSet_[i]->getData()->tracks_[rowIdx[r]].getPrevP());
              else
                Cp[r][c] = BIG_NEG_NUM;
            }
          }

          Murty murtyAlgo(Cp, nRows);
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
                permutation_log_likelihood += log(1 - this->particleSet_[i]->getData()->tracks_[rowIdx[a]].getP());
              }
            }
            for (int a = nRows; a < nRows + nCols; a++) { // outliers
              if (o[a] < nCols) {
                return 0;
              }
            }
            partition_likelihood += exp(permutation_log_likelihood);
            nPerm = pl.next(o);
          }

        } // End lexicographic ordering

      } // End non zero partition

      l *= partition_likelihood;

    } // End partitions

    return (l);
  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  double RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::rfsMapLikelihood(
      const int particleIdx, std::vector<unsigned int>& evalPtIdx) {

    TPose x;
    this->particleSet_[particleIdx]->getPose(x);
    // find all the tracks in the field of view of the robot in the updated map estimate.
    std::vector<unsigned int> tracksInFoVIdx;
    tracksInFoVIdx.reserve(evalPtIdx.size());

    for (int m = 0; m < this->particleSet_[particleIdx]->getData()->tracks_.size(); m++) {
      TLandmark* plm_temp;
      double w, w_prev;
      bool closeToSensingLim;

      if (this->particleSet_[particleIdx]->getData()->tracks_[m].getP() > 0) {

        for (int i = 0; i < this->particleSet_[particleIdx]->getData()->tracks_[m].getGaussianCount(); i++) {
          this->particleSet_[particleIdx]->getData()->tracks_[m].getGaussian(i, plm_temp, w, w_prev);
          if (plm_temp != NULL) {
            double Pd = this->pMeasurementModel_->probabilityOfDetection(x, *plm_temp, closeToSensingLim);
            if (Pd > 0) {
              tracksInFoVIdx.push_back(m);
              break;
            }
          }
        } // loop over gaussians
      } // loop over tracks

    }

    // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight;

    const int i = particleIdx;

    const int nM = evalPtIdx.size();
    const int nMM = tracksInFoVIdx.size();

    TLandmark* evalPt;
    TLandmark evalPt_copy;
    double const threshold = config.importanceWeightingMeasurementLikelihoodMDThreshold_
        * config.importanceWeightingMeasurementLikelihoodMDThreshold_;
    double md2; // Mahalanobis distance squared

    // Create and fill in likelihood table (nMM x nM)
    double** L;
    CostMatrixGeneral likelihoodMatrix(L, nMM, nM);

    for (int m = 0; m < nM; m++) {

      this->particleSet_[i]->getData()->tracks_[evalPtIdx[m]].getMaxGaussian(evalPt); // get location of m
      evalPt_copy = *evalPt; // so that we don't change the actual data //
      evalPt_copy.setCov(MeasurementModel::TLandmark::Mat::Zero()); //

      for (int n = 0; n < nMM; n++) {

        // calculate the likelihood with existence probability
        L[n][m] = 0;
        for (unsigned int g = 0; g < this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getGaussianCount(); g++) {
          TLandmark* lm = this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getGaussian(g);
          if (lm != NULL) {
            double weight =this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getWeight(g);
            double like = lm->evalGaussianLikelihood(*evalPt);
            L[n][m] += weight*like;
          }
        } // gaussian in track m for loop
        L[n][m] *= this->particleSet_[i]->getData()->tracks_[tracksInFoVIdx[n]].getP(); // new line

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

      isZeroPartition = !likelihoodMatrix.getPartition(p, Cp, nRows, nCols, rowIdx, colIdx, useMurtyAlgorithm);
      if (nRows < nCols){
    	  std::cerr << "one of the evalpoints does not exist!!\n";
    	  return 0;
      }

      if (isZeroPartition ) { // No tracks in this partition exist, nCols should be zero

    	  if (nCols != 0){
    		  std::cerr << "one of the evalpoints does not exist!!\n";
    		  std::cerr << "Error: Prev Map likelihood is zero \n";
    		  std::cerr << "ZeroP: " << isZeroPartition << "  nRows: " << nRows << "  nCols:  " << nCols << "\n";
    		  std::cerr << "nP:  " << nP << "  p: "<<p  << "  nMM: " << nMM << "  nM: " << nM <<"\n";
    		  return 0;
    	  }
    	  partition_likelihood = 1;
    	  for (int r = 0; r<nRows; r++){
    		  partition_likelihood *=1 - this->particleSet_[i]->getData()->tracks_[rowIdx[r]].getP();
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
            for (int c = nCols; c < nRows; c++) {

              Cp[r][c] = log(1 - this->particleSet_[i]->getData()->tracks_[rowIdx[r]].getP());

            }
          }

          Murty murtyAlgo(Cp, nRows );
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
                permutation_log_likelihood += log(1 - this->particleSet_[i]->getData()->tracks_[rowIdx[a]].getP());
              }
            }
            for (int a = nRows; a < nRows + nCols; a++) { // outliers
              if (o[a] < nCols) {
                return 0;
              }
            }
            partition_likelihood += exp(permutation_log_likelihood);
            nPerm = pl.next(o);
          }

        } // End lexicographic ordering

      } // End non zero partition

      l *= partition_likelihood;

    } // End partitions

    return (l);
  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  double RBCBMeMBerFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::rfsMeasurementLikelihood(
      const int particleIdx, std::vector<unsigned int> &evalPtIdx, std::vector<double> &evalPtPd) {

    // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight;

    const int i = particleIdx;
    TPose x;
    this->particleSet_[i]->getPose(x);
    const int nM = evalPtIdx.size();
    const int nZ = this->measurements_.size();

    TLandmark* evalPt;
    TLandmark evalPt_copy;
    TMeasurement expected_z;
    double const threshold = config.importanceWeightingMeasurementLikelihoodMDThreshold_
        * config.importanceWeightingMeasurementLikelihoodMDThreshold_;
    double md2; // Mahalanobis distance squared

    // Create and fill in likelihood table (nM x nZ)
    double** L;
    CostMatrixGeneral likelihoodMatrix(L, nM, nZ);

    for (int m = 0; m < nM; m++) {

      this->particleSet_[i]->getData()->tracks_.at(evalPtIdx[m]).getMaxGaussian(evalPt); // get location of m
      evalPt_copy = *evalPt; // so that we don't change the actual data //
      assert(!evalPt_copy.get().isZero());
      evalPt_copy.setCov(MeasurementModel::TLandmark::Mat::Zero()); //
      this->pMeasurementModel_->measure(x, evalPt_copy, expected_z); // get expected measurement for m
      double Pd = evalPtPd[m]; // get the prob of detection of m

      for (int n = 0; n < nZ; n++) {

        // calculate measurement likelihood with detection statistics
        L[m][n] = expected_z.evalGaussianLikelihood(this->measurements_[n], &md2) * Pd; // new line
        if (md2 > threshold) {
          L[m][n] = 0;
        }
      }
    }

    // Partition the Likelihood Table and turn into a log-likelihood table
    int nP = likelihoodMatrix.partition();
    double l = 1;
    double const BIG_NEG_NUM = -1000; // to represent log(0)
    double clutter[nZ];
    for (int n = 0; n < nZ; n++) {
      clutter[n] = this->pMeasurementModel_->clutterIntensity(this->measurements_[n], nZ);
    }

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
          partition_likelihood *= evalPtPd[rowIdx[r]];
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
                Cp[r][c] = log(1 - evalPtPd[rowIdx[r]]);
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
                permutation_log_likelihood += log(1 - evalPtPd[rowIdx[a]]);
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

    return (l / this->pMeasurementModel_->clutterIntensityIntegral(nZ));
  }

}  // rfs  name space

#endif
