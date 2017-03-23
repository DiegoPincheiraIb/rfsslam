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
#ifndef RBLMBFILTER_HPP
#define RBLMBFILTER_HPP

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
 *  \class RBLMBFilter
 *  \brief Rao-Blackwellized Labeled  Multi Bernoulli Filter class
 *  
 *  This class implements the Rao-Blackwellized Labeled  Multi Bernoulli Filter, as described in
 *  "The Labeled Multi-Bernoulli SLAM Filter" by Hendrik Deusch, Stephan Reuter, and Klaus Dietmayer.
 *  The constructor of this class will internally instantiate the
 *  process model for both the robot and landmarks, the measurement model, 
 *  and the Kalman filter. Users have access to these through pointers that 
 *  can be obtained by calling the appropraite get function.
 *
 *  \tparam RobotProcessModel A robot process model derived from ProcessModel
 *  \tparam LmkProcessModel A landmark process model derived from ProcessModel
 *  \tparam MeasurementModel A sensor model derived from MeasurementModel
 *  \tparam KalmanFilter A Kalman filter that uses LmkProcessModel and MeasurementModel
 *  \author  Felipe Inostroza
 */
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  class RBLMBFilter : public ParticleFilter<RobotProcessModel, MeasurementModel,
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
     * \brief Configurations for this RBLMBFilter
     */
    struct Config
    {

      /**  New birth Tracks are set this probability of existence */
      double birthTrackMaxPE_;
      double birthTrackMinPE_;
      double expectedBirths_;

      /**  New Gaussians are only created during map update if the innovation mahalanobis distance
       is less than this threshold */
      double newGaussianCreateInnovMDThreshold_;





      /** The threshold used to determine if a possible meaurement-landmark
       *  pairing is significant to worth considering
       */
      double MeasurementLikelihoodThreshold_;

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
    RBLMBFilter(int n);

    /** Destructor */
    ~RBLMBFilter();
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
    std::vector<std::vector<double> > birthPe_;

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
    double updateMap(const uint particleIdx);

    /**
     * Importance weighting does nothing, this is done together with the map update!
     * \param[in] idx particle index
     */
    void importanceWeighting(const uint idx){};

    /**
     * Random Finite Set measurement likelihood evaluation
     * \brief The current measurements in measurements_ are used to determine the
     * RFS measurement likelihood given a set of landmarks
     * \param[in] particleIdx particle for which the likelihood is calcuated
     * \param[in] evalPtIdx indices of evaluation points in maps_[particleIdx]
     * \param[in] evalPtPd probability of detection of evaluation point
     * \return measurement likelihood
     */
    double rfsMeasurementLikelihood( const int particleIdx, std::vector<unsigned int> &tracksinfov, std::vector<double> &Pd, double** L, CostMatrixGeneral &likelihoodMatrix);



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
  LmkProcessModel* RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getLmkProcessModel() {
    return lmkModelPtr_;
  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::predict(
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
  void RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::addBirthTracks() {

    for (int i = 0; i < this->nParticles_; i++) {

      while (unused_measurements_[i].size() > 0) {

        // get measurement
        int unused_idx = unused_measurements_[i].back();
        TMeasurement unused_z = this->measurements_[unused_idx];
        unused_measurements_[i].pop_back();
        double probE= birthPe_[i].back();
        birthPe_[i].pop_back();

        // use inverse measurement model to get landmark
        TPose robot_pose;
        TLandmark landmark_pos;
        robot_pose= * this->particleSet_[i]->getPose();
        this->pMeasurementModel_->inverseMeasure(robot_pose, unused_z, landmark_pos);

        // add birth landmark to Gaussian mixture (last param = true to allocate mem)
        this->particleSet_[i]->getData()->tracks_.resize(this->particleSet_[i]->getData()->tracks_.size() + 1);
        this->particleSet_[i]->getData()->tracks_[this->particleSet_[i]->getData()->tracks_.size() - 1].addGaussian(
            &landmark_pos, 1, true);
        this->particleSet_[i]->getData()->tracks_[this->particleSet_[i]->getData()->tracks_.size() - 1].setP(probE);

      }

    }

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::RBLMBFilter(int n) :
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
      birthPe_.push_back(std::vector<double>());
    }

    config.birthTrackMaxPE_ = 0.5;
    config.birthTrackMinPE_ = 0.01;

    config.gaussianMergingThreshold_ = 0.5;
    config.gaussianMergingCovarianceInflationFactor_ = 1.5;
    config.gaussianPruningThreshold_ = 0.2;
    config.trackPruningThreshold_ = 0.001;
    config.MeasurementLikelihoodThreshold_ = 1e-4;
    config.newGaussianCreateInnovMDThreshold_ = 0.2;
    config.minUpdatesBeforeResample_ = 1;

    nUpdatesSinceResample = 0;

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  KalmanFilter* RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getKalmanFilter() {
    return &(kfs_[0]);
  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::setParticlePose(int i,
                                                                                                             TPose &p) {

	*(this->particleSet_[i]) = p;

  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  bool RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getTrack(
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
  typename RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::TimingInfo* RBLMBFilter<
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
  int RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::getTrackNum(int i) {

    if (i >= 0 && i < this->nParticles_)
      return (this->particleSet_[i]->getData()->tracks_.size());
    else
      return -1;
  }

template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::~RBLMBFilter() {

    for (int i = 0; i < this->nParticles_; i++) {
      this->particleSet_[i]->deleteData();
    }
    delete lmkModelPtr_;

  }
template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
  void RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::update(
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
  double  RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::updateMap(
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
    birthPe_[i].clear();
    if (nM == 0) { // No existing landmark case -> flag all measurements as unused and go to next particles
      for (int z = 0; z < nZ; z++) {
        unused_measurements_[i].push_back(z);
        birthPe_[i].push_back(std::min(config.birthTrackMaxPE_,config.expectedBirths_/nZ));
      }
      return 1; // goto next particle
    }




    //----------  2. Kalman Filter map update ----------

    if (nThreads_ == 1)
      timer_mapUpdate_kf_.resume();
    const TPose *pose = this->particleSet_[i]->getPose();
    threshold_mahalanobisDistance2_mapUpdate_ = config.newGaussianCreateInnovMDThreshold_
        * config.newGaussianCreateInnovMDThreshold_;

    double innovationLikelihood;
    double innovationMahalanobisDist2;


    // determine the probability of detection for tracks in FoV


    std::vector<unsigned int> tracksInFovIdx;
    std::vector<double> tracksInFovPd,tracksInFovPe,tracksInFovUpdatedPe;
    std::vector<int> landmarkCloseToSensingLimit;

    // Calculate PD

    for (unsigned int m = 0; m < nM; m++) {

      bool isCloseToSensingLimit = false;
      double Pd = 0;
      for (unsigned int g = 0; g < this->particleSet_[i]->getData()->tracks_[m].getGaussianCount(); g++) {
        TLandmark* lm = this->particleSet_[i]->getData()->tracks_[m].getGaussian(g);

        bool isClose;
        Pd += this->particleSet_[i]->getData()->tracks_[m].getWeight(g)
            * this->pMeasurementModel_->probabilityOfDetection(*pose, *lm, isClose);
        isCloseToSensingLimit = isCloseToSensingLimit || isClose;
      } // gaussian in track m for loop
      if(Pd>0){
        tracksInFovIdx.push_back(m);
        tracksInFovPd.push_back(Pd);
        tracksInFovPe.push_back(this->particleSet_[i]->getData()->tracks_[m].getP());
        landmarkCloseToSensingLimit.push_back(isCloseToSensingLimit);
      }

      // store probability of existence on prevP. might be unnnecesary
      double prevP = this->particleSet_[i]->getData()->tracks_[m].getP();
      this->particleSet_[i]->getData()->tracks_[m].setPrevP(prevP);

    }
    tracksInFovUpdatedPe.resize(tracksInFovIdx.size(),0); // updated p starts at zero
    // Calculate Clutter
    double clutter[nZ];
    for (int n = 0; n < nZ; n++) {
      clutter[n] = this->pMeasurementModel_->clutterIntensity(this->measurements_[n], nZ);
    }

    // setup Likelihood matrix

    double** L;
    CostMatrixGeneral likelihoodMatrix (L, tracksInFovIdx.size (), nZ);

    GaussianMixture<TLandmark> **mixtureMatrix; // matrix of updated landmarks
    double **mixtureWeight; // weight of the updated gaussians in the updated track note: the "legacy" part is at the end of each array
    mixtureMatrix = new GaussianMixture<TLandmark>*[tracksInFovIdx.size ()];
    mixtureWeight = new double*[tracksInFovIdx.size ()];
    for (int m = 0; m < tracksInFovIdx.size (); m++) {
      mixtureMatrix[m] = new GaussianMixture<TLandmark>[nZ];
      mixtureWeight[m] = new double[nZ+1];
      for(int z=0;z<nZ+1;z++)
        mixtureWeight[m][z]=0;
    }




    // generate updated tracks

    for (unsigned int z = 0; z < nZ; z++) {


      for (unsigned int m =0; m < tracksInFovIdx.size();m++) {

        double p =tracksInFovPe[m];


        double gamma = 0;
        for (unsigned int g = 0; g < this->particleSet_[i]->getData()->tracks_[tracksInFovIdx[m]].getGaussianCount(); g++) {
          TLandmark* lm = this->particleSet_[i]->getData()->tracks_[tracksInFovIdx[m]].getGaussian(g);
          TLandmark lmNew;
          if (kfs_[threadnum].correct(*pose, this->measurements_[z], *lm, lmNew, &innovationLikelihood,
                                      &innovationMahalanobisDist2)) {
            double w=innovationLikelihood * this->particleSet_[i]->getData()->tracks_[tracksInFovIdx[m]].getWeight(g) * tracksInFovPd[m]/clutter[z];

            gamma += w;
            mixtureMatrix[m][z].addGaussian(
                            &lmNew,
                            w,
                            true);
          }

        } // gaussian in track m for loop

        gamma=gamma*p;
        if(gamma > config.MeasurementLikelihoodThreshold_){
          L[m][z]=gamma; // fill likelihood matrix
        }else{
          L[m][z]=0;
        }


      } // existing track for loop

    } // measurements for loop





    // Partition the Likelihood Table and turn into a log-likelihood table
    int nP = likelihoodMatrix.partition();
    double l = 1;
    double const BIG_NEG_NUM = -1000; // to represent log(0)


    // Go through each partition and determine the likelihood and add the corresponding gaussian mixtures to tracks
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

      if (isZeroPartition || nCols == 0 ) { // all landmarks in this partition are mis-detected. All measurements are outliers

        partition_likelihood = 1;
        for (int r = 0; r < nRows; r++) {
          partition_likelihood *= 1-tracksInFovPd[rowIdx[r]]*tracksInFovPe[rowIdx[r]];

        }

        for (int r = 0; r < nRows; r++) {


          this->particleSet_[i]->getData ()->tracks_[tracksInFovIdx[rowIdx[r]]].setP (
                        tracksInFovPe[rowIdx[r]] * (1 - tracksInFovPd[rowIdx[r]])
                        / (1 - tracksInFovPd[rowIdx[r]] * tracksInFovPe[rowIdx[r]]));
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
                Cp[r][c] = log(1 - tracksInFovPd[rowIdx[r]] * tracksInFovPe[rowIdx[r]]);
              else
                Cp[r][c] = BIG_NEG_NUM;
            }
          }

          // clutter
          for (int r = nRows; r < nRows + nCols; r++) {
            for (int c = 0; c < nCols; c++) {
              if (r - nRows == c)
                Cp[r][c] = 0;
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

            for(int r = 0; r<nRows; r++){

              if(a[r]<nCols){ // landmark r is assigned
                double labelLikelihood = exp(permutation_log_likelihood);
                mixtureWeight[rowIdx[r]][colIdx[a[r]]] +=labelLikelihood;
                tracksInFovUpdatedPe[rowIdx[r]] +=labelLikelihood;

              }else{  // landmark r unassigned
                double labelLikelihood = exp(permutation_log_likelihood)*tracksInFovPe[rowIdx[r]]*(1-tracksInFovPd[rowIdx[r]])/(1 - tracksInFovPd[rowIdx[r]] * tracksInFovPe[rowIdx[r]]);
                mixtureWeight[rowIdx[r]][nZ] +=labelLikelihood;
                tracksInFovUpdatedPe[rowIdx[r]] +=labelLikelihood;
              }

            }


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
              else { // mis-detection or non-existence
                permutation_log_likelihood += log(1 - tracksInFovPd[rowIdx[a]] * tracksInFovPe[rowIdx[a]]);
              }
            }

            partition_likelihood += exp(permutation_log_likelihood);

            for(int r = 0; r<nRows; r++){

              if(o[r]<nCols){ // landmark r is assigned
                double labelLikelihood = exp(permutation_log_likelihood);
                mixtureWeight[rowIdx[r]][colIdx[o[r]]] +=labelLikelihood;
                tracksInFovUpdatedPe[rowIdx[r]] +=labelLikelihood;

              }else{  // landmark r unassigned
                double labelLikelihood = exp(permutation_log_likelihood)*tracksInFovPe[rowIdx[r]]*(1-tracksInFovPd[rowIdx[r]])/(1 - tracksInFovPd[rowIdx[r]] * tracksInFovPe[rowIdx[r]]);
                mixtureWeight[rowIdx[r]][nZ] +=labelLikelihood;
                tracksInFovUpdatedPe[rowIdx[r]] +=labelLikelihood;
              }
            }


            nPerm = pl.next(o);
          }

        } // End lexicographic ordering

        // -----------   Generate updated tracks
        for(int r = 0; r<nRows; r++){
          tracksInFovUpdatedPe[rowIdx[r]]=tracksInFovUpdatedPe[rowIdx[r]]/partition_likelihood;
          this->particleSet_[i]->getData ()->tracks_[tracksInFovIdx[rowIdx[r]]].setP (tracksInFovUpdatedPe[rowIdx[r]]);

          // "legacy" part of the track
          this->particleSet_[i]->getData ()->tracks_[tracksInFovIdx[rowIdx[r]]].multiplyWeights(mixtureWeight[rowIdx[r]][nZ]);
          //updated parts
          for(int c = 0; c< nCols; c++){
            this->particleSet_[i]->getData ()->tracks_[tracksInFovIdx[rowIdx[r]]].addGM(mixtureMatrix[rowIdx[r]][colIdx[c]],mixtureWeight[rowIdx[r]][colIdx[c]]);
          }
          this->particleSet_[i]->getData ()->tracks_[tracksInFovIdx[rowIdx[r]]].normalizeWeights();
        }

      } // End non zero partition

      l *= partition_likelihood;

    } // End partitions



    if (nThreads_ == 1)
      timer_mapUpdate_kf_.stop();

    //----------  5. Identify unused measurements for adding birth Gaussians later, calculate the birth probability of existence based on the adaptive PE strategy(see LMB paper) ----------
    unused_measurements_[i].clear();
    birthPe_[i].clear();
    int num_unused=0;
    std::vector<double> associationP(nZ, 0);
    double noAssociationPsum=0;
    for (int z = 0; z < nZ; z++) {
      bool is_unused=true;

      for(int l =0; l < tracksInFovIdx.size(); l++ ){
        associationP[z]+=mixtureWeight[l][z];

      }
      noAssociationPsum+=1-associationP[z];
    }

    for(int z = 0; z < nZ; z++) {
      double p= std::min(config.birthTrackMaxPE_ , (1-associationP[z]));

      if (p>config.birthTrackMinPE_){

        unused_measurements_[i].push_back(z);
        birthPe_[i].push_back(p);
        num_unused++;
      }
    }

    // free memory

    for (int m = 0; m < tracksInFovIdx.size (); m++) {
      delete[] mixtureMatrix[m];
      delete[] mixtureWeight[m];
    }
    delete[] mixtureMatrix;
    delete[] mixtureWeight;

    // Update importance Weight

    double prev_weight = this->particleSet_[particleIdx]->getWeight();
    this->particleSet_[particleIdx]->setWeight(l * prev_weight);
    return l;


}



template<class RobotProcessModel, class LmkProcessModel, class MeasurementModel, class KalmanFilter>
double RBLMBFilter<RobotProcessModel, LmkProcessModel, MeasurementModel, KalmanFilter>::rfsMeasurementLikelihood(
    const int particleIdx, std::vector<unsigned int> &tracksinfov, std::vector<double> &Pd, double** L, CostMatrixGeneral &likelihoodMatrix) {

  // eval points are first nEvalPoints elements of maps_[i], which are already ordered by weight;

  const int i = particleIdx;

  const int nM = Pd.size();
  const int nZ = this->measurements_.size();







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
        partition_likelihood *= Pd[rowIdx[r]];
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
              Cp[r][c] = log(1 - Pd[rowIdx[r]]);
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
              permutation_log_likelihood += log(1 - Pd[rowIdx[a]]);
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
