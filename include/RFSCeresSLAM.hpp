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
#ifndef RFSCERESSLAM_HPP
#define RFSCERESSLAM_HPP

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
#include "ceres/ceres.h"
#include <unordered_map>

#ifdef _PERFTOOLS_CPU
#include <gperftools/profiler.h>
#endif
#ifdef _PERFTOOLS_HEAP
#include <gperftools/heap-profiler.h>
#endif

namespace rfs {

  /**
   *  \class RFSCeresSLAM
   *  \brief Random Finite Set  optimization using ceres solver for  feature based SLAM
   *
   *
   *
   *  \tparam RobotProcessModel A robot process model derived from ProcessModel
   *  \tparam MeasurementModel A sensor model derived from MeasurementModel
   *  \author  Felipe Inostroza
   */
  template<class RobotProcessModel, class MeasurementModel>
    class RFSCeresSLAM  : public ceres::FirstOrderFunction{

    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      typedef typename RobotProcessModel::TState TPose;
      typedef typename RobotProcessModel::TInput TInput;
      typedef typename MeasurementModel::TLandmark TLandmark;
      typedef typename MeasurementModel::TMeasurement TMeasurement;
      typedef RFSPSOParticle<RobotProcessModel, MeasurementModel> TParticle;
      typedef std::vector<TParticle> TParticleSet;
      static const int PoseDim =  TPose::Vec::RowsAtCompileTime;
      static const int LandmarkDim =  TLandmark::Vec::RowsAtCompileTime;
      /**
       * \brief Configurations for this RFSBatchPSO optimizer
       */
      struct Config {



        double ospa_c_; /**< ospa-like c value for calculating set differences and speed*/

        /** The threshold used to determine if a possible meaurement-landmark
         *  pairing is significant to worth considering
         */
        double MeasurementLikelihoodThreshold_;

        double mapFromMeasurementProb_; /**< probability that each measurement will initialize a landmark on map initialization*/


      } config;

      /**
       * Constructor
       */
      RFSCeresSLAM ();

      /** Destructor */
      ~RFSCeresSLAM ();

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
      initTrajectory ();

      /**
       * Generate random maps based on the already initialized trajectories
       */
      void
      initMap ();

      /**
       * return an initialization of the optimizer
       */
      double *
      init ();




      /**
       * Calculates the measurement likelihood of particle particleIdx at time k
       * @param k the time for which to calculate the likelihood
       * @return the measurement likelihood
       */

      double
      rfsMeasurementLikelihood (const int k, Eigen::Matrix<double , PoseDim , 1> &pose_gradient, Eigen::VectorXd &landmark_gradient) ;

      /**
       * Calculates the measurement likelihood of particle particleIdx  including all available times
       * @return the measurement likelihood
       */

      double
      rfsMeasurementLikelihood (Eigen::VectorXd &trajectory_gradient, Eigen::VectorXd &landmark_gradient) ;



      virtual bool Evaluate(const double* parameters,
                            double* cost,
                            double* gradient) const ;

      virtual int NumParameters() const ;



      MeasurementModel *mModelPtr_;
      RobotProcessModel *robotProcessModelPtr_;

      mutable std::vector<TPose> trajectory_;
      mutable std::vector<TLandmark> landmarks_;
    private:

      int nThreads_; /**< Number of threads  */




      std::vector<TInput> inputs_; /**< vector containing all odometry inputs */
      std::vector<std::vector<TMeasurement> > Z_; /**< vector containing all feature measurements */
      std::vector<TimeStamp> time_;


    };

  //////////////////////////////// Implementation ////////////////////////


  template<class RobotProcessModel, class MeasurementModel>
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::RFSCeresSLAM () {
      nThreads_ = 1;

#ifdef _OPENMP
      nThreads_ = omp_get_max_threads();
#endif
      mModelPtr_ = new MeasurementModel();
      robotProcessModelPtr_ = new RobotProcessModel();

    }

  template<class RobotProcessModel, class MeasurementModel>
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::~RFSCeresSLAM () {

      delete mModelPtr_;
      delete robotProcessModelPtr_;
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::addInput (TInput u) {

      inputs_.push_back(u);
      time_.push_back(u.getTime());
      std::sort(inputs_.begin(), inputs_.begin());
      std::sort(time_.begin(), time_.begin());
      Z_.resize(inputs_.size() + 1);

    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::setInputs (std::vector<TInput> U) {

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
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::addMeasurement (TMeasurement z) {

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
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::addMeasurement (std::vector<TMeasurement> Z) {

      for (int i = 0; i < Z.size(); i++) {
        this->addMeasurement(Z[i]);
      }
    }

  template<class RobotProcessModel, class MeasurementModel>
    double *
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::init () {

      initTrajectory();
      initMap();
      double * parameters = new double[NumParameters()];
      for (int k = 0; k < trajectory_.size(); k++) {
        for (int d = 0; d < PoseDim; d++) {
          parameters[k * PoseDim + d] = trajectory_[k][d];
        }
      }
      int trajectoryDim = trajectory_.size() * PoseDim;
      for (int m = 0; m < landmarks_.size(); m++) {
        for (int d = 0; d < LandmarkDim; d++) {
          parameters[trajectoryDim + m * LandmarkDim + d] = landmarks_[m][d];
        }
      }
      return parameters;
    }

  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::initTrajectory () {
      trajectory_.resize(inputs_.size() + 1);
      for (int k = 0; k < inputs_.size(); k++) {
        TimeStamp dT = inputs_[k].getTime() - trajectory_[k].getTime();
        robotProcessModelPtr_->sample(trajectory_[k + 1], trajectory_[k], inputs_[k], dT, false, true);
      }

    }
  template<class RobotProcessModel, class MeasurementModel>
    void
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::initMap () {



        for (int k = 0; k < trajectory_.size(); k++) {
          for (int nz = 0; nz < Z_[k].size(); nz++) {

            if (drand48() < config.mapFromMeasurementProb_) {

              TLandmark lm;
              this->mModelPtr_->inverseMeasure(trajectory_[k], Z_[k][nz], lm);
              landmarks_.push_back(lm);
            }
          }
        }




    }


  template<class RobotProcessModel, class MeasurementModel>
        int
        RFSCeresSLAM<RobotProcessModel, MeasurementModel>::NumParameters() const {
    return trajectory_.size()*PoseDim + landmarks_.size()*LandmarkDim;
  }

  template<class RobotProcessModel, class MeasurementModel>
      bool
      RFSCeresSLAM<RobotProcessModel, MeasurementModel>::Evaluate(const double* parameters,
                                                                  double* cost,
                                                                  double* gradient) const {



    RFSCeresSLAM<RobotProcessModel, MeasurementModel>  copy = *this;
    copy.config = this->config;
    copy.mModelPtr_ = new MeasurementModel();
    typename MeasurementModel::TMeasurement::Mat R;
    this->mModelPtr_->getNoise(R);
    copy.mModelPtr_->setNoise(R);
    copy.mModelPtr_->config = this->mModelPtr_->config;
    copy.robotProcessModelPtr_ = new RobotProcessModel();
    *copy.robotProcessModelPtr_ = *this->robotProcessModelPtr_;

    copy.trajectory_.resize(trajectory_.size());
    for (int k = 0; k < copy.trajectory_.size() ; k++ ){
      for (int d = 0; d < copy.PoseDim ; d++){

        copy.trajectory_[k][d] = parameters[k*PoseDim + d];
      }
    }
    int trajectoryDim = copy.trajectory_.size()*copy.PoseDim;
    copy.landmarks_.resize(landmarks_.size());
    for(int m = 0; m <  copy.landmarks_.size() ;  m++){
      for(int d = 0; d < copy.LandmarkDim ; d++){

        copy.landmarks_[m][d] = parameters[trajectoryDim + m*LandmarkDim +d];
      }
    }
    Eigen::VectorXd trajectory_gradient , landmark_gradient;
    (*cost) = -copy.rfsMeasurementLikelihood(trajectory_gradient , landmark_gradient);

    if(gradient !=  NULL){
    for( int i = 0; i < trajectory_gradient.size() ;  i++){

        gradient[i] = -trajectory_gradient[i];


    }

    for(int i = 0; i <  landmark_gradient.size() ;  i++){
        gradient[trajectoryDim + i] = -landmark_gradient[i];


    }

    }
    return true;
  }



  template<class RobotProcessModel, class MeasurementModel>
    double
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::rfsMeasurementLikelihood (Eigen::VectorXd &trajectory_gradient, Eigen::VectorXd &landmark_gradient) {
    trajectory_gradient.resize(PoseDim * trajectory_.size(),1);
    landmark_gradient.resize(LandmarkDim * landmarks_.size(),1);
    trajectory_gradient.setZero();
    landmark_gradient.setZero();
    typename TPose::Vec pose_gradient ;
      double l =  log(rfsMeasurementLikelihood( 0, pose_gradient , landmark_gradient));
      trajectory_gradient.segment<PoseDim>(0) = pose_gradient;
      TimeStamp dT;
      for (int k = 1; k < trajectory_.size(); k++) {
        Eigen::VectorXd lmgrad;
        pose_gradient.setZero();
        l += log(rfsMeasurementLikelihood(k, pose_gradient , lmgrad));
        trajectory_gradient.segment<PoseDim>(PoseDim * k ) = pose_gradient;
        dT = time_[k] - time_[k - 1];
        typename TPose::Vec pose_process_gradient;
        l += log(robotProcessModelPtr_->likelihood(trajectory_[k], trajectory_[k-1], inputs_[k - 1],  dT, &pose_process_gradient));

        landmark_gradient += lmgrad;
      }

      //std::cout << "likelihood   " << l << "\n";
      return l;
    }

  template<class RobotProcessModel, class MeasurementModel>
    double
    RFSCeresSLAM<RobotProcessModel, MeasurementModel>::rfsMeasurementLikelihood ( const int k, Eigen::Matrix<double , PoseDim , 1> &pose_gradient ,  Eigen::VectorXd &landmark_gradient) {



      const TPose &pose = trajectory_[k];
      const int nZ = this->Z_[k].size();
      const unsigned int mapSize = landmarks_.size();


      pose_gradient.setZero();
      landmark_gradient.resize(mapSize * LandmarkDim , 1 );
      landmark_gradient.setZero();
      //landmark_gradient.setZero();

      // Find map points within field of view and their probability of detection
      std::vector<unsigned int> lmInFovIdx;
      lmInFovIdx.reserve(mapSize);
      std::vector<double> lmInFovPd;
      lmInFovPd.reserve(mapSize);
      std::vector<int> landmarkCloseToSensingLimit;
      landmarkCloseToSensingLimit.reserve(mapSize);

      for (unsigned int m = 0; m < mapSize; m++) {

        bool isCloseToSensingLimit = false;

        const TLandmark &lm = this->landmarks_[m];

        bool isClose;
        double Pd = this->mModelPtr_->probabilityOfDetection(pose, lm, isCloseToSensingLimit);

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

      Eigen::Matrix<double,  MeasurementModel::TMeasurement::Vec::RowsAtCompileTime ,  MeasurementModel::TPose::Vec::RowsAtCompileTime >   jacobian_wrt_pose;
      Eigen::Matrix<double, MeasurementModel::TMeasurement::Vec::RowsAtCompileTime ,  MeasurementModel::TLandmark::Vec::RowsAtCompileTime >  jacobian_wrt_lmk;

      std::unordered_map< int , typename MeasurementModel::TLandmark::Vec> landmark_gradients;
      std::unordered_map< int , typename MeasurementModel::TPose::Vec> pose_gradients;




      for (int m = 0; m < nM; m++) {

        evalPt = &this->landmarks_[lmInFovIdx[m]]; // get location of m
        evalPt_copy = *evalPt; // so that we don't change the actual data //
        evalPt_copy.setCov(MeasurementModel::TLandmark::Mat::Zero()); //

        this->mModelPtr_->measure(pose, evalPt_copy, expected_z, &(jacobian_wrt_lmk) , &(jacobian_wrt_pose) ); // get expected measurement for m
        double Pd = lmInFovPd[m]; // get the prob of detection of m

        for (int n = 0; n < nZ; n++) {
          typename MeasurementModel::TMeasurement::Vec n_error;
          // calculate measurement likelihood with detection statistics
          L[m][n] = expected_z.evalGaussianLikelihood(this->Z_[k][n], n_error, &md2) * Pd; // new line
          if (L[m][n] < config.MeasurementLikelihoodThreshold_) {
            L[m][n] = 0;
          }else{

            typename MeasurementModel::TLandmark::Vec lm_grad ;
            typename MeasurementModel::TPose::Vec pose_grad ;
            lm_grad = jacobian_wrt_lmk.transpose() * n_error;
            pose_grad = jacobian_wrt_pose.transpose() * n_error;

            landmark_gradients.insert(std::make_pair(m*nZ+n ,lm_grad ));
            pose_gradients.insert(std::make_pair(m*nZ+n,pose_grad));

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
        Eigen::Matrix<double , PoseDim , 1> partition_pose_gradient;
        partition_pose_gradient.setZero();

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
            double permutation_likelihood = 0;
            murtyAlgo.setRealAssignmentBlock(nRows, nCols);
            for (int k = 0; k < 200; k++) {
              int rank = murtyAlgo.findNextBest(a, permutation_log_likelihood);
              if (rank == -1 || permutation_log_likelihood < BIG_NEG_NUM)
                break;
              permutation_likelihood = exp(permutation_log_likelihood);
              partition_likelihood += permutation_likelihood;

              // find the gradients
              if (permutation_likelihood>0){
              for (int r = 0; r < nRows; r++) {
                if (a[r] < nCols) {

                  typename MeasurementModel::TLandmark::Vec lm_grad = landmark_gradients.at(rowIdx[r]*nZ+ colIdx[a[r]]);
                  typename MeasurementModel::TPose::Vec pose_grad = pose_gradients.at(rowIdx[r]*nZ + colIdx[a[r]]);
                  landmark_gradient.segment<MeasurementModel::TLandmark::Vec::RowsAtCompileTime>(rowIdx[r] * MeasurementModel::TLandmark::Vec::RowsAtCompileTime) += lm_grad * permutation_likelihood;
                  partition_pose_gradient += pose_grad * permutation_likelihood;
                }
              }
              }
            }


          }
          else { // use lexicographic ordering

            partition_likelihood = 0;
            double permutation_log_likelihood = 0;
            double permutation_likelihood = 0;

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
              permutation_likelihood = exp(permutation_log_likelihood);
              partition_likelihood += permutation_likelihood;

              // find the gradients
              if (permutation_likelihood>0){
              for (int a = 0; a < nRows; a++) {
                if (o[a] < nCols) {
                  typename MeasurementModel::TLandmark::Vec lm_grad = landmark_gradients.at(rowIdx[a]*nZ+ colIdx[o[a]]);
                  typename MeasurementModel::TPose::Vec pose_grad  = pose_gradients.at(rowIdx[a]*nZ+ colIdx[o[a]]);
                  landmark_gradient.segment< MeasurementModel::TLandmark::Vec::RowsAtCompileTime >(rowIdx[a] * MeasurementModel::TLandmark::Vec::RowsAtCompileTime ) += lm_grad * permutation_likelihood;
                  partition_pose_gradient += pose_grad * permutation_likelihood;


                }
              }
              }



              nPerm = pl.next(o);
            }


          } // End lexicographic ordering



        } // End non zero partition

        // normalize landmark gradients


        pose_gradient += partition_pose_gradient / partition_likelihood;
        for (int r = 0; r < nRows; r++) {
          landmark_gradient.segment< MeasurementModel::TLandmark::Vec::RowsAtCompileTime >(rowIdx[r] * MeasurementModel::TLandmark::Vec::RowsAtCompileTime ) /= partition_likelihood;
        }


        l *= partition_likelihood;

      } // End partitions

      return (l / this->mModelPtr_->clutterIntensityIntegral(nZ));
    }

}

#endif
