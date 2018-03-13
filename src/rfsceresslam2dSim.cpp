/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
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

#define BOOST_NO_CXX11_SCOPED_ENUMS // required for boost/filesystem to work with C++11
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "ProcessModel_Odometry2D.hpp"
#include "RFSCeresSLAM.hpp"
#include "MeasurementModel_RngBrg.hpp"
#include <stdio.h>
#include <string>
#include <sys/ioctl.h>

#include "ceres/ceres.h"
#include "ceres/dynamic_autodiff_cost_function.h"

#ifdef _PERFTOOLS_CPU
#include <gperftools/profiler.h>
#endif
#ifdef _PERFTOOLS_HEAP
#include <gperftools/heap-profiler.h>
#endif

using namespace rfs;

struct ObsevationError1d {
  ObsevationError1d (double obs, double std_dev) :
      obs(obs), std_dev(std_dev) {
  }

  template<typename T>
    bool
    operator() (const T* const pose, const T* const landmark, T* residuals) const {

      residuals[0] = (landmark[0] - pose[0] - T(obs)) / std_dev;
      return true;
    }

  // Factory to hide the construction of the CostFunction object from
  // the client code.
  static ceres::CostFunction*
  Create (const double obs, const double std_dev) {
    return (new ceres::AutoDiffCostFunction<ObsevationError1d, 1, 1, 1>(new ObsevationError1d(obs, std_dev)));
  }

private:
  double obs, std_dev;
};

/**
 * \class Simulator_RFSPSOSLAM_2d
 * \brief A 2d SLAM Simulator using PSO
 * \author Felipe Inostroza
 */
class Simulator_RFSCERESSLAM_2d {

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  ;

  Simulator_RFSCERESSLAM_2d () {
    ceresslam_ = new RFSCeresSLAM<MotionModel_Odometry2d, MeasurementModel_RngBrg>();
  }

  ~Simulator_RFSCERESSLAM_2d () {

    if (ceresslam_ != NULL) {
      delete ceresslam_;
    }

  }

  /** Read the simulator configuration file */
  bool
  readConfigFile (const char* fileName) {

    cfgFileName_ = fileName;

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(fileName, pt);

    logToFile_ = false;
    if (pt.get("config.logging.logToFile", 0) == 1)
      logToFile_ = true;
    logDirPrefix_ = pt.get<std::string>("config.logging.logDirPrefix", "./");

    kMax_ = pt.get<int>("config.timesteps");
    dT_ = pt.get<double>("config.sec_per_timestep");
    dTimeStamp_ = TimeStamp(dT_);

    nSegments_ = pt.get<int>("config.trajectory.nSegments");
    max_dx_ = pt.get<double>("config.trajectory.max_dx_per_sec");
    max_dy_ = pt.get<double>("config.trajectory.max_dy_per_sec");
    max_dz_ = pt.get<double>("config.trajectory.max_dz_per_sec");
    min_dx_ = pt.get<double>("config.trajectory.min_dx_per_sec");
    vardx_ = pt.get<double>("config.trajectory.vardx");
    vardy_ = pt.get<double>("config.trajectory.vardy");
    vardz_ = pt.get<double>("config.trajectory.vardz");

    nLandmarks_ = pt.get<int>("config.landmarks.nLandmarks");

    rangeLimitMax_ = pt.get<double>("config.measurements.rangeLimitMax");
    rangeLimitMin_ = pt.get<double>("config.measurements.rangeLimitMin");
    rangeLimitBuffer_ = pt.get<double>("config.measurements.rangeLimitBuffer");
    Pd_ = pt.get<double>("config.measurements.probDetection");
    c_ = pt.get<double>("config.measurements.clutterIntensity");
    varzr_ = pt.get<double>("config.measurements.varzr");
    varzb_ = pt.get<double>("config.measurements.varzb");

    maxiter_ = pt.get<int>("config.optimizer.iterations");
    ospa_c_ = pt.get<double>("config.optimizer.ospaC");
    if (pt.get("config.optimizer.useDataAssociation", 0) == 1)
      useDataAssociation_ = true;
    else
      useDataAssociation_ = false;

    initMapProb_ = pt.get<double>("config.optimizer.initMapProb");

    pNoiseInflation_ = pt.get<double>("config.optimizer.predict.processNoiseInflationFactor");

    zNoiseInflation_ = pt.get<double>("config.optimizer.update.measurementNoiseInflationFactor");

    MeasurementLikelihoodThreshold_ = pt.get<double>("config.optimizer.MLThreshold");

    return true;
  }

  /** Generate a trajectory in 2d space 
     *  \param[in] randSeed random seed for generating trajectory
     */
    void generateTrajectory(int randSeed = 0){

      srand48( randSeed);
      initializeGaussianGenerators();

      TimeStamp t;
      int seg = 0;
      MotionModel_Odometry2d::TState::Mat Q;
      Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
      MotionModel_Odometry2d motionModel(Q);
      MotionModel_Odometry2d::TInput input_k(t);
      MotionModel_Odometry2d::TState pose_k(t);
      MotionModel_Odometry2d::TState pose_km(t);
      groundtruth_displacement_.reserve( kMax_ );
      groundtruth_pose_.reserve( kMax_ );
      groundtruth_displacement_.push_back(input_k);
      groundtruth_pose_.push_back(pose_k);

      for( int k = 1; k < kMax_; k++ ){

        t += dTimeStamp_;

        if( k <= 0 ){
          double dx = 0;
          double dy = 0;
          double dz = 0;
          MotionModel_Odometry2d::TInput::Vec d;
          MotionModel_Odometry2d::TInput::Vec dCovDiag;
          d << dx, dy, dz;
          dCovDiag << 0, 0, 0;
          input_k = MotionModel_Odometry2d::TInput(d, dCovDiag.asDiagonal(), k);
        }else if( k >= kMax_ / nSegments_ * seg ){
          seg++;
          double dx = drand48() * max_dx_ * dT_;
          while( dx < min_dx_ * dT_ ){
            dx = drand48() * max_dx_ * dT_;
          }
          double dy = (drand48() * max_dy_ * 2 - max_dy_) * dT_;
          double dz = (drand48() * max_dz_ * 2 - max_dz_) * dT_;
          MotionModel_Odometry2d::TInput::Vec d;
          MotionModel_Odometry2d::TInput::Vec dCovDiag;
          d << dx, dy, dz;
          dCovDiag << Q(0,0), Q(1,1), Q(2,2);
          input_k = MotionModel_Odometry2d::TInput(d, dCovDiag.asDiagonal(), k);
        }

        groundtruth_displacement_.push_back(input_k);
        groundtruth_displacement_.back().setTime(t);

        MotionModel_Odometry2d::TState x_k;
        motionModel.step(x_k, groundtruth_pose_[k - 1], input_k, dTimeStamp_);
        groundtruth_pose_.push_back(x_k);
        groundtruth_pose_.back().setTime(t);

      }

    }

    /** Generate odometry measurements */
    void generateOdometry() {

      odometry_.reserve(kMax_);
      MotionModel_Odometry2d::TInput zero;
      MotionModel_Odometry2d::TInput::Vec u0;
      u0.setZero();
      zero.set(u0, 0);
      //odometry_.push_back(zero);

      MotionModel_Odometry2d::TState::Mat Q;
      Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
      MotionModel_Odometry2d motionModel(Q);
      deadReckoning_pose_.reserve(kMax_);
      deadReckoning_pose_.push_back(groundtruth_pose_[0]);

      TimeStamp t;

      for( int k = 1; k < kMax_; k++){

        t += dTimeStamp_;
        double dt = dTimeStamp_.getTimeAsDouble();

        MotionModel_Odometry2d::TInput in = groundtruth_displacement_[k];
        MotionModel_Odometry2d::TState::Mat Qk = Q * dt * dt;
        in.setCov(Qk);
        MotionModel_Odometry2d::TInput out;
        in.sample(out);

        odometry_.push_back( out );

        MotionModel_Odometry2d::TState p;
        motionModel.step(p, deadReckoning_pose_[k - 1], odometry_[k], dTimeStamp_);
        p.setTime(t);
        deadReckoning_pose_.push_back( p );
      }

    }

    /** Generate landmarks */
    void generateLandmarks() {

      MeasurementModel_RngBrg measurementModel(varzr_, varzb_);
      MeasurementModel_RngBrg::TPose pose;

      groundtruth_landmark_.reserve(nLandmarks_);

      int nLandmarksCreated = 0;
      for (int k = 1; k < kMax_; k++) {

        if (k >= kMax_ / nLandmarks_ * nLandmarksCreated) {

          MeasurementModel_RngBrg::TPose pose;
          MeasurementModel_RngBrg::TMeasurement measurementToCreateLandmark;
          MeasurementModel_RngBrg::TMeasurement::Vec z;
          double r = drand48() * rangeLimitMax_;
          double b = drand48() * 2 * PI;
          z << r, b;
          measurementToCreateLandmark.set(z);
          MeasurementModel_RngBrg::TLandmark lm;

          measurementModel.inverseMeasure(groundtruth_pose_[k], measurementToCreateLandmark, lm);

          groundtruth_landmark_.push_back(lm);

          nLandmarksCreated++;

        }

      }

    }

    /** Generate landmark measurements */
    void generateMeasurements() {

      MeasurementModel_RngBrg measurementModel(varzr_, varzb_);
      MeasurementModel_RngBrg::TMeasurement::Mat R;
      measurementModel.getNoise(R);
      measurementModel.config.rangeLimMax_ = rangeLimitMax_;
      measurementModel.config.rangeLimMin_ = rangeLimitMin_;
      measurementModel.config.probabilityOfDetection_ = Pd_;
      measurementModel.config.uniformClutterIntensity_ = c_;
      double meanClutter = measurementModel.clutterIntensityIntegral();

      double expNegMeanClutter = exp(-meanClutter);
      double poissonPmf[100];
      double poissonCmf[100];
      double mean_pow_i = 1;
      double i_factorial = 1;
      poissonPmf[0] = expNegMeanClutter;
      poissonCmf[0] = poissonPmf[0];
      for (int i = 1; i < 100; i++) {
        mean_pow_i *= meanClutter;
        i_factorial *= i;
        poissonPmf[i] = mean_pow_i / i_factorial * expNegMeanClutter;
        poissonCmf[i] = poissonCmf[i - 1] + poissonPmf[i];
      }

      lmkFirstObsTime_.resize(groundtruth_landmark_.size());
      for (int m = 0; m < lmkFirstObsTime_.size(); m++) {
        lmkFirstObsTime_[m] = -1;
      }

      TimeStamp t;

      for (int k = 1; k < kMax_; k++) {

        t += dTimeStamp_;

        groundtruth_pose_[k];

        // Real detections
        for (int m = 0; m < groundtruth_landmark_.size(); m++) {

          bool success;
          MeasurementModel_RngBrg::TMeasurement z_m_k;
          success = measurementModel.sample(groundtruth_pose_[k], groundtruth_landmark_[m], z_m_k);
          if (success) {

            if (z_m_k.get(0) <= rangeLimitMax_ && z_m_k.get(0) >= rangeLimitMin_ && drand48() <= Pd_) {
              z_m_k.setTime(t);
              // z_m_k.setCov(R);
              measurements_.push_back(z_m_k);
              groundtruthDataAssociation_.push_back(m);


            if (lmkFirstObsTime_[m] == -1) {
              lmkFirstObsTime_[m] = t.getTimeAsDouble();
            }
            }
          }

        }

        // False alarms
        double randomNum = drand48();
        int nClutterToGen = 0;
        while (randomNum > poissonCmf[nClutterToGen]) {
          nClutterToGen++;
        }
        for (int i = 0; i < nClutterToGen; i++) {

          double r = drand48() * rangeLimitMax_;
          while (r < rangeLimitMin_)
            r = drand48() * rangeLimitMax_;
          double b = drand48() * 2 * PI - PI;
          MeasurementModel_RngBrg::TMeasurement z_clutter;
          MeasurementModel_RngBrg::TMeasurement::Vec z;
          z << r, b;
          z_clutter.set(z, t);
          measurements_.push_back(z_clutter);
          groundtruthDataAssociation_.push_back(-2);
        }

      }

    }

    /** Data Logging */
    void exportSimData() {

      if (!logToFile_)
        return;

      boost::filesystem::path dir(logDirPrefix_);
      boost::filesystem::create_directories(dir);

      boost::filesystem::path cfgFilePathSrc(cfgFileName_);
      std::string cfgFileDst(logDirPrefix_);
      cfgFileDst += "simSettings.cfg";
      boost::filesystem::path cfgFilePathDst(cfgFileDst.data());
      boost::filesystem::copy_file(cfgFilePathSrc, cfgFilePathDst, boost::filesystem::copy_option::overwrite_if_exists);

      TimeStamp t;

      FILE* pGTPoseFile;
      std::string filenameGTPose(logDirPrefix_);
      filenameGTPose += "gtPose.dat";
      pGTPoseFile = fopen(filenameGTPose.data(), "w");
      MotionModel_Odometry2d::TState::Vec x;
      for (int i = 0; i < groundtruth_pose_.size(); i++) {
        groundtruth_pose_[i].get(x, t);
        fprintf(pGTPoseFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), x(0), x(1), x(2));
      }
      fclose(pGTPoseFile);

      FILE* pGTLandmarkFile;
      std::string filenameGTLandmark(logDirPrefix_);
      filenameGTLandmark += "gtLandmark.dat";
      pGTLandmarkFile = fopen(filenameGTLandmark.data(), "w");
      MeasurementModel_RngBrg::TLandmark::Vec m;
      for (int i = 0; i < groundtruth_landmark_.size(); i++) {
        groundtruth_landmark_[i].get(m);
        fprintf(pGTLandmarkFile, "%f   %f   %f\n", m(0), m(1), lmkFirstObsTime_[i]);
      }
      fclose(pGTLandmarkFile);

      FILE* pOdomFile;
      std::string filenameOdom(logDirPrefix_);
      filenameOdom += "odometry.dat";
      pOdomFile = fopen(filenameOdom.data(), "w");
      MotionModel_Odometry2d::TInput::Vec u;
      for (int i = 0; i < odometry_.size(); i++) {
        odometry_[i].get(u, t);
        fprintf(pOdomFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), u(0), u(1), u(2));
      }
      fclose(pOdomFile);

      FILE* pMeasurementFile;
      std::string filenameMeasurement(logDirPrefix_);
      filenameMeasurement += "measurement.dat";
      pMeasurementFile = fopen(filenameMeasurement.data(), "w");
      MeasurementModel_RngBrg::TMeasurement::Vec z;
      for (int i = 0; i < measurements_.size(); i++) {
        measurements_[i].get(z, t);
        fprintf(pMeasurementFile, "%f   %f   %f\n", t.getTimeAsDouble(), z(0), z(1));
      }
      fclose(pMeasurementFile);

      FILE* pDeadReckoningFile;
      std::string filenameDeadReckoning(logDirPrefix_);
      filenameDeadReckoning += "deadReckoning.dat";
      pDeadReckoningFile = fopen(filenameDeadReckoning.data(), "w");
      MotionModel_Odometry2d::TState::Vec odo;
      for (int i = 0; i < deadReckoning_pose_.size(); i++) {
        deadReckoning_pose_[i].get(odo, t);
        fprintf(pDeadReckoningFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), odo(0), odo(1), odo(2));
      }
      fclose(pDeadReckoningFile);

    }






  /** PSO optimizer Setup */
  void
  setup () {

    double dt = dTimeStamp_.getTimeAsDouble();

    // configure robot motion model (only need to set once since timesteps are constant)
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= (pNoiseInflation_ * dt * dt);
    ceresslam_->robotProcessModelPtr_->setNoise(Q);



    // configure measurement model
    MeasurementModel_RngBrg::TMeasurement::Mat R;
    R << varzr_, 0, 0, varzb_;
    R *= zNoiseInflation_;
    ceresslam_->mModelPtr_->setNoise(R);
    ceresslam_->mModelPtr_->config.probabilityOfDetection_ = Pd_;
    ceresslam_->mModelPtr_->config.uniformClutterIntensity_ = c_;
    ceresslam_->mModelPtr_->config.rangeLimMax_ = rangeLimitMax_;
    ceresslam_->mModelPtr_->config.rangeLimMin_ = rangeLimitMin_;
    ceresslam_->mModelPtr_->config.rangeLimBuffer_ = rangeLimitBuffer_;


    // configure the filter

    ceresslam_->config.ospa_c_ = ospa_c_;
    ceresslam_->config.mapFromMeasurementProb_ = initMapProb_;
    ceresslam_->config.MeasurementLikelihoodThreshold_ = MeasurementLikelihoodThreshold_;

    ceresslam_->setInputs(odometry_);
    if (useDataAssociation_)
      ceresslam_->addMeasurement(measurements_, groundtruthDataAssociation_);
    else
      ceresslam_->addMeasurement(measurements_);

  }

  /** Run the simulator */
  void
  run () {

    printf("Running simulation\n\n");

#ifdef _PERFTOOLS_CPU
    std::string perfCPU_file = logDirPrefix_ + "rbcbmemberslam2dSim_cpu.prof";
    ProfilerStart(perfCPU_file.data());
#endif
#ifdef _PERFTOOLS_HEAP
    std::string perfHEAP_file = logDirPrefix_ + "rbcbmemberslam2dSim_heap.prof";
    HeapProfilerStart(perfHEAP_file.data());
#endif
    //////// Initialization by sampling motion model //////////

    FILE* pParticlePoseFile;
    if (logToFile_) {
      std::string filenameParticlePoseFile(logDirPrefix_);
      filenameParticlePoseFile += "particlePose.dat";
      pParticlePoseFile = fopen(filenameParticlePoseFile.data(), "w");
    }
    FILE* pBestParticlePoseFile;
    if (logToFile_) {
      std::string filenameBestParticlePoseFile(logDirPrefix_);
      filenameBestParticlePoseFile += "bestParticlePose.dat";
      pBestParticlePoseFile = fopen(filenameBestParticlePoseFile.data(), "w");
    }
    FILE* pLandmarkEstFile;
    if (logToFile_) {
      std::string filenameLandmarkEstFile(logDirPrefix_);
      filenameLandmarkEstFile += "landmarkEst.dat";
      pLandmarkEstFile = fopen(filenameLandmarkEstFile.data(), "w");
    }
    FILE* pBestLandmarkEstFile;
    if (logToFile_) {
      std::string filenameBestLandmarkEstFile(logDirPrefix_);
      filenameBestLandmarkEstFile += "landmarkEst.dat";
      pBestLandmarkEstFile = fopen(filenameBestLandmarkEstFile.data(), "w");
    }
    MotionModel_Odometry2d::TState x_i;
    int zIdx = 0;


    // configure robot motion model

    MotionModel_Odometry2d::TState::Mat Q;
    double dt = dTimeStamp_.getTimeAsDouble();
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= (pNoiseInflation_ * dt * dt);
    ceresslam_->robotProcessModelPtr_->setNoise(Q);

    // calculate ground truth likelihood!!!

    ceresslam_->landmarks_.resize(groundtruth_landmark_.size());
    //ceresslam_->landmarks_ = groundtruth_landmark_;
    ceresslam_->trajectory_.resize( groundtruth_pose_.size());
    //ceresslam_->trajectory_ = groundtruth_pose_;
    int nparams = ceresslam_->NumParameters();
    double * gtparams = new double[nparams];

    for (int i = 1; i < groundtruth_pose_.size(); i++) {
      for(int d =0; d < ceresslam_->PoseDim ;d++){
        gtparams[(i  - 1) * ceresslam_->PoseDim + d] = groundtruth_pose_[i][d];
      }
    }
    for (int i = 0; i < groundtruth_landmark_.size(); i++) {
      for(int d = 0 ; d < ceresslam_->LandmarkDim ; d++){
      gtparams[(i)*ceresslam_->LandmarkDim + d + (groundtruth_pose_.size()-1)*ceresslam_->PoseDim ] = groundtruth_landmark_[i][d];
      }
    }
    double cost;
    ceresslam_->Evaluate(gtparams, &cost, NULL);
    std::cout << "gt likelihood:  " << cost << " nparams:  " << ceresslam_->NumParameters() <<"\n";

/*
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= ( dt * dt);
    ceresslam_->robotProcessModelPtr_->setNoise(Q);
*/

    ceresslam_->Evaluate(gtparams, &cost, NULL);
        std::cout << "gt likelihood:  " << cost << " nparams:  " << ceresslam_->NumParameters() <<"\n";

    double * initParams = ceresslam_->init();
    //ceresslam_->landmarks_.resize(groundtruth_landmark_.size());

/*
    for (int i = 0; i < (groundtruth_pose_.size()-1)*ceresslam_->PoseDim; i++){

      gtparams[i] =initParams[i];
    }*/
    //initParams= gtparams;
    std::cout << "initparams:   ";
    for (int i = 0; i < ceresslam_->NumParameters(); i++) {
      std::cout << initParams[i] << "   ";
    }
    std::cout << "\n";

    ceresslam_->Evaluate(gtparams, &cost, NULL);
    std::cout << "gt likelihood:  " << cost << " nparams:  " << ceresslam_->NumParameters() << "\n";

    ceresslam_->Evaluate(initParams, &cost, NULL);
    std::cout << "init likelihood:  " << cost << " nparams:  " << ceresslam_->NumParameters() << "\n";

    ceresslam_->Evaluate(gtparams, &cost, NULL);
    std::cout << "gt likelihood:  " << cost << " nparams:  " << ceresslam_->NumParameters() << "\n";
    /*
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= ( pNoiseInflation_ * dt * dt);
    ceresslam_->robotProcessModelPtr_->setNoise(Q);
*/
    int iteration = 0;

    if (logToFile_) {
      /////////// log particle trajectories  //////////////
      ///

      fprintf(pParticlePoseFile, "%d   ", iteration);

      fprintf(pParticlePoseFile, "%f   ", 0.0);
      fprintf(pParticlePoseFile, "%f   ", 0.0);
      fprintf(pParticlePoseFile, "%f   ", 0.0);
      fprintf(pParticlePoseFile, "%f   ", 0.0);
      for (int k = 0; k < (ceresslam_->trajectory_.size() - 1)*ceresslam_->PoseDim; k++) {

        fprintf(pParticlePoseFile, "%f   ", initParams[k]);
      }
      fprintf(pParticlePoseFile, "\n");

      ///////////// log landmarks ///////////////////

      fprintf(pLandmarkEstFile, "%d    ", iteration);
      for (int l = 0; l < ceresslam_->landmarks_.size()*ceresslam_->LandmarkDim; l++) {
        MeasurementModel_RngBrg::TLandmark landmark;

        fprintf(pLandmarkEstFile, "%f   ", initParams[(ceresslam_->trajectory_.size() - 1)*ceresslam_->PoseDim + l]);
      }
      fprintf(pLandmarkEstFile, "\n");

    }



    // run the optimization process
    ceres::GradientProblem *problem_general_unconstrained;

    if (useDataAssociation_) {
      ceres::Problem problem;
      double * params = new double[ceresslam_->NumParameters() + 1];
      params[0] = 0;
      for (int i = 0; i < ceresslam_->NumParameters(); i++) {
        params[i + 1] = initParams[i];
      }
      for (int i = 0; i < odometry_.size(); ++i) {
        // Each Residual block takes a point and a camera as input and outputs a 2
        // dimensional residual. Internally, the cost function stores the observed
        // image location and compares the reprojection against the observation.

        ceres::CostFunction* cost_function = ObsevationError1d::Create(odometry_[i][0], sqrt(vardx_));
        problem.AddResidualBlock(cost_function,
        NULL /* squared loss */,
                                 params + i - 1, params + i);
      }
      for (int k = 0; k < ceresslam_->Z_.size(); k++) {
        for (int j = 0; j < ceresslam_->Z_[k].size(); j++) {
          if (ceresslam_->DA_[k][j] < 0)
            continue;

          ceres::CostFunction* cost_function = ObsevationError1d::Create(ceresslam_->Z_[k][j][0], sqrt(varzr_));
          problem.AddResidualBlock(cost_function,
          NULL /* squared loss */,
                                   params + k, params + ceresslam_->trajectory_.size() + ceresslam_->DA_[k][j]);
        }
      }
      problem.SetParameterBlockConstant(params);
      ceres::Solver::Options options;
      options.linear_solver_type = ceres::DENSE_SCHUR;
      options.minimizer_progress_to_stdout = true;
      ceres::Solver::Summary summary;
      ceres::Solve(options, &problem, &summary);
      std::cout << summary.FullReport() << "\n";
      for (int i = 0; i < ceresslam_->NumParameters(); i++) {
        initParams[i] = params[i + 1];
      }
    }
    else {

      problem_general_unconstrained = new ceres::GradientProblem(ceresslam_);

      ceres::GradientProblemSolver::Options options;
      options.minimizer_progress_to_stdout = true;
      //options.line_search_direction_type = ceres::NONLINEAR_CONJUGATE_GRADIENT;
      options.max_num_iterations = 10000;
      options.function_tolerance = 1e-10;

      ceres::GradientProblemSolver::Summary summary;

      ceres::Solve(options, *problem_general_unconstrained, initParams, &summary);

      std::cout << summary.FullReport() << "\n";



    }
    iteration = 1;

    if (logToFile_) {
      /////////// log particle trajectories  //////////////
      ///

      fprintf(pParticlePoseFile, "%d   ", iteration);

      fprintf(pParticlePoseFile, "%f   ", 0.0);
      fprintf(pParticlePoseFile, "%f   ", 0.0);
      fprintf(pParticlePoseFile, "%f   ", 0.0);
      fprintf(pParticlePoseFile, "%f   ", 0.0);
      for (int k = 0; k < (ceresslam_->trajectory_.size() - 1)*ceresslam_->PoseDim; k++) {

        fprintf(pParticlePoseFile, "%f   ", initParams[k]);
      }
      fprintf(pParticlePoseFile, "\n");

      ///////////// log landmarks ///////////////////

      fprintf(pLandmarkEstFile, "%d    ", iteration);
      for (int l = 0; l < ceresslam_->landmarks_.size()*ceresslam_->LandmarkDim; l++) {
        MeasurementModel_RngBrg::TLandmark landmark;

        fprintf(pLandmarkEstFile, "%f   ", initParams[(ceresslam_->trajectory_.size() - 1)*ceresslam_->PoseDim + l]);
      }
      fprintf(pLandmarkEstFile, "\n");

    }

    if (logToFile_) {
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }

delete problem_general_unconstrained; //< deletes ceresslam_ object
ceresslam_ = NULL;
  }

private:

  const char* cfgFileName_;

  int kMax_; /**< number of timesteps */
  double dT_; /**< duration of timestep in seconds */
  TimeStamp dTimeStamp_; /**< duration of timestep in timestamp */

  // Trajectory
  int nSegments_;
  double max_dx_;
  double max_dy_;
  double max_dz_;
  double min_dx_;
  double vardx_;
  double vardy_;
  double vardz_;
  std::vector<MotionModel_Odometry2d::TInput> groundtruth_displacement_;
  std::vector<MotionModel_Odometry2d::TState> groundtruth_pose_;
  std::vector<MotionModel_Odometry2d::TInput> odometry_;
  std::vector<MotionModel_Odometry2d::TState> deadReckoning_pose_;

  // Landmarks 
  int nLandmarks_;
  std::vector<MeasurementModel_RngBrg::TLandmark> groundtruth_landmark_;
  double varlmx_;
  double varlmy_;

  std::vector<double> lmkFirstObsTime_;

  // Range-Bearing Measurements
  double rangeLimitMax_;
  double rangeLimitMin_;
  double rangeLimitBuffer_;
  double Pd_;
  double c_;
  double varzr_;
  double varzb_;
  std::vector<MeasurementModel_RngBrg::TMeasurement> measurements_;
  std::vector<int> groundtruthDataAssociation_;

  // Filters
  RFSCeresSLAM<MotionModel_Odometry2d, MeasurementModel_RngBrg> *ceresslam_;
  bool useDataAssociation_;

  double ospa_c_;
  double initMapProb_;
  int maxiter_;

  double pNoiseInflation_;
  double zNoiseInflation_;

  double MeasurementLikelihoodThreshold_;

  bool logToFile_;

public:
  std::string logDirPrefix_;
};

int
main (int argc, char* argv[]) {

  Simulator_RFSCERESSLAM_2d sim;
  int seed = time(NULL);
  srand(seed);
  int trajNum = rand();
  std::string cfgFileName;
  boost::program_options::options_description desc("Options");
  desc.add_options()("help,h", "produce this help message")(
      "cfg,c", boost::program_options::value<std::string>(&cfgFileName)->default_value("cfg/rfsceresslam2dSim.xml"), "configuration xml file")(
      "trajectory,t", boost::program_options::value<int>(&trajNum), "trajectory number (default: a random integer)")(
      "seed,s", boost::program_options::value<int>(&seed), "random seed for running the simulation (default: based on current system time)");
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("cfg")) {
    cfgFileName = vm["cfg"].as<std::string>();
  }
  std::cout << "Configuration file: " << cfgFileName << std::endl;
  if (!sim.readConfigFile(cfgFileName.data())) {
    return -1;
  }

  if (vm.count("trajectory")) {
    trajNum = vm["trajectory"].as<int>();
  }
  std::cout << "Trajectory: " << trajNum << std::endl;
  sim.generateTrajectory(trajNum);

  sim.generateLandmarks();
  sim.generateOdometry();
  sim.generateMeasurements();

  sim.exportSimData();

  sim.setup();
  if (vm.count("seed")) {
    seed = vm["seed"].as<int>();
    std::cout << "Simulation random seed manually set to: " << seed << std::endl;
  }

  srand48(seed);
  initializeGaussianGenerators();

  // boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Simulation run time: %ws\n");

  sim.run();

  // std::cout << "mem use: " << MemProfile::getCurrentRSS() << "(" << MemProfile::getPeakRSS() << ")\n";
  //delete timer;

  return 0;

}
