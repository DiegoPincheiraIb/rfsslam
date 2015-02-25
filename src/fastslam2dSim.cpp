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

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "FastSLAM.hpp"
#include "KalmanFilter_RngBrg.hpp"
#include "MeasurementModel_RngBrg.hpp"
#include "ProcessModel_Odometry2D.hpp"
#include <stdio.h>
#include <string>

using namespace rfs;

/**
 * \class Simulator_FastSLAM_2d
 * \brief A 2d FastSLAM Simulator
 * \author Keith Leung
 */
class Simulator_FastSLAM_2d{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Simulator_FastSLAM_2d(){
    pFilter_ = NULL;
  }
  
  ~Simulator_FastSLAM_2d(){
    
    if(pFilter_ != NULL){
      delete pFilter_;
    }

  }

  /** Read the simulator configuration file */
  bool readConfigFile(const char* fileName){

     cfgFileName_ = fileName;

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(fileName, pt);

    logToFile_ = false;
    if( pt.get("config.logging.logToFile", 0) == 1 )
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
    varlmx_ = pt.get<double>("config.landmarks.varlmx");
    varlmy_ = pt.get<double>("config.landmarks.varlmy");

    rangeLimitMax_ = pt.get<double>("config.measurements.rangeLimitMax");
    rangeLimitMin_ = pt.get<double>("config.measurements.rangeLimitMin");
    rangeLimitBuffer_ = pt.get<double>("config.measurements.rangeLimitBuffer");
    Pd_ = pt.get<double>("config.measurements.probDetection");
    c_ = pt.get<double>("config.measurements.clutterIntensity");
    varzr_ = pt.get<double>("config.measurements.varzr");
    varzb_ = pt.get<double>("config.measurements.varzb");
   
    nParticles_ = pt.get("config.filter.nParticles", 200);

    pNoiseInflation_ = pt.get("config.filter.predict.processNoiseInflationFactor", 1.0);

    zNoiseInflation_ = pt.get("config.filter.update.measurementNoiseInflationFactor", 1.0);
    maxNDataAssocHypotheses_ = pt.get<int>("config.filter.update.maxNDataAssocHypotheses", 1);
    maxDataAssocLogLikelihoodDiff_ = pt.get("config.filter.update.maxDataAssocLogLikelihoodDiff", 3.0);

    innovationRangeThreshold_ = pt.get<double>("config.filter.update.KalmanFilter.innovationThreshold.range");
    innovationBearingThreshold_ = pt.get<double>("config.filter.update.KalmanFilter.innovationThreshold.bearing");

    minLogMeasurementLikelihood_ = pt.get("config.filter.weighting.minLogMeasurementLikelihood",-10.0);

    effNParticleThreshold_ = pt.get("config.filter.resampling.effNParticle", nParticles_);
    minUpdatesBeforeResample_ = pt.get("config.filter.resampling.minTimesteps", 1);
    
    landmarkExistencePruningThreshold_ = pt.get("config.filter.prune.threshold", -5.0);

    return true;   
  }

  /** Generate a random trajectory in 2d space */
  void generateTrajectory(int randSeed = 0){

    printf("Generating trajectory with random seed = %d\n", randSeed);
    srand48( randSeed );

    TimeStamp t;
    int seg = 0;
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q = dT_ * Q * dT_;
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

      if( k <= 50 ){
	// No motion
	MotionModel_Odometry2d::TInput::Vec d;
	MotionModel_Odometry2d::TInput::Vec dCovDiag;
	d << 0, 0, 0;
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
	input_k = MotionModel_Odometry2d::TInput(d, dCovDiag.asDiagonal(), t);  
      }

      groundtruth_displacement_.push_back(input_k);
      groundtruth_displacement_.back().setTime(t);

      MotionModel_Odometry2d::TState x_k;
      motionModel.step(x_k, groundtruth_pose_[k-1], input_k, dTimeStamp_);
      groundtruth_pose_.push_back( x_k );
      groundtruth_pose_.back().setTime(t);

    }

  }
  
  /** Generate odometry measurements */
  void generateOdometry(){

    odometry_.reserve( kMax_ );
    MotionModel_Odometry2d::TInput zero;
    MotionModel_Odometry2d::TInput::Vec u0;
    u0.setZero();
    zero.set(u0, 0);
    odometry_.push_back( zero );


    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    MotionModel_Odometry2d motionModel(Q);
    deadReckoning_pose_.reserve( kMax_ );
    deadReckoning_pose_.push_back( groundtruth_pose_[0] );

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
      motionModel.step(p, deadReckoning_pose_[k-1], odometry_[k], dTimeStamp_);
      p.setTime(t);
      deadReckoning_pose_.push_back( p );
    }

  }

  /** Generate landmarks */
  void generateLandmarks(){

    MeasurementModel_RngBrg measurementModel( varzr_, varzb_);
    MeasurementModel_RngBrg::TPose pose;

    groundtruth_landmark_.reserve(nLandmarks_);

    int nLandmarksCreated = 0;
    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nLandmarks_ * nLandmarksCreated){

	MeasurementModel_RngBrg::TPose pose;
	MeasurementModel_RngBrg::TMeasurement measurementToCreateLandmark;
	MeasurementModel_RngBrg::TMeasurement::Vec z;
	double r = drand48() * rangeLimitMax_;
	double b = drand48() * 2 * PI;
	z << r, b;
	measurementToCreateLandmark.set(z);
	MeasurementModel_RngBrg::TLandmark lm;
	
	measurementModel.inverseMeasure( groundtruth_pose_[k], 
					 measurementToCreateLandmark, 
					 lm);

	groundtruth_landmark_.push_back(lm);

	nLandmarksCreated++;
	
      }

    }

  }

  /** Generate landmark measurements */
  void generateMeasurements(){

    MeasurementModel_RngBrg measurementModel( varzr_, varzb_ );
    MeasurementModel_RngBrg::TMeasurement::Mat R;
    measurementModel.getNoise(R);
    measurementModel.config.rangeLimMax_ = rangeLimitMax_;
    measurementModel.config.rangeLimMin_ = rangeLimitMin_;
    measurementModel.config.probabilityOfDetection_ = Pd_;
    measurementModel.config.uniformClutterIntensity_ = c_;
    double meanClutter = measurementModel.clutterIntensityIntegral();
    
    double expNegMeanClutter = exp( -meanClutter );
    double poissonPmf[100];
    double poissonCmf[100];
    double mean_pow_i = 1;
    double i_factorial = 1;
    poissonPmf[0] = expNegMeanClutter;
    poissonCmf[0] = poissonPmf[0];
    for( int i = 1; i < 100; i++){
      mean_pow_i *= meanClutter;
      i_factorial *= i;
      poissonPmf[i] = mean_pow_i / i_factorial * expNegMeanClutter;
      poissonCmf[i] = poissonCmf[i-1] + poissonPmf[i]; 
    }

    lmkFirstObsTime_.resize( groundtruth_landmark_.size());
    for( int m = 0; m < lmkFirstObsTime_.size(); m++ ){
      lmkFirstObsTime_[m] = -1;
    }

    TimeStamp t;

    for( int k = 1; k < kMax_; k++ ){
      
      t += dTimeStamp_;

      groundtruth_pose_[k];
      
      // Real detections
      for( int m = 0; m < groundtruth_landmark_.size(); m++){
	
	bool success;
	MeasurementModel_RngBrg::TMeasurement z_m_k;
	success = measurementModel.sample( groundtruth_pose_[k],
					   groundtruth_landmark_[m],
					   z_m_k);
	if(success){

	  if(z_m_k.get(0) <= rangeLimitMax_ && z_m_k.get(0) >= rangeLimitMin_ && drand48() <= Pd_){
	    z_m_k.setTime(t);
	    z_m_k.setCov(R);
	    measurements_.push_back( z_m_k );
	  }

	  // Observe the first time in which a landmark is observed 
	  if(lmkFirstObsTime_[m] == -1){
	    lmkFirstObsTime_[m] = t.getTimeAsDouble();
	  }
	}

      }

      // False alarms
      double randomNum = drand48();
      int nClutterToGen = 0;
      while( randomNum > poissonCmf[ nClutterToGen ] ){
	nClutterToGen++;
      }
      for( int i = 0; i < nClutterToGen; i++ ){
	
	double r = drand48() * rangeLimitMax_;
	while(r < rangeLimitMin_)
	  r = drand48() * rangeLimitMax_;
	double b = drand48() * 2 * PI - PI;
	MeasurementModel_RngBrg::TMeasurement z_clutter;
	MeasurementModel_RngBrg::TMeasurement::Vec z;
	z << r, b;
	z_clutter.set(z, t);
	measurements_.push_back(z_clutter);
	
      }
      
    }
    
  }

  /** Data Logging */
  void exportSimData(){

    if(!logToFile_)
      return;

    boost::filesystem::path dir(logDirPrefix_);
    boost::filesystem::create_directories(dir);

    boost::filesystem::path cfgFilePathSrc( cfgFileName_ );
    std::string cfgFileDst( logDirPrefix_ );
    cfgFileDst += "simSettings.cfg";
    boost::filesystem::path cfgFilePathDst( cfgFileDst.data() );
    boost::filesystem::copy_file( cfgFilePathSrc, cfgFilePathDst, boost::filesystem::copy_option::overwrite_if_exists);

    TimeStamp t;


    FILE* pGTPoseFile;
    std::string filenameGTPose( logDirPrefix_ );
    filenameGTPose += "gtPose.dat";
    pGTPoseFile = fopen(filenameGTPose.data(), "w");
    MotionModel_Odometry2d::TState::Vec x;
    for(int i = 0; i < groundtruth_pose_.size(); i++){
      groundtruth_pose_[i].get(x, t);
      fprintf( pGTPoseFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), x(0), x(1), x(2));
    }
    fclose(pGTPoseFile);

    FILE* pGTLandmarkFile;
    std::string filenameGTLandmark( logDirPrefix_ );
    filenameGTLandmark += "gtLandmark.dat";
    pGTLandmarkFile = fopen(filenameGTLandmark.data(), "w");
    MeasurementModel_RngBrg::TLandmark::Vec m;
    for(int i = 0; i < groundtruth_landmark_.size(); i++){
      groundtruth_landmark_[i].get(m);
      fprintf( pGTLandmarkFile, "%f   %f   %f\n", m(0), m(1), lmkFirstObsTime_[i]);
    }
    fclose(pGTLandmarkFile);

    FILE* pOdomFile;
    std::string filenameOdom( logDirPrefix_ );
    filenameOdom += "odometry.dat";
    pOdomFile = fopen(filenameOdom.data(),"w");
    MotionModel_Odometry2d::TInput::Vec u;
    for(int i = 0; i < odometry_.size(); i++){
      odometry_[i].get(u, t);
      fprintf( pOdomFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), u(0), u(1), u(2));
    }
    fclose(pOdomFile);

    FILE* pMeasurementFile;
    std::string filenameMeasurement( logDirPrefix_ );
    filenameMeasurement += "measurement.dat";
    pMeasurementFile = fopen(filenameMeasurement.data(), "w");
    MeasurementModel_RngBrg::TMeasurement::Vec z;
    for(int i = 0; i < measurements_.size(); i++){
      measurements_[i].get(z, t);
      fprintf( pMeasurementFile, "%f   %f   %f\n", t.getTimeAsDouble(), z(0), z(1) );
    }
    fclose(pMeasurementFile);

    FILE* pDeadReckoningFile;
    std::string filenameDeadReckoning( logDirPrefix_ );
    filenameDeadReckoning += "deadReckoning.dat";
    pDeadReckoningFile = fopen(filenameDeadReckoning.data(), "w");
    MotionModel_Odometry2d::TState::Vec odo;
    for(int i = 0; i < deadReckoning_pose_.size(); i++){
      deadReckoning_pose_[i].get(odo, t);
      fprintf( pDeadReckoningFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), odo(0), odo(1), odo(2));
    }
    fclose(pDeadReckoningFile);

  }

  /** FastSLAM Setup */
  void setupRBPHDFilter(){
    
    pFilter_ = new FastSLAM<MotionModel_Odometry2d,
			    StaticProcessModel<Landmark2d>,
			    MeasurementModel_RngBrg,
			    KalmanFilter_RngBrg>( nParticles_ );

    double dt = dTimeStamp_.getTimeAsDouble();

    // configure robot motion model (only need to set once since timesteps are constant)
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= (pNoiseInflation_ * dt * dt);
    pFilter_->getProcessModel()->setNoise(Q);

    // configure landmark process model (only need to set once since timesteps are constant)
    Landmark2d::Mat Q_lm;
    Q_lm << varlmx_, 0, 0, varlmy_;
    Q_lm = Q_lm * dt * dt;
    pFilter_->getLmkProcessModel()->setNoise(Q_lm);

    // configure measurement model
    MeasurementModel_RngBrg::TMeasurement::Mat R;
    R << varzr_, 0, 0, varzb_;
    R *= zNoiseInflation_;
    pFilter_->getMeasurementModel()->setNoise(R);
    pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
    pFilter_->getMeasurementModel()->config.uniformClutterIntensity_ = c_;
    pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
    pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
    pFilter_->getMeasurementModel()->config.rangeLimBuffer_ = rangeLimitBuffer_;

    // configure the Kalman filter for landmark updates
    pFilter_->getKalmanFilter()->config.rangeInnovationThreshold_ = innovationRangeThreshold_;
    pFilter_->getKalmanFilter()->config.bearingInnovationThreshold_ = innovationBearingThreshold_;

    // configure the filter
    pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
    pFilter_->config.minUpdatesBeforeResample_ = minUpdatesBeforeResample_;
    pFilter_->config.minLogMeasurementLikelihood_ = minLogMeasurementLikelihood_;
    pFilter_->config.maxNDataAssocHypotheses_ = maxNDataAssocHypotheses_;
    pFilter_->config.maxDataAssocLogLikelihoodDiff_ = maxDataAssocLogLikelihoodDiff_;
    pFilter_->config.mapExistencePruneThreshold_ = landmarkExistencePruningThreshold_;
    pFilter_->config.landmarkExistencePrior_ = 0.5;
  }

  /** Run the simulator */
  void run(){
    
    printf("Running simulation\n\n");

    //////// Initialization at first timestep //////////

    FILE* pParticlePoseFile;
    if(logToFile_){
      std::string filenameParticlePoseFile( logDirPrefix_ );
      filenameParticlePoseFile += "particlePose.dat";
      pParticlePoseFile = fopen(filenameParticlePoseFile.data(), "w");
    }
    FILE* pLandmarkEstFile;
    if(logToFile_){
      std::string filenameLandmarkEstFile( logDirPrefix_ );
      filenameLandmarkEstFile += "landmarkEst.dat";
      pLandmarkEstFile = fopen(filenameLandmarkEstFile.data(), "w");
    }
    MotionModel_Odometry2d::TState x_i;
    int zIdx = 0;

    if(logToFile_){
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	x_i = *(pFilter_->getParticleSet()->at(i));
	fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   1.0\n", 0.0, i, x_i.get(0), x_i.get(1), x_i.get(2));
      }  
    }


    /////////// Run simulator from k = 1 to kMax_ /////////

    TimeStamp time;

    for(int k = 1; k < kMax_; k++){

      time += dTimeStamp_;
 
      if( k % 100 == 0)
	printf("k = %d\n", k);
 
      
      ////////// Prediction Step //////////

      // configure robot motion model ( not necessary since in simulation, timesteps are constant)
      // MotionModel_Odometry2d::TState::Mat Q;
      // Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
      // Q *= (pNoiseInflation_ * dt * dt);
      // pFilter_->getProcessModel()->setNoise(Q);

      // configure landmark process model ( not necessary since in simulation, timesteps are constant)
      // Landmark2d::Mat Q_lm;
      // Q_lm << varlmx_, 0, 0, varlmy_;
      // Q_lm = Q_lm * dt * dt;
      // pFilter_->getLmkProcessModel()->setNoise(Q_lm);

      pFilter_->predict( odometry_[k], dTimeStamp_ );
      
      if( k <= 100){ // hold still for the first 100 timesteps
	for( int i = 0; i < nParticles_; i++)
	  pFilter_->setParticlePose(i, groundtruth_pose_[k]);
      }

      // Prepare measurement vector for update
      std::vector<MeasurementModel_RngBrg::TMeasurement> Z;
      TimeStamp kz = measurements_[ zIdx ].getTime();
      while( kz == time ){
	Z.push_back( measurements_[zIdx] );
	zIdx++;
	if(zIdx >= measurements_.size())
	  break;
	kz = measurements_[ zIdx ].getTime();
      }

      ////////// Update Step //////////
      pFilter_->update(Z);

      // Log particle poses
      if(logToFile_){
	int i_w_max = 0;
	double w_max = 0;
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  x_i = *(pFilter_->getParticleSet()->at(i));
	  double w = pFilter_->getParticleSet()->at(i)->getWeight();
	  if(w > w_max){
	    i_w_max = i;
	    w_max = w;
	  }
	  fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   %f\n", time.getTimeAsDouble(), i, x_i.get(0), x_i.get(1), x_i.get(2), w);
	}
	fprintf( pParticlePoseFile, "\n");
      }

      // Log landmark estimates
      if(logToFile_){

	int mapSize = pFilter_->getGMSize(i_w_max);
	for( int m = 0; m < mapSize; m++ ){
	  MeasurementModel_RngBrg::TLandmark::Vec u;
	  MeasurementModel_RngBrg::TLandmark::Mat S;
	  double w;
	  pFilter_->getLandmark(i_w_max, m, u, S, w);
	    
	  fprintf( pLandmarkEstFile, "%f   %d   ", time.getTimeAsDouble(), i_w_max);
	  fprintf( pLandmarkEstFile, "%f   %f      ", u(0), u(1));
	  fprintf( pLandmarkEstFile, "%f   %f   %f", S(0,0), S(0,1), S(1,1));
	  fprintf( pLandmarkEstFile, "   %f\n", 1 - 1/(1 + exp(w)) );
	  
	}
      }
    }

        printf("Elapsed Timing Information [nsec]\n");
    printf("Prediction    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->predict_wall, pFilter_->getTimingInfo()->predict_cpu);
    printf("Map Update    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapUpdate_wall, pFilter_->getTimingInfo()->mapUpdate_cpu);
    printf("Resampling    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->particleResample_wall, pFilter_->getTimingInfo()->particleResample_cpu);
    printf("Total         -- wall: %lld   cpu: %lld\n",
	   pFilter_->getTimingInfo()->predict_wall +
	   pFilter_->getTimingInfo()->mapUpdate_wall +
	   pFilter_->getTimingInfo()->particleResample_wall,
	   pFilter_->getTimingInfo()->predict_cpu +
	   pFilter_->getTimingInfo()->mapUpdate_cpu +
	   pFilter_->getTimingInfo()->particleResample_cpu);
    printf("\n");

    if(logToFile_){
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }
  }

private:

  const char* cfgFileName_;

  int kMax_; /**< number of timesteps */
  double dT_; /**< duration of timestep */
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

  // Filters
  KalmanFilter_RngBrg kf_;
  FastSLAM<MotionModel_Odometry2d, 
	   StaticProcessModel<Landmark2d>,
	   MeasurementModel_RngBrg,  
	   KalmanFilter_RngBrg> *pFilter_; 
  int nParticles_;
  double pNoiseInflation_;
  double zNoiseInflation_;
  double innovationRangeThreshold_;
  double innovationBearingThreshold_;
  double effNParticleThreshold_;
  int minUpdatesBeforeResample_;
  double minLogMeasurementLikelihood_;
  int maxNDataAssocHypotheses_;
  double maxDataAssocLogLikelihoodDiff_;
  double landmarkExistencePruningThreshold_;
  
  bool logToFile_;
  std::string logDirPrefix_;
};




int main(int argc, char* argv[]){   

  int initRandSeed = 0;
  const char* cfgFileName = "cfg/fastslam2dSim.cfg";
  if( argc >= 2 ){
    initRandSeed = boost::lexical_cast<int>(argv[1]);
  }
  if( argc >= 3 ){
    cfgFileName = argv[2];
  }

  Simulator_FastSLAM_2d sim;
  if( !sim.readConfigFile( cfgFileName ) ){
    return -1;
  }
  sim.generateTrajectory( initRandSeed );
  sim.generateLandmarks();

  sim.generateOdometry();
  sim.generateMeasurements();
  sim.exportSimData();
  sim.setupRBPHDFilter();

  srand48( time(NULL) );

  boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Simulation run time: %ws\n");

  sim.run();

  delete timer;

  return 0;

}
