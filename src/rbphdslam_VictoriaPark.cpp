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
#include "ProcessModel_Odometry2D.hpp"
#include "RBPHDFilter.hpp"
#include "KalmanFilter_RngBrg.hpp"
#include <stdio.h>
#include <string>

using namespace rfs;

/**
 * \class RBPHDSLAM_VictoriaPark
 * \brief Run the RB-PHD-SLAM Algorithm on the Victoria Park Dataset.
 * \author Keith Leung
 */
class RBPHDSLAM_VictoriaPark{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef RBPHDFilter<MotionModel_Odometry2d, 
		      StaticProcessModel<Landmark2d>,
		      MeasurementModel_RngBrg,  
		      KalmanFilter_RngBrg> SLAM_Filter;

  RBPHDSLAM_VictoriaPark(){
    pFilter_ = NULL;
  }
  
  ~RBPHDSLAM_VictoriaPark(){
    
    if(pFilter_ != NULL){
      delete pFilter_;
    }

  }

  /** Read the simulator configuration file */
  bool readConfigFile(const char* fileName){
    
    cfgFileName_ = fileName;

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(fileName, pt);

    std::string dataDir = pt.get<std::string>("config.dataset.directory");
    if(*(dataDir.rbegin()) != '/'){
      dataDir += '/';
    }
    dataFileGPS_ = dataDir + pt.get<std::string>("config.dataset.filename.gps");
    dataFileLidar_ = dataDir + pt.get<std::string>("config.dataset.filename.lidar");
    dataFileInput_ = dataDir + pt.get<std::string>("config.dataset.filename.input");
    dataFileSensorManager_ = dataDir + pt.get<std::string>("config.dataset.filename.manager");

    logToFile_ = false;
    if( pt.get("config.logging.logToFile", 0) == 1 )
      logToFile_ = true;
    logDirPrefix_ = pt.get<std::string>("config.logging.logDirPrefix", "./");

    vardx_ = pt.get<double>("config.process.vardx");
    vardy_ = pt.get<double>("config.process.vardy");
    vardz_ = pt.get<double>("config.process.vardz");    

    varlmx_ = pt.get<double>("config.landmarks.varlmx");
    varlmy_ = pt.get<double>("config.landmarks.varlmy");

    rangeLimitMax_ = pt.get<double>("config.measurements.rangeLimitMax");
    rangeLimitMin_ = pt.get<double>("config.measurements.rangeLimitMin");
    bearingLimitMax_ = pt.get<double>("config.measurements.bearingLimitMax");
    bearingLimitMin_ = pt.get<double>("config.measurements.bearingLimitMin");
    rangeLimitBuffer_ = pt.get<double>("config.measurements.rangeLimitBuffer");
    Pd_ = pt.get<double>("config.measurements.probDetection");
    c_ = pt.get<double>("config.measurements.clutterIntensity");
    varzr_ = pt.get<double>("config.measurements.varzr");
    varzb_ = pt.get<double>("config.measurements.varzb");
    varzd_ = pt.get<double>("config.measurements.varzd");

    nParticles_ = pt.get("config.filter.nParticles", 200);

    pNoiseInflation_ = pt.get("config.filter.predict.processNoiseInflationFactor", 1.0);
    birthGaussianWeight_ = pt.get("config.filter.predict.birthGaussianWeight", 0.01);

    zNoiseInflation_ = pt.get("config.filter.update.measurementNoiseInflationFactor", 1.0);
    innovationRangeThreshold_ = pt.get<double>("config.filter.update.KalmanFilter.innovationThreshold.range");
    innovationBearingThreshold_ = pt.get<double>("config.filter.update.KalmanFilter.innovationThreshold.bearing");
    newGaussianCreateInnovMDThreshold_ = pt.get<double>("config.filter.update.GaussianCreateInnovMDThreshold");

    importanceWeightingEvalPointCount_ = pt.get("config.filter.weighting.nEvalPt", 15);
    importanceWeightingEvalPointGuassianWeight_ = pt.get("config.filter.weighting.minWeight", 0.75);
    importanceWeightingMeasurementLikelihoodMDThreshold_ = pt.get("config.filter.weighting.threshold", 3.0);
    useClusterProcess_ = false;
    if( pt.get("config.filter.weighting.useClusterProcess", 0) == 1 )
      useClusterProcess_ = true;

    effNParticleThreshold_ = pt.get("config.filter.resampling.effNParticle", nParticles_);
    minUpdatesBeforeResample_ = pt.get("config.filter.resampling.minTimesteps", 1);
    
    gaussianMergingThreshold_ = pt.get<double>("config.filter.merge.threshold");
    gaussianMergingCovarianceInflationFactor_ = pt.get("config.filter.merge.covInflationFactor", 1.0);
    
    gaussianPruningThreshold_ = pt.get("config.filter.prune.threshold", birthGaussianWeight_);

    // Copy config file to logDir
    if(logToFile_){
      boost::filesystem::path dir(logDirPrefix_);
      boost::filesystem::create_directories(dir);
      boost::filesystem::path cfgFilePathSrc( cfgFileName_ );
      std::string cfgFileDst( logDirPrefix_ );
      cfgFileDst += "settings.cfg";
      boost::filesystem::path cfgFilePathDst( cfgFileDst.data() );
      boost::filesystem::copy_file( cfgFilePathSrc, cfgFilePathDst, boost::filesystem::copy_option::overwrite_if_exists);
    }

    return true;   
  }

  /** \brief Import dataset from files */
  void readData(){

    // Read sensor manager log
    // \TODO define log
    std::cout << "Reading process input file: " << dataFileSensorManager_ << std::endl;

    // Read Process Model Inputs
    // \TODO define motion model
    // odometry_.push_back()
    std::cout << "Reading process input file: " << dataFileInput_ << std::endl;

    // Read Lidar measurements
    // \TODO define measurement model 
    // measurements_.push_back( z_m_k );
    std::cout << "Reading process input file: " << dataFileLidar_ << std::endl;
    
    // Ground truth GPS
    //groundtruth_pose_.push_back( x_k );
    //groundtruth_pose_.back().setTime(t);
    std::cout << "Reading process input file: " << dataFileGPS_ << std::endl;

  }

  /** \brief Peform dead reckoning with process input data */
  void deadReckoning(){
    /*
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    MotionModel_Odometry2d motionModel(Q);
    deadReckoning_pose_.reserve( kMax_ );
    deadReckoning_pose_.push_back( groundtruth_pose_[0] );
    */
    
    if(logToFile_){
      FILE* pDeadReckoningFile;
      std::string filenameDeadReckoning( logDirPrefix_ );
      filenameDeadReckoning += "deadReckoning.dat";
      pDeadReckoningFile = fopen(filenameDeadReckoning.data(), "w");
      MotionModel_Odometry2d::TState::Vec odo;
      TimeStamp t;
      for(int i = 0; i < deadReckoning_pose_.size(); i++){
	deadReckoning_pose_[i].get(odo, t);
	fprintf( pDeadReckoningFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), odo(0), odo(1), odo(2));
      }
      fclose(pDeadReckoningFile);
    }

  }

  /** RB-PHD Filter Setup */
  void setupRBPHDFilter(){
    
    pFilter_ = new RBPHDFilter<MotionModel_Odometry2d,
			       StaticProcessModel<Landmark2d>,
			       MeasurementModel_RngBrg,
			       KalmanFilter_RngBrg>( nParticles_ );

    //double dt = dTimeStamp_.getTimeAsDouble();

    // configure robot motion model (only need to set once since timesteps are constant)
    MotionModel_Odometry2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q *= (pNoiseInflation_);
    //Q *= (pNoiseInflation_ * dt * dt);
    pFilter_->getProcessModel()->setNoise(Q);

    // configure landmark process model (only need to set once since timesteps are constant)
    Landmark2d::Mat Q_lm;
    Q_lm << varlmx_, 0, 0, varlmy_;
    // Q_lm = Q_lm * dt * dt;
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
    pFilter_->config.birthGaussianWeight_ = birthGaussianWeight_;
    pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
    pFilter_->config.minUpdatesBeforeResample_ = minUpdatesBeforeResample_;
    pFilter_->config.newGaussianCreateInnovMDThreshold_ = newGaussianCreateInnovMDThreshold_;
    pFilter_->config.importanceWeightingMeasurementLikelihoodMDThreshold_ = importanceWeightingMeasurementLikelihoodMDThreshold_;
    pFilter_->config.importanceWeightingEvalPointCount_ = importanceWeightingEvalPointCount_;
    pFilter_->config.importanceWeightingEvalPointGuassianWeight_ = importanceWeightingEvalPointGuassianWeight_;
    pFilter_->config.gaussianMergingThreshold_ = gaussianMergingThreshold_;
    pFilter_->config.gaussianMergingCovarianceInflationFactor_ = gaussianMergingCovarianceInflationFactor_;
    pFilter_->config.gaussianPruningThreshold_ = gaussianPruningThreshold_;
    pFilter_->config.useClusterProcess_ = useClusterProcess_;
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

    
    if(logToFile_){
      MotionModel_Odometry2d::TState x_i;
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	pFilter_->getParticleSet()->at(i)->getPose(x_i);
	fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   1.0\n", 0.0, i, x_i.get(0), x_i.get(1), x_i.get(2));
      }   
    }
   
    int zIdx = 0;


  
    /////////// Run simulator from k = 1 to kMax_ /////////
    /*
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
      
      if( k <= 100){
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
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  pFilter_->getParticleSet()->at(i)->getPose(x_i);
	  double w = pFilter_->getParticleSet()->at(i)->getWeight();
	  fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   %f\n", time.getTimeAsDouble(), i, x_i.get(0), x_i.get(1), x_i.get(2), w);
	}
	fprintf( pParticlePoseFile, "\n");
      }

      // Log landmark estimates
      if(logToFile_){
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  int mapSize = pFilter_->getGMSize(i);
	  for( int m = 0; m < mapSize; m++ ){
	    MeasurementModel_RngBrg::TLandmark::Vec u;
	    MeasurementModel_RngBrg::TLandmark::Mat S;
	    double w;
	    pFilter_->getLandmark(i, m, u, S, w);
	    
	    fprintf( pLandmarkEstFile, "%f   %d   ", time.getTimeAsDouble(), i);
	    fprintf( pLandmarkEstFile, "%f   %f      ", u(0), u(1));
	    fprintf( pLandmarkEstFile, "%f   %f   %f", S(0,0), S(0,1), S(1,1));
	    fprintf( pLandmarkEstFile, "   %f\n", w );
	  }
	}
      }

    }*/
    /*
    printf("Elapsed Timing Information [nsec]\n");
    printf("Prediction    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->predict_wall, pFilter_->getTimingInfo()->predict_cpu);
    printf("Map Update    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapUpdate_wall, pFilter_->getTimingInfo()->mapUpdate_cpu);
    printf("Map Update KF -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapUpdate_kf_wall, pFilter_->getTimingInfo()->mapUpdate_kf_cpu);
    printf("Weighting     -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->particleWeighting_wall, pFilter_->getTimingInfo()->particleWeighting_cpu);
    printf("Map Merge     -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapMerge_wall, pFilter_->getTimingInfo()->mapMerge_cpu);
    printf("Map Prune     -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->mapPrune_wall, pFilter_->getTimingInfo()->mapPrune_cpu);
    printf("Resampling    -- wall: %lld   cpu: %lld\n", 
	   pFilter_->getTimingInfo()->particleResample_wall, pFilter_->getTimingInfo()->particleResample_cpu);
    printf("Total         -- wall: %lld   cpu: %lld\n",
	   pFilter_->getTimingInfo()->predict_wall +
	   pFilter_->getTimingInfo()->mapUpdate_wall +
	   pFilter_->getTimingInfo()->particleWeighting_wall +
	   pFilter_->getTimingInfo()->mapMerge_wall +
	   pFilter_->getTimingInfo()->mapPrune_wall +
	   pFilter_->getTimingInfo()->particleResample_wall,
	   pFilter_->getTimingInfo()->predict_cpu +
	   pFilter_->getTimingInfo()->mapUpdate_cpu +
	   pFilter_->getTimingInfo()->particleWeighting_cpu +
	   pFilter_->getTimingInfo()->mapMerge_cpu +
	   pFilter_->getTimingInfo()->mapPrune_cpu + 
	   pFilter_->getTimingInfo()->particleResample_cpu);
    printf("\n");
    */

    if(logToFile_){
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }
  }

private:

  const char* cfgFileName_;

  // Dataset
  std::string dataFileGPS_;
  std::string dataFileLidar_;
  std::string dataFileInput_;
  std::string dataFileSensorManager_;

  // Process model
  double vardx_;
  double vardy_;
  double vardz_;
  std::vector<MotionModel_Odometry2d::TState> groundtruth_pose_;
  std::vector<MotionModel_Odometry2d::TInput> odometry_;
  std::vector<MotionModel_Odometry2d::TState> deadReckoning_pose_;

  // Landmarks 
  double varlmx_;
  double varlmy_;

  // Range-Bearing Measurements
  double rangeLimitMax_;
  double rangeLimitMin_;
  double bearingLimitMax_;
  double bearingLimitMin_;
  double rangeLimitBuffer_;
  double Pd_;
  double c_;
  double varzr_;
  double varzb_;
  double varzd_;
  std::vector<MeasurementModel_RngBrg::TMeasurement> measurements_;

  // Filters
  KalmanFilter_RngBrg kf_;
  SLAM_Filter *pFilter_; 
  int nParticles_;
  double pNoiseInflation_;
  double zNoiseInflation_;
  double innovationRangeThreshold_;
  double innovationBearingThreshold_;
  double birthGaussianWeight_;
  double newGaussianCreateInnovMDThreshold_;
  double importanceWeightingMeasurementLikelihoodMDThreshold_;
  double importanceWeightingEvalPointGuassianWeight_;
  double effNParticleThreshold_;
  int minUpdatesBeforeResample_;
  double gaussianMergingThreshold_;
  double gaussianMergingCovarianceInflationFactor_;
  double gaussianPruningThreshold_;
  int importanceWeightingEvalPointCount_;
  bool useClusterProcess_;

  bool logToFile_;

public:
  std::string logDirPrefix_;
};



int main(int argc, char* argv[]){

  int initRandSeed = 0;
  const char* cfgFileName = "cfg/rbphdslam_VictoriaPark.xml";
  if( argc >= 2 ){
    initRandSeed = boost::lexical_cast<int>(argv[1]);
  }
  if( argc >= 3 ){
    cfgFileName = argv[2];
  }
  std::cout << "rbphdslam_VictoraPark [randomSeed] [cfgFile]\n";
  std::cout << "[randomSeed] = " << initRandSeed << std::endl;
  std::cout << "[cfgFile] = " << cfgFileName << std::endl;

  RBPHDSLAM_VictoriaPark slam;

  // Read config file
  if( !slam.readConfigFile( cfgFileName ) ){
    std::cout << "[Error] Unable to read config file: " << cfgFileName << std::endl;
    return -1;
  }

  // Read files from the dataset
  slam.readData();

  /*

  sim.setupRBPHDFilter();

  srand48( time(NULL) );

  boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Run time: %ws\n");

  sim.run(); 

  delete timer;
  */

  return 0;

}
