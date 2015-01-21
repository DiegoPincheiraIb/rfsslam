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

#include <assert.h>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "ProcessModel_Ackerman2D.hpp"
#include "RBPHDFilter.hpp"
#include "KalmanFilter_RngBrg.hpp"
#include <stdio.h>
#include <string>
#include <sstream>

using namespace rfs;

/**
 * \class RBPHDSLAM_VictoriaPark
 * \brief Run the RB-PHD-SLAM Algorithm on the Victoria Park Dataset.
 * \author Keith Leung
 */
class RBPHDSLAM_VictoriaPark{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  typedef RBPHDFilter<MotionModel_Ackerman2d, 
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
    dataFileDetection_ = dataDir + pt.get<std::string>("config.dataset.filename.detection");
    dataFileLidar_ = dataDir + pt.get<std::string>("config.dataset.filename.lidar");
    dataFileInput_ = dataDir + pt.get<std::string>("config.dataset.filename.input");
    dataFileSensorManager_ = dataDir + pt.get<std::string>("config.dataset.filename.manager");

    logToFile_ = false;
    if( pt.get("config.logging.logToFile", 0) == 1 )
      logToFile_ = true;
    logDirPrefix_ = pt.get<std::string>("config.logging.logDirPrefix", "./");


    ackerman_h_ = pt.get<double>("config.process.AckermanModel.rearWheelOffset");
    ackerman_l_ = pt.get<double>("config.process.AckermanModel.frontToRearDist");
    ackerman_dx_ = pt.get<double>("config.process.AckermanModel.sensorOffset_x");
    ackerman_dy_ = pt.get<double>("config.process.AckermanModel.sensorOffset_y");
    var_uv_ = pt.get<double>("config.process.varuv");
    var_ur_ = pt.get<double>("config.process.varur");

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

  struct SensorManagerMsg{
    TimeStamp t;
    enum Type {GPS=1, Input=2, Lidar=3};
    Type sensorType;
    uint idx;
  };

  struct LidarScanMsg{
    TimeStamp t;
    double scan[361];
  };

  /** \brief Import dataset from files */
  void readData(){

    std::string msgLine;

    // Read sensor manager log
    std::cout << "Reading input file: " << dataFileSensorManager_ << std::endl;
    std::ifstream file_sensorManager( dataFileSensorManager_.c_str() );
    assert( file_sensorManager.is_open() );
    while( std::getline( file_sensorManager, msgLine ) ){
      SensorManagerMsg msg;
      double time;
      int type;
      std::stringstream ss( msgLine );
      ss >> time >> type >> msg.idx;
      msg.t.setTime(time);
      msg.sensorType = (SensorManagerMsg::Type)type;
      msg.idx--;
      sensorManagerMsgs_.push_back(msg);
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << msg.t.getTimeAsDouble() 
      //	<< std::setw(10) << (int)(msg.sensorType) 
      //	<< std::setw(10) << msg.idx << std::endl;
    }
    file_sensorManager.close();

    // Read Process Model Inputs    
    std::cout << "Reading input file: " << dataFileInput_ << std::endl;
    std::ifstream file_input( dataFileInput_.c_str() );
    assert( file_input.is_open() );
    while( std::getline( file_input, msgLine ) ){
      double time, vel, steer;
      std::stringstream ss( msgLine );
      ss >> time >> vel >> steer;
      SLAM_Filter::TInput::Vec uVec;
      uVec << vel, steer;
      SLAM_Filter::TInput u(uVec, TimeStamp(time)); 
      motionInputs_.push_back( u );
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << u.getTime().getTimeAsDouble() 
      //	<< std::setw(10) << u[0] 
      //	<< std::setw(10) << std::fixed << std::setprecision(4) << u[1] << std::endl;
    }
    file_input.close();

    // Read Lidar detections
    std::cout << "Reading input file: " << dataFileDetection_ << std::endl;
    std::ifstream file_measurements( dataFileDetection_.c_str() );
    assert( file_measurements.is_open() );
    while( std::getline( file_measurements, msgLine ) ){
      double time, range, bearing, diameter;
      std::stringstream ss(msgLine);
      ss >> time >> range >> bearing >> diameter;
      // \TODO
      /* uncomment after new measurement model is defined 
      SLAM_Filter::TMeasurement::Vec zVec;
      zVec << range, bearing, diameter;
      SLAM_Filter::TMeasurement z(zVec, TimeStamp(time));
      measurements_.push_back(z);
      std::cout << std::setw(10) << std::fixed << std::setprecision(3) << z.getTime().getTimeAsDouble() 
      	<< std::setw(10) << std::fixed << std::setprecision(5) << z[0] 
	<< std::setw(10) << std::fixed << std::setprecision(5) << z[1]
      	<< std::setw(10) << std::fixed << std::setprecision(5) << z[2] << std::endl;
      */
    }
    file_measurements.close();

    // Read Lidar raw scans
    std::cout << "Reading input file: " << dataFileLidar_ << std::endl;
    std::ifstream file_lidar( dataFileLidar_.c_str() );
    assert( file_lidar.is_open() );
    double t = -1;
    LidarScanMsg msg;
    file_lidar >> t; 
    while(t >= 0){
      msg.t.setTime(t);
      for(int i = 0; i < 361; i++){
	file_lidar >> msg.scan[i];
      }
      lidarScans_.push_back(msg);
      t = -1;
      file_lidar >> t;
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << msg.t.getTimeAsDouble() 
      //		<< std::setw(10) << std::fixed << std::setprecision(5) << msg.scan[0] 
      //		<< std::setw(10) << std::fixed << std::setprecision(5) << msg.scan[360] << std::endl;
    }
    file_lidar.close();

    
    // Ground truth GPS
    //groundtruth_pose_.push_back( x_k );
    //groundtruth_pose_.back().setTime(t);
    std::cout << "Reading input file: " << dataFileGPS_ << std::endl;
    std::ifstream file_gps( dataFileGPS_.c_str() );
    assert( file_gps.is_open() );
    while( std::getline( file_gps, msgLine ) ){
      double time, x, y;
      std::stringstream ss(msgLine);
      ss >> time >> x >> y;
      Position2d::Vec pVec;
      pVec << x, y;
      Position2d p(pVec, TimeStamp(time));
      groundtruth_pos_.push_back(p);
      //std::cout << std::setw(10) << std::fixed << std::setprecision(3) << p.getTime().getTimeAsDouble() 
      //	<< std::setw(10) << std::fixed << std::setprecision(2) << p[0] 
      //	<< std::setw(10) << std::fixed << std::setprecision(2) << p[1] << std::endl;
    }
    file_gps.close();

  }

  /** \brief Peform dead reckoning with process input data */
  void deadReckoning(){
    /*
    MotionModel_Ackerman2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    MotionModel_Ackerman2d motionModel(Q);
    deadReckoning_pose_.reserve( kMax_ );
    deadReckoning_pose_.push_back( groundtruth_pose_[0] );
   
    
    if(logToFile_){
      FILE* pDeadReckoningFile;
      std::string filenameDeadReckoning( logDirPrefix_ );
      filenameDeadReckoning += "deadReckoning.dat";
      pDeadReckoningFile = fopen(filenameDeadReckoning.data(), "w");
      MotionModel_Ackerman2d::TState::Vec odo;
      TimeStamp t;
      for(int i = 0; i < deadReckoning_pose_.size(); i++){
	deadReckoning_pose_[i].get(odo, t);
	fprintf( pDeadReckoningFile, "%f   %f   %f   %f\n", t.getTimeAsDouble(), odo(0), odo(1), odo(2));
      }
      fclose(pDeadReckoningFile);
    }

    */

  }

  /** RB-PHD Filter Setup */
  void setupRBPHDFilter(){

    pFilter_ = new SLAM_Filter( nParticles_ );

    // Configure process model
    pFilter_->getProcessModel()->setAckermanParams( ackerman_h_, ackerman_l_, ackerman_dx_, ackerman_dy_);
  
    // configure measurement model
    /*
    MeasurementModel_RngBrg::TMeasurement::Mat R;
    R << varzr_, 0, 0, varzb_;
    R *= zNoiseInflation_;
    pFilter_->getMeasurementModel()->setNoise(R);
    pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
    pFilter_->getMeasurementModel()->config.uniformClutterIntensity_ = c_;
    pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
    pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
    pFilter_->getMeasurementModel()->config.rangeLimBuffer_ = rangeLimitBuffer_;
    */

    // configure the filter

    pFilter_->getKalmanFilter()->config.rangeInnovationThreshold_ = innovationRangeThreshold_;
    pFilter_->getKalmanFilter()->config.bearingInnovationThreshold_ = innovationBearingThreshold_;

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

  /** \brief Process the data and peform SLAM */
  void run(){

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
      MotionModel_Ackerman2d::TState x_i;
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	pFilter_->getParticleSet()->at(i)->getPose(x_i);
	fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   1.0\n", 
		 sensorManagerMsgs_[0].t.getTimeAsDouble(), i, x_i[0], x_i[1], x_i[2]);
      }   
    }
   
    // Process all sensor messages sequentially
    bool isInInitialStationaryState = true;
    TimeStamp t_km(0); // last sensor msg time
    SLAM_Filter::TInput u_km; // last process input
    SLAM_Filter::TInput::Cov u_km_cov;
    u_km_cov << var_uv_, 0, 0, var_ur_;
    u_km.setCov( u_km_cov );
    u_km.setTime( t_km );
    Landmark2d::Cov Q_m_k; // landmark process model additive noise
    for(uint k = 0; k < sensorManagerMsgs_.size(); k++ ){

      if( k % 100 == 0){
	std::cout << "Sensor messages processed: " << k << "/" << sensorManagerMsgs_.size()-1 << std::endl;
      }

      if(sensorManagerMsgs_[k].sensorType == SensorManagerMsg::Input){

	TimeStamp t_k = sensorManagerMsgs_[k].t;
	TimeStamp dt = t_k - t_km;

	Q_m_k << varlmx_, 0, 0, varlmy_;
	Q_m_k = Q_m_k * dt.getTimeAsDouble() * dt.getTimeAsDouble();
	pFilter_->getLmkProcessModel()->setNoise(Q_m_k);

	if(isInInitialStationaryState){
	  pFilter_->predict( u_km, dt, false, false ); // this basically makes all initial particles sit still
	}else{
	  pFilter_->predict( u_km, dt, false, true ); // true for use noise from u_km
	}
	
	u_km = motionInputs_[ sensorManagerMsgs_[k].idx ];
	if(u_km[0] != 0){
	  isInInitialStationaryState = false;
	}

      }else if(sensorManagerMsgs_[k].sensorType == SensorManagerMsg::Lidar){

	TimeStamp t_k = sensorManagerMsgs_[k].t;
	TimeStamp dt = t_k - t_km;

	// Propagate particles up to lidar scan time

	Q_m_k << varlmx_, 0, 0, varlmy_;
	Q_m_k = Q_m_k * dt.getTimeAsDouble() * dt.getTimeAsDouble();
	pFilter_->getLmkProcessModel()->setNoise(Q_m_k);

	if(isInInitialStationaryState){
	  pFilter_->predict( u_km, dt, false, false ); // this basically makes all initial particles sit still
	}else{
	  pFilter_->predict( u_km, dt, false, true ); // true for use noise from u_km
	}

	// Update particles with lidar scan data

	// \todo

	// Log data

	// \todo
	
      }

      t_km = sensorManagerMsgs_[k].t;

    }
    
    
  
    /////////// Run simulator from k = 1 to kMax_ /////////
    /*
    TimeStamp time;

    for(int k = 1; k < kMax_; k++){


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

    if(logToFile_){
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }
  }

private:

  const char* cfgFileName_;

  // Dataset
  std::string dataFileGPS_;
  std::string dataFileDetection_;
  std::string dataFileLidar_;
  std::string dataFileInput_;
  std::string dataFileSensorManager_;
  std::vector<SensorManagerMsg> sensorManagerMsgs_;

  // Process model
  double ackerman_h_;
  double ackerman_l_;
  double ackerman_dx_;
  double ackerman_dy_;
  double var_uv_;
  double var_ur_;
  std::vector<Position2d> groundtruth_pos_;
  std::vector<MotionModel_Ackerman2d::TInput> motionInputs_;
  std::vector<MotionModel_Ackerman2d::TState> deadReckoning_pose_;

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
  std::vector<LidarScanMsg> lidarScans_;

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

  slam.readData();
  slam.setupRBPHDFilter();

  srand48( time(NULL) );
  boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Run time: %ws\n");

  //sim.run(); 

  delete timer;
 
  return 0;

}
