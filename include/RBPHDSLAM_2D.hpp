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
#include "RBPHDFilter.hpp"
#include "KalmanFilter.hpp"
#include <stdio.h>
#include <string>
#include <sstream>
#include <sys/ioctl.h>
#include <type_traits>

#ifdef _PERFTOOLS_CPU
#include <gperftools/profiler.h>
#endif
#ifdef _PERFTOOLS_HEAP
#include <gperftools/heap-profiler.h>
#endif

using namespace rfs;

/**
 * \class RBPHDSLAM_2D
 * \brief A 2d SLAM class using the RB-PHD Filter, it can either simulate a dataset or read an isam style dataset (Support for other 2d datasets should be added here if possible).
 * \author Keith Leung, Felipe Inostroza
 */
template<class MeasurementModel_2D>
class RBPHDSLAM_2D{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  RBPHDSLAM_2D(){
    pFilter_ = NULL;
  }
  
  ~RBPHDSLAM_2D(){
    
    if(pFilter_ != NULL){
      delete pFilter_;
    }

  }

  /** Read the simulator configuration file */
  bool readConfigFile(const char* fileName){
    
    cfgFileName_ = fileName;

    boost::property_tree::ptree pt;
    boost::property_tree::xml_parser::read_xml(fileName, pt);

    logResultsToFile_ = false;
    if( pt.get("config.logging.logResultsToFile", 0) == 1 )
      logResultsToFile_ = true;
    logTimingToFile_ = false;
    if( pt.get("config.logging.logTimingToFile", 0) == 1 )
      logTimingToFile_ = true;
    logDirPrefix_ = pt.get<std::string>("config.logging.logDirPrefix", "./");
    if( *logDirPrefix_.rbegin() != '/')
      logDirPrefix_ += '/';

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
    if (std::is_same<MeasurementModel_2D, MeasurementModel_RngBrg>::value){
      varzr_ = pt.get<double>("config.measurements.varzr");
      varzb_ = pt.get<double>("config.measurements.varzb");
    }else{
      varzy_ = pt.get<double>("config.measurements.varzy");
      varzx_ = pt.get<double>("config.measurements.varzx");
    }


    
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
    
    return true;   
  }
  /** 
   *   This function reads isam style datasets. 
   * 
   * @param fileName The filename. 
   * 
   * @return True on  success false on failure.
   */
  bool readISAMFile(const char* fileName){
    std::ifstream isamFile(fileName);
    
    
    
    if( isamFile.fail() ){
      std::cerr << "Error while opening " << fileName << std::endl;
      exit(1);
    }
    TimeStamp t;
    MotionModel_Odometry2d::TState pose_k(t);
    deadReckoning_pose_.push_back(pose_k);
    std::string line;
    while(std::getline(isamFile, line )){
      
      parseISAMLine(line.c_str());
    }
    return true;
    
  }
  /** 
   *   This function reads an isam style dataset and uses it as ground truth. 
   * 
   * @param fileName The filename. 
   * 
   * @return True on  success false on failure.
   */
  bool readISAMGroundTruth(const char* fileName){
    std::ifstream isamFile(fileName);
    
    
    
    if( isamFile.fail() ){
      std::cerr << "Error while opening " << fileName << std::endl;
      exit(1);
    }
    TimeStamp t;
    MotionModel_Odometry2d::TState pose_k(t);
    groundtruth_pose_.push_back(pose_k);
    
    std::string line;
    while(std::getline(isamFile, line )){
      
      parseISAMLine(line.c_str(), true);
    }
    return true;
  }


  
  /** 
   *  Parse a single 2D isam style line
   * 
   * @param line 
   * 
   * @param t  
   *
   * @param isGroundTruth False if parsing a regular file, true if parsing the ground truth
   *
   * @return True on success
   */
  bool parseISAMLine(const char* line, bool isGroundTruth=false){
    char linestart[500];
    int linestartLength;
    sscanf(line, "%s%n", linestart, &linestartLength);
    std::string linestart_s(linestart);
    //----this is an odometry reading----------
    if( linestart_s == "ODOMETRY" || linestart_s == "EDGE2" || linestart_s == "Constraint2"){
      unsigned int i_x0,i_x1;
      double x, y, theta, info_xx, info_xy, info_xtheta, info_yy, info_ytheta, info_thetatheta;
      int nread = sscanf(line+linestartLength, "%i %i %lg %lg %lg %lg %lg %lg %lg %lg %lg", &i_x0, &i_x1, &x, &y, &theta, &info_xx, &info_xy, &info_xtheta, &info_yy, &info_ytheta, &info_thetatheta);
      if(nread!=11){
	std::cerr << "Error while reading odometry\n";
	 exit(1);
      }
      if(i_x1!=i_x0+1){
	std::cerr << "Constraint on nonconsecutive poses, it cannot be processed by filtering solutions.\n Ignore and continue?[Y/n]\n";
	std::string ans;
	std::cin >> ans;
	if( ans == "n" || ans == "N")
	  exit(1);
	else
	  std::cerr << "Continuing\n";
      }

      // get Covariance from information matrix

      MotionModel_Odometry2d::TInput::Mat information_sqrt, information, covariance;
      MotionModel_Odometry2d::TInput::Vec  odom_mean;

      odom_mean << x,y,theta;
      information_sqrt <<
	info_xx, info_xy, info_xtheta,
	0.,  info_yy, info_ytheta,
	0.,   0., info_thetatheta;
      information = information_sqrt.transpose()*information_sqrt;
      covariance = information.inverse();
      TimeStamp t(i_x1);
      MotionModel_Odometry2d::TInput odom(odom_mean, covariance, t);
      if (isGroundTruth){
	MotionModel_Odometry2d::TState::Mat Q;
	Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
	MotionModel_Odometry2d motionModel(Q);
    
	groundtruth_displacement_.push_back(odom);
	MotionModel_Odometry2d::TState x_k;
	motionModel.step(x_k, groundtruth_pose_.back(), odom, dTimeStamp_);
	groundtruth_pose_.push_back( x_k );
	groundtruth_pose_.back().setTime(odom.getTime());
	
      }
      else{
	// deadReckoning
	MotionModel_Odometry2d::TState::Mat Q;
	Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
	MotionModel_Odometry2d motionModel(Q);
    
        
	MotionModel_Odometry2d::TState x_k;
	motionModel.step(x_k, deadReckoning_pose_.back(), odom, dTimeStamp_);
	deadReckoning_pose_.push_back( x_k );
	deadReckoning_pose_.back().setTime(odom.getTime());
        
	vardx_= covariance(0,0);
	vardy_= covariance(1,1);
	vardz_= covariance(2,2);
	odometry_.push_back( odom);
      }
      
    }else if ( linestart_s == "LANDMARK"){
       unsigned int i_x,i_landmark;
       double x,y,info_xx,info_xy,info_yy;
       int nread  = sscanf(line+linestartLength, "%i %i %lg %lg %lg %lg %lg", &i_x, &i_landmark, &x, &y, &info_xx, &info_xy, &info_yy);
       if (nread != 7){
	 std::cerr << "Error while reading landmark\n";
	 exit(1);
       }
      MeasurementModel_RngBrg::TMeasurement::Mat information_sqrt, information, covariance;
      typename MeasurementModel_2D::TMeasurement::Vec  measurement_mean;

      measurement_mean << x,y;
      information_sqrt <<
	info_xx, info_xy,
	0.,  info_yy;
      information = information_sqrt.transpose()*information_sqrt;
      covariance = information.inverse();
      TimeStamp t(i_x);
      typename  MeasurementModel_2D::TMeasurement measurement(measurement_mean, covariance, t);
      if (isGroundTruth){
	MeasurementModel_2D measurementModel( varzx_, varzy_ );
	typename  MeasurementModel_2D::TLandmark lm;
	measurementModel.inverseMeasure( groundtruth_pose_[i_x], 
					 measurement, 
					 lm);
	groundtruth_landmark_.push_back(lm);
	lmkFirstObsTime_.push_back(i_x);
      }else {
	varzx_ = covariance(0,0);
	varzy_ = covariance(1,1);
	measurements_.push_back(measurement);
      }
    }
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

      if( k <= 50 ){
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
    MeasurementModel_2D measurementModel;
    typename MeasurementModel_2D::TMeasurement::Mat covZ;
    if (std::is_same<MeasurementModel_2D, MeasurementModel_RngBrg>::value)
      covZ << varzr_, 0, 0, varzb_;
    else
      covZ << varzx_, 0, 0, varzy_;
    measurementModel.setNoise(covZ);
    typename MeasurementModel_2D::TPose pose;

    groundtruth_landmark_.reserve(nLandmarks_);

    int nLandmarksCreated = 0;
    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nLandmarks_ * nLandmarksCreated){

	typename MeasurementModel_2D::TPose pose;
	typename MeasurementModel_2D::TMeasurement measurementToCreateLandmark;
	typename MeasurementModel_2D::TMeasurement::Vec z;
	double r = drand48() * rangeLimitMax_;
	double b = drand48() * 2 * PI;
	z << r, b;
	measurementToCreateLandmark.set(z);
	typename MeasurementModel_2D::TLandmark lm;
	
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
    MeasurementModel_2D measurementModel;
    typename MeasurementModel_2D::TMeasurement::Mat R;
    if (std::is_same<MeasurementModel_2D, MeasurementModel_RngBrg>::value)
               R << varzr_,0 ,0,  varzb_ ;
            else
              R << varzx_,0 ,0,  varzy_ ;

    measurementModel.setNoise(R);
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
	typename MeasurementModel_2D::TMeasurement z_m_k;
	success = measurementModel.sample( groundtruth_pose_[k],
					   groundtruth_landmark_[m],
					   z_m_k);
        if(success){
          double range;
          if (std::is_same<MeasurementModel_2D, MeasurementModel_RngBrg>::value){
            range = z_m_k.get(0);
          }else{
            range = sqrt(z_m_k.get(0)*z_m_k.get(0) + z_m_k.get(1)*z_m_k.get(1));
          }
          if(range <= rangeLimitMax_ && range >= rangeLimitMin_ && drand48() <= Pd_){
            z_m_k.setTime(t);
            //z_m_k.setCov(R);
            measurements_.push_back( z_m_k );
          }

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
	typename MeasurementModel_2D::TMeasurement z_clutter;
	typename MeasurementModel_2D::TMeasurement::Vec z;
	z << r, b;
	z_clutter.set(z, t);
	measurements_.push_back(z_clutter);
	
      }
      
    }
    
  }

  /** Data Logging */
  void exportSimData(){

    if(logResultsToFile_ || logTimingToFile_ ){
      boost::filesystem::path dir(logDirPrefix_);
      boost::filesystem::create_directories(dir);
      boost::filesystem::path cfgFilePathSrc( cfgFileName_ );
      std::string cfgFileDst( logDirPrefix_ );
      cfgFileDst += "simSettings.xml";
      boost::filesystem::path cfgFilePathDst( cfgFileDst.data() );
      boost::filesystem::copy_file( cfgFilePathSrc, cfgFilePathDst, boost::filesystem::copy_option::overwrite_if_exists);
    }

    if(!logResultsToFile_)
      return;

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
    typename MeasurementModel_2D::TLandmark::Vec m;
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
    typename MeasurementModel_2D::TMeasurement::Vec z;
    for(int i = 0; i < measurements_.size(); i++){
      measurements_[i].get(z, t);
      fprintf( pMeasurementFile, "%f   %f   %f\n", t.getTimeAsDouble(), sqrt(z(0)*z(0)+z(1)*z(1)), atan2(z(1),z(0)) );
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

  /** RB-PHD Filter Setup */
  void setupRBPHDFilter(){
    
    pFilter_ = new RBPHDFilter<MotionModel_Odometry2d,
			       StaticProcessModel<Landmark2d>,
			       MeasurementModel_2D,
			       KalmanFilter<StaticProcessModel<Landmark2d>,MeasurementModel_2D >>( nParticles_ );

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
    typename MeasurementModel_2D::TMeasurement::Mat R;

    if (std::is_same<MeasurementModel_2D, MeasurementModel_RngBrg>::value)
           R << varzr_,0 ,0,  varzb_ ;
        else
          R << varzx_,0 ,0,  varzy_ ;

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

#ifdef _PERFTOOLS_CPU
    std::string perfCPU_file = logDirPrefix_ + "rbphdslam2dSim_cpu.prof";
    ProfilerStart(perfCPU_file.data());
#endif
#ifdef _PERFTOOLS_HEAP
    std::string perfHEAP_file = logDirPrefix_ + "rbphdslam2dSim_heap.prof";
    HeapProfilerStart(perfHEAP_file.data());
#endif


    //////// Initialization at first timestep //////////

    if(!logResultsToFile_){
      std::cout << "Note: results are NOT being logged to file (see config xml file)\n";
    }
    FILE* pParticlePoseFile;
    if(logResultsToFile_){
      std::string filenameParticlePoseFile( logDirPrefix_ );
      filenameParticlePoseFile += "particlePose.dat";
      pParticlePoseFile = fopen(filenameParticlePoseFile.data(), "w");
    }
    FILE* pLandmarkEstFile;
    if(logResultsToFile_){
      std::string filenameLandmarkEstFile( logDirPrefix_ );
      filenameLandmarkEstFile += "landmarkEst.dat";
      pLandmarkEstFile = fopen(filenameLandmarkEstFile.data(), "w");
    }
    MotionModel_Odometry2d::TState x_i;
    int zIdx = 0;

    if(logResultsToFile_){
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	x_i = *(pFilter_->getParticleSet()->at(i));
	fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   1.0\n", 0.0, i, x_i.get(0), x_i.get(1), x_i.get(2));
      }   
    }
  
    /////////// Run simulator from k = 1 to kMax_ /////////

    TimeStamp time;

    

    for(int k = 0; k < kMax_; k++){

      time=odometry_[k].getTime();
      
      if( k % 100 == 0 || k == kMax_ - 1){
	float progressPercent = float(k+1) / float(kMax_);
	int progressBarW = 50;
	struct winsize ws;
	if(ioctl(1, TIOCGWINSZ, &ws) >= 0)
	  progressBarW = ws.ws_col - 30;
	int progressPos = progressPercent * progressBarW;
	if(progressBarW >= 50){
	  std::cout << "["; 
	  for(int i = 0; i < progressBarW; i++){
	    if(i < progressPos)
	      std::cout << "=";
	    else if(i == progressPos)
	      std::cout << ">";
	    else
	      std::cout << " ";
	  }
	  std::cout << "] ";
	}
	std::cout << "k = " << k << " (" << int(progressPercent * 100.0) << " %)\r"; 
	std::cout.flush();
      }
      if(k == kMax_ - 1)
	std::cout << std::endl << std::endl;

#ifdef _PERFTOOLS_HEAP
      if( k % 20 == 0)
	HeapProfilerDump("Timestep interval dump");
#endif
      
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
      
      if( k <= -1100){
	for( int i = 0; i < nParticles_; i++)
	  pFilter_->setParticlePose(i, groundtruth_pose_[k]);
      }

      // Prepare measurement vector for update
      std::vector<typename MeasurementModel_2D::TMeasurement> Z;
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
      int i_w_max = 0;
      double w_max = 0;
      if(logResultsToFile_){
	for(int i = 0; i < pFilter_->getParticleCount(); i++){
	  x_i = *(pFilter_->getParticleSet()->at(i));
	  double w = pFilter_->getParticleSet()->at(i)->getWeight();
	  if(w > w_max){
	    i_w_max = i;
	    w_max = w;
	  }
	  fprintf( pParticlePoseFile, "%f   %d   %f   %f   %f   %f\n", time.getTimeAsDouble(), i, x_i.get(0), x_i.get(1), x_i.get(2), w);
	}
	//fprintf( pParticlePoseFile, "\n");
      }

      // Log landmark estimates
      if(logResultsToFile_){

	int mapSize = pFilter_->getGMSize(i_w_max);
	for( int m = 0; m < mapSize; m++ ){
	  typename MeasurementModel_2D::TLandmark::Vec u;
	  typename MeasurementModel_2D::TLandmark::Mat S;
	  double w;
	  pFilter_->getLandmark(i_w_max, m, u, S, w);
	    
	  fprintf( pLandmarkEstFile, "%f   %d   ", time.getTimeAsDouble(), i_w_max);
	  fprintf( pLandmarkEstFile, "%f   %f      ", u(0), u(1));
	  fprintf( pLandmarkEstFile, "%f   %f   %f", S(0,0), S(0,1), S(1,1));
	  fprintf( pLandmarkEstFile, "   %f\n", w );
	  
	}
      }

    }

      
#ifdef _PERFTOOLS_HEAP
    HeapProfilerStop();
#endif
#ifdef _PERFTOOLS_CPU
    ProfilerStop();
#endif


    std::cout << "Elapsed Timing Information [nsec]\n";
    std::cout << std::setw(15) << std::left << "Prediction" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->predict_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->predict_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Update" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Update (KF)" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Weighting" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Merge" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Map Prune" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_wall
	      << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Resampling" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_wall
	      << std::setw(6) << std::left << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
    std::cout << std::setw(15) << std::left << "Total" << std::setw(15)
	      << std::setw(6) << std::right << "wall:" << std::setw(15)
	      << pFilter_->getTimingInfo()->predict_wall +
                 pFilter_->getTimingInfo()->mapUpdate_wall +
                 pFilter_->getTimingInfo()->particleWeighting_wall +
                 pFilter_->getTimingInfo()->mapMerge_wall +
                 pFilter_->getTimingInfo()->mapPrune_wall +
                 pFilter_->getTimingInfo()->particleResample_wall 
	      << std::setw(6) << std::right << "cpu:" << std::setw(15)
	      << pFilter_->getTimingInfo()->predict_cpu +
                 pFilter_->getTimingInfo()->mapUpdate_cpu +
                 pFilter_->getTimingInfo()->particleWeighting_cpu +
                 pFilter_->getTimingInfo()->mapMerge_cpu +
                 pFilter_->getTimingInfo()->mapPrune_cpu + 
                 pFilter_->getTimingInfo()->particleResample_cpu << std::endl;

    if(logTimingToFile_){
      std::ofstream timingFile( (logDirPrefix_ + "timing.dat").data() );
      timingFile << "Elapsed Timing Information [nsec]\n";
      timingFile << std::setw(15) << std::left << "Prediction" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->predict_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->predict_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Update" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Update (KF)" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Weighting" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleWeighting_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Merge" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapMerge_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Map Prune" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_wall
		 << std::setw(6) << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->mapPrune_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Resampling" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_wall
		 << std::setw(6) << std::left << std::right << "cpu:" << std::setw(15) << pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
      timingFile << std::setw(15) << std::left << "Total" << std::setw(15)
		 << std::setw(6) << std::right << "wall:" << std::setw(15)
		 << pFilter_->getTimingInfo()->predict_wall +
	pFilter_->getTimingInfo()->mapUpdate_wall +
	pFilter_->getTimingInfo()->particleWeighting_wall +
	pFilter_->getTimingInfo()->mapMerge_wall +
	pFilter_->getTimingInfo()->mapPrune_wall +
	pFilter_->getTimingInfo()->particleResample_wall 
		 << std::setw(6) << std::right << "cpu:" << std::setw(15)
		 << pFilter_->getTimingInfo()->predict_cpu +
	pFilter_->getTimingInfo()->mapUpdate_cpu +
	pFilter_->getTimingInfo()->particleWeighting_cpu +
	pFilter_->getTimingInfo()->mapMerge_cpu +
	pFilter_->getTimingInfo()->mapPrune_cpu + 
	pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
      timingFile.close();
    }
    
    if(logResultsToFile_){
      fclose(pParticlePoseFile);
      fclose(pLandmarkEstFile);
    }
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
  std::vector<typename MeasurementModel_2D::TLandmark> groundtruth_landmark_;
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
  double varzx_;
  double varzy_;
  std::vector<typename MeasurementModel_2D::TMeasurement> measurements_;

  // Filters
  KalmanFilter_RngBrg kf_;
  RBPHDFilter<MotionModel_Odometry2d, 
	      StaticProcessModel<Landmark2d>,
	      MeasurementModel_2D,  
	      KalmanFilter<StaticProcessModel<Landmark2d>,MeasurementModel_2D >> *pFilter_; 
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

  bool logResultsToFile_;
  bool logTimingToFile_;
  
public:
  std::string logDirPrefix_;
};









