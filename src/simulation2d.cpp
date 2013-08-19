// Simulator for testing the RBPHDFilter
// Keith Leung 2013
// Requires libconfig to be installed

#include <boost/lexical_cast.hpp>
#include <libconfig.h++>
#include "RBPHDFilter.hpp"
#include <stdio.h>

class Simulator2d{

public:

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  Simulator2d(){
    
    pFilter_ = NULL;
   
  }
  
  ~Simulator2d(){
    
    if(pFilter_ != NULL){
      delete pFilter_;
    }

  }

  /** Read the simulator configuration file */
  bool readConfigFile(){
    
    libconfig::Config cfg;
    const char* fileName= "cfg/simulator2d.cfg";
    try{
      cfg.readFile( fileName );
    }catch( libconfig::FileIOException &ex){
      printf("\nCannot read file: %s\n\n", fileName);
      return false;
    }catch( libconfig::ParseException &ex){
      const char* error = ex.getError();
      int line = ex.getLine();
      printf("\n%s LINE %d\n\n", error, line);
      return false;
    }
    kMax_ = cfg.lookup("timesteps");
    dT_ = cfg.lookup("sec_per_timestep");
    
    nSegments_ = cfg.lookup("Trajectory.nSegments");
    max_dx_ = cfg.lookup("Trajectory.max_dx_per_sec");
    max_dy_ = cfg.lookup("Trajectory.max_dy_per_sec");
    max_dz_ = cfg.lookup("Trajectory.max_dz_per_sec");
    min_dx_ = cfg.lookup("Trajectory.min_dx_per_sec");
    vardx_ = cfg.lookup("Trajectory.vardx");
    vardy_ = cfg.lookup("Trajectory.vardy");
    vardz_ = cfg.lookup("Trajectory.vardz");

    nLandmarks_ = cfg.lookup("Landmarks.nLandmarks");
    varlmx_ = cfg.lookup("Landmarks.varlmx");
    varlmy_ = cfg.lookup("Landmarks.varlmy");

    rangeLimitMax_ = cfg.lookup("Measurement.rangeLimitMax");
    rangeLimitMin_ = cfg.lookup("Measurement.rangeLimitMin");
    rangeLimitBuffer_ = cfg.lookup("Measurement.rangeLimitBuffer");
    Pd_ = cfg.lookup("Measurement.probDetection");
    c_ = cfg.lookup("Measurement.clutterIntensity");
    varzr_ = cfg.lookup("Measurement.varzr");
    varzb_ = cfg.lookup("Measurement.varzb");

    nParticles_ = cfg.lookup("Filter.nParticles");
    zNoiseInflation_ = cfg.lookup("Filter.measurementNoiseInflationFactor");
    birthGaussianWeight_ = cfg.lookup("Filter.birthGaussianWeight");
    newGaussianCreateInnovMDThreshold_ = cfg.lookup("Filter.newGaussianCreateInnovMDThreshold");
    importanceWeightingMeasurementLikelihoodMDThreshold_ = cfg.lookup("Filter.importanceWeightingMeasurementLikelihoodMDThreshold");
    effNParticleThreshold_ = cfg.lookup("Filter.effectiveNumberOfParticlesThreshold");
    minInterSampleTimesteps_ = cfg.lookup("Filter.minInterSampleTimesteps");
    gaussianMergingThreshold_ = cfg.lookup("Filter.gaussianMergingThreshold");
    gaussianMergingCovarianceInflationFactor_ = cfg.lookup("Filter.gaussianMergingCovarianceInflationFactor");
    gaussianPruningThreshold_ = cfg.lookup("Filter.gaussianPruningThreshold");
    return true;   
  }

  /** Generate a random trajectory in 2d space */
  void generateTrajectory(int randSeed = 0){

    printf("Generating trajectory with random seed = %d\n", randSeed);
    srand48( randSeed);

    int seg = 0;
    OdometryMotionModel2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    Q = dT_ * Q * dT_;
    // std::cout << "\n\n" << Q << "\n\n";
    OdometryMotionModel2d motionModel(Q);
    OdometryMotionModel2d::TInput input_k(0, 0, 0, 0, 0, 0, 0);
    OdometryMotionModel2d::TState pose_k(0, 0, 0, 0, 0, 0, 0);
    OdometryMotionModel2d::TState pose_km(0, 0, 0, 0, 0, 0, 0);
    groundtruth_displacement_.reserve( kMax_ );
    groundtruth_pose_.reserve( kMax_ );
    groundtruth_displacement_.push_back(input_k);
    groundtruth_pose_.push_back(pose_k);

    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nSegments_ * seg ){
	seg++;
	double dx = drand48() * max_dx_ * dT_;
	while( dx < min_dx_ * dT_ ){
	  dx = drand48() * max_dx_ * dT_;
	}
	double dy = (drand48() * max_dy_ * 2 - max_dy_) * dT_;
	double dz = (drand48() * max_dz_ * 2 - max_dz_) * dT_; 
	input_k = OdometryMotionModel2d::TInput(dx, dy, dz, 
						Q(0,0), Q(1,1), Q(2,2), k);  
      }

      groundtruth_displacement_.push_back(input_k);
      groundtruth_displacement_.back().setTime(k);

      OdometryMotionModel2d::TState x_k;
      motionModel.step(x_k, groundtruth_pose_[k-1], input_k);
      groundtruth_pose_.push_back( x_k );
      groundtruth_pose_.back().setTime(k);
      
      /*OdometryMotionModel2d::TState::Vec x;
      OdometryMotionModel2d::TInput::Vec u;
      double t;
      groundtruth_pose_[k].get(x);
      groundtruth_displacement_[k].get(u, t);
      printf("x[%d] = [%f %f %f]  u[%f] = [%f %f %f]\n", 
           k, x(0), x(1), x(2), t, u(0), u(1), u(2)); */

    }

  }
  
  void generateOdometry(){

    odometry_.reserve( kMax_ );
    OdometryMotionModel2d::TInput zero;
    OdometryMotionModel2d::TInput::Vec u0;
    u0.setZero();
    zero.set(u0, 0);
    odometry_.push_back( zero );


    OdometryMotionModel2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    OdometryMotionModel2d motionModel(Q);
    deadReckoning_pose_.reserve( kMax_ );
    deadReckoning_pose_.push_back( groundtruth_pose_[0] );

    for( int k = 1; k < kMax_; k++){
      
      OdometryMotionModel2d::TInput in = groundtruth_displacement_[k];
      in.setCov(Q);
      OdometryMotionModel2d::TInput out;
      RandomVecMathTools<OdometryMotionModel2d::TInput>::sample(in, out);
      /*printf("u[%d] = [%f %f %f]  u_[%d] = [%f %f %f]\n", 
	     k, in.get(0),  in.get(1),  in.get(2), 
	     k, out.get(0), out.get(1), out.get(2) );*/
      
      odometry_.push_back( out );

      OdometryMotionModel2d::TState p;
      motionModel.step(p, deadReckoning_pose_[k-1], odometry_[k]);
      p.setTime(k);
      deadReckoning_pose_.push_back( p );
    }

  }

  void generateLandmarks(){

    RangeBearingModel measurementModel( varzr_, varzb_);
    RangeBearingModel::TPose pose;

    groundtruth_landmark_.reserve(nLandmarks_);

    int nLandmarksCreated = 0;
    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nLandmarks_ * nLandmarksCreated){

	RangeBearingModel::TPose pose;
	RangeBearingModel::TMeasurement measurementToCreateLandmark;
	RangeBearingModel::TMeasurement::Vec z;
	double r = drand48() * rangeLimitMax_;
	double b = drand48() * 2 * PI;
	z << r, b;
	measurementToCreateLandmark.set(z);
	RangeBearingModel::TLandmark lm;
	
	measurementModel.inverseMeasure( groundtruth_pose_[k], 
					 measurementToCreateLandmark, 
					 lm);

	groundtruth_landmark_.push_back(lm);

	nLandmarksCreated++;
	
      }

    }

  }

  void generateMeasurements(){


    RangeBearingModel measurementModel( varzr_, varzb_ );
    RangeBearingModel::TMeasurement::Mat R;
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

    for( int k = 1; k < kMax_; k++ ){
      
      groundtruth_pose_[k];
      
      // Real detections
      for( int m = 0; m < groundtruth_landmark_.size(); m++){
	
	bool success;
	RangeBearingModel::TMeasurement z_m_k;
	success = measurementModel.sample( groundtruth_pose_[k],
					   groundtruth_landmark_[m],
					   z_m_k);
	/*success = measurementModel.measure( groundtruth_pose_[k],
					   groundtruth_landmark_[m],
					   z_m_k);*/
	if(success){
   
	  if(z_m_k.get(0) <= rangeLimitMax_ && z_m_k.get(0) >= rangeLimitMin_ && drand48() <= Pd_){
	    /*printf("Measurement[%d] = [%f %f]\n", int(measurements_.size()),
	      z_m_k.get(0), z_m_k.get(1)); */
	    z_m_k.setTime(k);
	    z_m_k.setCov(R);
	    measurements_.push_back( z_m_k );
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
	double b = drand48() * 2 * PI - PI;
	RangeBearingModel::TMeasurement z_clutter;
	RangeBearingModel::TMeasurement::Vec z;
	z << r, b;
	z_clutter.set(z, k);
	measurements_.push_back(z_clutter);
	
      }
      
    }
    
  }

  void exportSimData(){

    double t;

    FILE* pGTPoseFile;
    pGTPoseFile = fopen("data/gtPose.dat", "w");
    OdometryMotionModel2d::TState::Vec x;
    for(int i = 0; i < groundtruth_pose_.size(); i++){
      groundtruth_pose_[i].get(x, t);
      fprintf( pGTPoseFile, "%f   %f   %f   %f\n", t, x(0), x(1), x(2));
    }
    fclose(pGTPoseFile);

    FILE* pGTLandmarkFile;
    pGTLandmarkFile = fopen("data/gtLandmark.dat", "w");
    RangeBearingModel::TLandmark::Vec m;
    for(int i = 0; i < groundtruth_landmark_.size(); i++){
      groundtruth_landmark_[i].get(m);
      fprintf( pGTLandmarkFile, "%f   %f\n", m(0), m(1));
    }
    fclose(pGTLandmarkFile);

    FILE* pOdomFile;
    pOdomFile = fopen("data/odometry.dat","w");
    OdometryMotionModel2d::TInput::Vec u;
    for(int i = 0; i < odometry_.size(); i++){
      odometry_[i].get(u, t);
      fprintf( pOdomFile, "%f   %f   %f   %f\n", t, u(0), u(1), u(2));
    }
    fclose(pOdomFile);

    FILE* pMeasurementFile;
    pMeasurementFile = fopen("data/measurement.dat", "w");
    RangeBearingModel::TMeasurement::Vec z;
    for(int i = 0; i < measurements_.size(); i++){
      measurements_[i].get(z, t);
      fprintf( pMeasurementFile, "%f   %f   %f\n", t, z(0), z(1) );
    }
    fclose(pMeasurementFile);

    FILE* pDeadReckoningFile;
    pDeadReckoningFile = fopen("data/deadReckoning.dat", "w");
    OdometryMotionModel2d::TState::Vec odo;
    for(int i = 0; i < deadReckoning_pose_.size(); i++){
      deadReckoning_pose_[i].get(odo, t);
      fprintf( pDeadReckoningFile, "%f   %f   %f   %f\n", t, odo(0), odo(1), odo(2));
    }
    fclose(pDeadReckoningFile);

  }

  void setupRBPHDFilter(){
    
    pFilter_ = new RBPHDFilter<OdometryMotionModel2d,
			       StaticProcessModel<Landmark2d>,
			       RangeBearingModel,
			       RangeBearingKalmanFilter>( nParticles_ );

    // configure robot motion model
    OdometryMotionModel2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    pFilter_->getProcessModel()->setNoise(Q);

    // configure landmark process model
    Landmark2d::Mat Q_lm;
    Q_lm << varlmx_, 0, 0, varlmy_;
    pFilter_->getLmkProcessModel()->setNoise(Q_lm);

    // configure measurement model
    RangeBearingModel::TMeasurement::Mat R;
    R << varzr_, 0, 0, varzb_;
    R *= zNoiseInflation_;
    pFilter_->getMeasurementModel()->setNoise(R);
    pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
    pFilter_->getMeasurementModel()->config.uniformClutterIntensity_ = c_;
    pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
    pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
    pFilter_->getMeasurementModel()->config.rangeLimBuffer_ = rangeLimitBuffer_;

    // configure the filter
    pFilter_->config.birthGaussianWeight_ = birthGaussianWeight_;
    pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
    pFilter_->config.minInterSampleTimesteps_ = minInterSampleTimesteps_;
    pFilter_->config.newGaussianCreateInnovMDThreshold_ = newGaussianCreateInnovMDThreshold_;
    pFilter_->config.importanceWeightingMeasurementLikelihoodMDThreshold_ = importanceWeightingMeasurementLikelihoodMDThreshold_;
    pFilter_->config.importanceWeightingEvalPointCount_ = importanceWeightingEvalPointCount_;
    pFilter_->config.gaussianMergingThreshold_ = gaussianMergingThreshold_;
    pFilter_->config.gaussianMergingCovarianceInflationFactor_ = gaussianMergingCovarianceInflationFactor_;
    pFilter_->config.gaussianPruningThreshold_ = gaussianPruningThreshold_;
  }

  void run(){

    printf("Running simulation\n\n");
    
    FILE* pParticlePoseFile;
    pParticlePoseFile = fopen("data/particlePose.dat", "w");
    fprintf( pParticlePoseFile, "Timesteps: %d\n", kMax_);
    fprintf( pParticlePoseFile, "nParticles: %d\n\n", pFilter_->getParticleCount());

    FILE* pLandmarkEstFile;
    pLandmarkEstFile = fopen("data/landmarkEst.dat", "w");
    fprintf( pLandmarkEstFile, "Timesteps: %d\n", kMax_);
    fprintf( pLandmarkEstFile, "nParticles: %d\n\n", pFilter_->getParticleCount());

    OdometryMotionModel2d::TState x_i;
    int zIdx = 0;

    for(int i = 0; i < pFilter_->getParticleCount(); i++){
      pFilter_->getParticleSet()->at(i)->getPose(x_i);
      fprintf( pParticlePoseFile, "%f   %f   %f   1.0\n", x_i.get(0), x_i.get(1), x_i.get(2));
    }

    for(int k = 1; k < kMax_; k++){

      if( k % 1 == 0)
	printf("k = %d\n", k);
      fprintf( pParticlePoseFile, "k = %d\n", k);

      pFilter_->predict( odometry_[k], k );
      //for( int i = 0; i < nParticles_; i++)
      //pFilter_->setParticlePose(0, groundtruth_pose_[k]);

      // Prepare measurement vector for update
      std::vector<RangeBearingModel::TMeasurement> Z;
      double kz = measurements_[ zIdx ].getTime();
      while( kz == k ){
	Z.push_back( measurements_[zIdx] );
	zIdx++;
	if(zIdx >= measurements_.size())
	  break;
	kz = measurements_[ zIdx ].getTime();
      }

      pFilter_->update(Z, k);
      
      // Log particle poses
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	pFilter_->getParticleSet()->at(i)->getPose(x_i);
	double w = pFilter_->getParticleSet()->at(i)->getWeight();
	fprintf( pParticlePoseFile, "%f   %f   %f   %f\n", x_i.get(0), x_i.get(1), x_i.get(2), w);
      }
      fprintf( pParticlePoseFile, "\n");

      // Log landmark estimates
      for(int i = 0; i < pFilter_->getParticleCount(); i++){
	int mapSize = pFilter_->getGMSize(i);
	fprintf( pLandmarkEstFile, "Timestep: %d\tParticle: %d\tMap Size: %d\n", k, i, mapSize);
	for( int m = 0; m < mapSize; m++ ){
	  RangeBearingModel::TLandmark::Vec u;
	  RangeBearingModel::TLandmark::Mat S;
	  double w;
	  pFilter_->getLandmark(i, m, u, S, w);
	  
	  fprintf( pLandmarkEstFile, "%f   %f      ", u(0), u(1));
	  fprintf( pLandmarkEstFile, "%f   %f   ", S(0,0), S(0,1));
	  fprintf( pLandmarkEstFile, "%f   %f   ", S(1,0), S(1,1));
	  fprintf( pLandmarkEstFile, "   %f\n", w );
	}
	fprintf( pLandmarkEstFile, "\n");
      }
      fprintf( pLandmarkEstFile, "\n");


    }
    fclose(pParticlePoseFile);
    fclose(pLandmarkEstFile);

  }

private:

  int kMax_; /**< number of timesteps */
  double dT_;

  // Trajectory
  int nSegments_;
  double max_dx_;
  double max_dy_;
  double max_dz_;
  double min_dx_;
  double vardx_;
  double vardy_;
  double vardz_;
  std::vector<OdometryMotionModel2d::TInput> groundtruth_displacement_;
  std::vector<OdometryMotionModel2d::TState> groundtruth_pose_;
  std::vector<OdometryMotionModel2d::TInput> odometry_;
  std::vector<OdometryMotionModel2d::TState> deadReckoning_pose_;

  // Landmarks 
  int nLandmarks_;
  std::vector<RangeBearingModel::TLandmark> groundtruth_landmark_;
  double varlmx_;
  double varlmy_;

  // Range-Bearing Measurements
  double rangeLimitMax_;
  double rangeLimitMin_;
  double rangeLimitBuffer_;
  double Pd_;
  double c_;
  double varzr_;
  double varzb_;
  std::vector<RangeBearingModel::TMeasurement> measurements_;

  // Filters
  RangeBearingKalmanFilter kf_;
  RBPHDFilter<OdometryMotionModel2d, 
	      StaticProcessModel<Landmark2d>,
	      RangeBearingModel,  
	      RangeBearingKalmanFilter> *pFilter_; 
  int nParticles_;
  double zNoiseInflation_;
  double birthGaussianWeight_;
  double newGaussianCreateInnovMDThreshold_;
  double importanceWeightingMeasurementLikelihoodMDThreshold_;
  double effNParticleThreshold_;
  int minInterSampleTimesteps_;
  double gaussianMergingThreshold_;
  double gaussianMergingCovarianceInflationFactor_;
  double gaussianPruningThreshold_;
  int importanceWeightingEvalPointCount_;
};

int main(int argc, char* argv[]){

  int initRandSeed = 0;
  if( argc == 2 ){
    initRandSeed = boost::lexical_cast<int>(argv[1]);
  }

  OdometryMotionModel2d g;
  RangeBearingModel h;

  Simulator2d sim;
  if( !sim.readConfigFile() ){
    return -1;
  }
  sim.generateTrajectory( initRandSeed );
  sim.generateOdometry();
  sim.generateLandmarks();
  sim.generateMeasurements();
  sim.exportSimData();
 
  sim.setupRBPHDFilter();
  sim.run();

  return 0;

}
