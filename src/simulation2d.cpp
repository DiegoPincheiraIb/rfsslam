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

    rangeLimit_ = cfg.lookup("Measurement.rangeLimit");
    Pd_ = cfg.lookup("Measurement.probDetection");
    Pfa_ = cfg.lookup("Measurement.probFalseAlarm");
    varzr_ = cfg.lookup("Measurement.varzr");
    varzb_ = cfg.lookup("Measurement.varzb");

    nParticles_ = cfg.lookup("Filter.nParticles");

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
      groundtruth_displacement_.back().set(k);

      OdometryMotionModel2d::TState x_k;
      motionModel.step(x_k, groundtruth_pose_[k-1], input_k);
      groundtruth_pose_.push_back( x_k );
      groundtruth_pose_.back().set(k);
      
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

    for( int k = 1; k < kMax_; k++){
      
      OdometryMotionModel2d::TInput &in = groundtruth_displacement_[k];
      OdometryMotionModel2d::TInput out;
      RandomVecMathTools<OdometryMotionModel2d::TInput>::sample(in, out);
      /*printf("u[%d] = [%f %f %f]  u_[%d] = [%f %f %f]\n", 
	     k, in.get(0),  in.get(1),  in.get(2), 
	     k, out.get(0), out.get(1), out.get(2) );*/
      
      odometry_.push_back( out );

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
	double r = drand48() * rangeLimit_;
	double b = drand48() * 2 * PI;
	z << r, b;
	measurementToCreateLandmark.set(z);
	RangeBearingModel::TLandmark lm;
	
	measurementModel.inverseMeasure( groundtruth_pose_[k], 
					 measurementToCreateLandmark, 
					 lm);

	groundtruth_landmark_.push_back(lm);
	//groundtruth_landmark_[nLandmarksCreated] = lm;
	
	/*printf("Landmark[%d] = [%f %f]\n", nLandmarksCreated,
	       groundtruth_landmark_[nLandmarksCreated].State::get(0),
	       groundtruth_landmark_[nLandmarksCreated].State::get(1));*/

	nLandmarksCreated++;
	
      }

    }

  }

  void generateMeasurements(){

    std::vector<RangeBearingModel::TMeasurement> z;
    for(int m = 0; m < 100000; m++){
      RangeBearingModel::TMeasurement a;
      z.push_back(a);
    }

    RangeBearingModel measurementModel( varzr_, varzb_);
    
    for( int k = 1; k < kMax_; k++ ){
      
      groundtruth_pose_[k];
      
      for( int m = 0; m < groundtruth_landmark_.size(); m++){
	
	RangeBearingModel::TMeasurement z_m_k;
	measurementModel.sample( groundtruth_pose_[k],
				 groundtruth_landmark_[m],
				 z_m_k);
	
	z_m_k.set(k);
	
	/*printf("Measurement[%d] = [%f %f]\n", int(measurements_.size()),
	  z_m_k.get(0), z_m_k.get(1)); */
	
        if(z_m_k.get(0) <= rangeLimit_)
	  measurements_.push_back( z_m_k );

      }
      
    }
    
  }

  void exportSimData(){

    double t;

    FILE* pGTPoseFile;
    pGTPoseFile = fopen("data/gtPose.dat", "w");
    OdometryMotionModel2d::TInput::Vec x;
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
    // Landmark2d::Mat Q_lm;
    // Q_lm.setZero();
    // pFilter_->getLmkProcessModel()->setNoise(Q_lm);

    // configure measurement model
    RangeBearingModel::TMeasurement::Mat R;
    R << varzr_, 0, 0, varzb_;
    pFilter_->getMeasurementModel()->setNoise(R);
    
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

  // Landmarks 
  int nLandmarks_;
  std::vector<RangeBearingModel::TLandmark> groundtruth_landmark_;

  // Range-Bearing Measurements
  double rangeLimit_;
  double Pd_;
  double Pfa_;
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

  return 0;

}
