// Simulator for testing the RBPHDFilter
// Keith Leung 2013
// Requires libconfig to be installed

#include <libconfig.h++>
#include "RBPHDFilter.hpp"

class Simulator2d{

public:

  Simulator2d(){}
  ~Simulator2d(){}

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
    vardx_ = cfg.lookup("Trajectory.vardx");
    vardy_ = cfg.lookup("Trajectory.vardy");
    vardz_ = cfg.lookup("Trajectory.vardz");
    
    printf("kMax_ %d\n", kMax_);
    printf("dT_ %f\n", dT_);
    printf("nSegments_ %d\n", nSegments_);
    printf("max_dx_ %f\n", max_dx_);
    printf("max_dy_ %f\n", max_dy_);
    printf("max_dz_ %f\n", max_dz_);
    printf("vardx_ %f\n", vardx_);
    printf("vardy_ %f\n", vardy_);
    printf("vardz_ %f\n", vardz_);

    return true;
  }

  /** Generate a random trajectory in 2d space */
  void generateTrajectory(){

    int seg = 0;
    OdometryMotionModel2d::TState::Mat Q;
    Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
    OdometryMotionModel2d motionModel(Q);
    OdometryMotionModel2d::TInput input_k(0, 0, 0, 0, 0, 0, 0);
    OdometryMotionModel2d::TState pose_k(0, 0, 0, 0, 0, 0, 0);
    OdometryMotionModel2d::TState pose_km(0, 0, 0, 0, 0, 0, 0);
    groundtruth_displacement_.resize( kMax_ );
    groundtruth_pose_.resize( kMax_ );
    groundtruth_displacement_[0] = input_k;

    for( int k = 1; k < kMax_; k++ ){

      if( k >= kMax_ / nSegments_ * seg ){
	seg++;
	double dx = drand48() * max_dx_ * dT_;
	double dy = (drand48() * max_dy_ * 2 - max_dy_) * dT_;
	double dz = (drand48() * max_dz_ * 2 - max_dz_) * dT_; 
	input_k = OdometryMotionModel2d::TInput(dx, dy, dz, 0, 0, 0, 1);  
      }

      groundtruth_displacement_[k] = input_k;
      groundtruth_displacement_[k].set(k);

      motionModel.step(groundtruth_pose_[k],
		       groundtruth_pose_[k-1], 
		       input_k);
      groundtruth_pose_[k].set(k);
      
      OdometryMotionModel2d::TState::Vec x;
      OdometryMotionModel2d::TInput::Vec u;
      double t;
      groundtruth_pose_[k].get(x);
      groundtruth_displacement_[k].get(u, t);
      printf("x[%d] = [%f %f %f]  u[%f] = [%f %f %f]\n", 
	     k, x(0), x(1), x(2), t, u(0), u(1), u(2));

    }

  }
  
  void generateOdometry(){

  }

  void generateLandmarks(){

  }

  void generateMeasurements(){

  }

private:

  int kMax_; /**< number of timesteps */
  double dT_;

  // Trajectory
  int nSegments_;
  double max_dx_;
  double max_dy_;
  double max_dz_;
  double vardx_;
  double vardy_;
  double vardz_;
  std::vector<OdometryMotionModel2d::TInput> groundtruth_displacement_;
  std::vector<OdometryMotionModel2d::TState> groundtruth_pose_;

};

int main(int argc, char* argv[]){

  OdometryMotionModel2d g;
  RangeBearingModel h;

  Simulator2d sim;
  if( !sim.readConfigFile() ){
    return -1;
  }
  sim.generateTrajectory();

  return 0;

}
