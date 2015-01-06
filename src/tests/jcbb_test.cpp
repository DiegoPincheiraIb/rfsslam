#include "MeasurementModel_XY.hpp"
#include "JCBB.hpp"

using namespace rfs;

int main(int argc, char* args[]){

  std::cout << "JCBB Test\n";

  typedef MeasurementModel_XY::TPose TPose;
  typedef MeasurementModel_XY::TLandmark TLandmark;
  typedef MeasurementModel_XY::TMeasurement TZ;

  MeasurementModel_XY z_model(1,1);

  Eigen::MatrixXd cov;
  unsigned int cov_dim = TPose::Vec::RowsAtCompileTime + TLandmark::Vec::RowsAtCompileTime * 5;
  cov.setZero(cov_dim, cov_dim);

  // Robot Pose
  TPose::Vec p;
  TPose::Vec pCovDiag;
  p << 0, 0, 0;
  pCovDiag << 1, 1, 0;
  TPose robot(p, pCovDiag.asDiagonal(), 0);
  cov.topLeftCorner(TPose::Vec::RowsAtCompileTime, TPose::Vec::RowsAtCompileTime) = robot.getCov();

  // Landmark positions // todo change to pointers
  std::vector<TLandmark> landmarks;
  std::vector<TLandmark*> lmkPtrs;
  TLandmark::Vec m;
  TLandmark::Mat P;
  P << 0.1, 0, 0, 0.1;
  m << -10, 0;
  landmarks.push_back(TLandmark(m, P));
  m << -8, 0;
  landmarks.push_back(TLandmark(m, P));
  m << 0, 10;
  landmarks.push_back(TLandmark(m, P));
  m << 10, -1;
  landmarks.push_back(TLandmark(m, P));
  m << 10, 1;
  landmarks.push_back(TLandmark(m, P));
  for(int i = 0; i < 5; i++){
    lmkPtrs.push_back( &(landmarks[i]) );
    int pos = TPose::Vec::RowsAtCompileTime + i * TLandmark::Vec::RowsAtCompileTime;
    cov.block(pos, pos, TLandmark::Vec::RowsAtCompileTime, TLandmark::Vec::RowsAtCompileTime) = P;
  }
  // std::cout << "Covariance:\n" << cov << std::endl;
  
  // Measurements
  std::vector<TZ> Z(6);
  TZ::Vec z;
  z << 0, -15;
  Z[0].set(z);

  z << 10, 0.9;
  Z[1].set(z);

  z << -9, 0;
  Z[2].set(z);

  z << -10.1, 0;
  Z[3].set(z);

  z << 0, -10;
  Z[4].set(z);

  z << 10, 0;
  Z[5].set(z);

  JCBB<MeasurementModel_XY> jcbb(0.95, &z_model, &Z, &robot, &lmkPtrs);

  for(int n = 0; n < 6; n++){
    std::cout << n << " --- " << jcbb.getAssociation(n) << "\n";
  } 

  JCBB<MeasurementModel_XY> jcbb2(0.95, &z_model, &Z, &robot, &lmkPtrs, &cov );

  for(int n = 0; n < 6; n++){
    std::cout << n << " --- " << jcbb2.getAssociation(n) << "\n";
  } 

} 
