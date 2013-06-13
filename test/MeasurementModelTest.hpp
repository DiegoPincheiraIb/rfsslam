// Test for classes derived from MeasurementModel
// Keith Leung 2013

#include "gtest/gtest.h"
#include "MeasurementModel.hpp"

/**
 * \class MeasurementModelTest
 * \brief Unit testing fixture for MeasurementModel derived classes
 * \author Keith Leung 
 */

class MeasurementModelTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  MeasurementModelTest(){}
  
  /** Destructor */
  virtual ~MeasurementModelTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //

};

////////// Test Cases //////////


// Test rangeBeringModel constructor, and set / get functions
TEST_F(MeasurementModelTest, rangeBearingModelConstructorTest){
  
  RangeBearingModel model;
  Eigen::Matrix2d cov_in, cov_out1, cov_out2, cov_out3, cov_out4;
  double Sr, Sb;
  
  Sr = 1.1;
  Sb = 0.1;
  cov_in << Sr, 0, 0, Sb;
  
  model.setCov(Sr, Sb);
  model.getCov(cov_out1);
  EXPECT_EQ( cov_in, cov_out1 );
  
  model.getCov(cov_in);
  model.getCov(cov_out2);
  EXPECT_EQ( cov_in, cov_out2 );

  RangeBearingModel model2(Sr, Sb);
  model2.getCov(cov_out3);
  EXPECT_EQ( cov_in, cov_out3 );

  RangeBearingModel model3( cov_in );
  model3.getCov(cov_out4);
  EXPECT_EQ( cov_in, cov_out4 );
}

TEST_F(MeasurementModelTest, rangeBearingModelPredictTest){

  RangeBearingModel model(1, 1);
  RangeBearingModel::tPose x;
  Eigen::Vector3d xPose;
  RangeBearingModel::tLandmark m;
  Eigen::Vector2d mPos;
  Eigen::Matrix2d mCov;
  RangeBearingModel::tMeasurement z;
  Eigen::Vector2d zVec;
  Eigen::Matrix2d zCov;
  double t;
  
  xPose << 0, 0, 0;
  x.set(xPose);
  mPos << -1, -1;
  mCov << 1, 0, 0, 1;
  m.set(mPos, mCov);
  model.predict( x, m, z);
  z.get(zVec, zCov, t);
  EXPECT_EQ( sqrt(2), zVec(0) );
  EXPECT_EQ( -0.75 * PI , zVec(1) );
}

TEST_F(MeasurementModelTest, rangeBearingModelInvPredictTest){

  RangeBearingModel model(1, 1);
  RangeBearingModel::tPose x;
  Eigen::Vector3d xPose;
  RangeBearingModel::tLandmark m;
  Eigen::Vector2d mPos;
  Eigen::Matrix2d mCov;
  RangeBearingModel::tMeasurement z;
  Eigen::Vector2d zVec;
  Eigen::Matrix2d zCov;
  double t;

  xPose << 0, 0, 0;
  x.set(xPose);
  zVec << sqrt(2) , -0.75 * PI;
  zCov << 1, 0, 0, 1;
  z.set( zVec, zCov, 0.4565);
  model.inversePredict(x, z, m);
  m.get(mPos, mCov);
  EXPECT_DOUBLE_EQ(-1, mPos(0));
  EXPECT_DOUBLE_EQ(-1, mPos(1)); 

  // Is there a good way of checking mCov?

}
