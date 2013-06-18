// Test for classes derived from ProcessModel
// Keith Leung 2013

#include "gtest/gtest.h"
#include "ProcessModel.hpp"

/**
 * \class ProcessModelTest
 * \brief Unit testing fixture for ProcessModel derived classes
 * \author Keith Leung 
 */

class ProcessModelTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  ProcessModelTest(){}
  
  /** Destructor */
  virtual ~ProcessModelTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //
  OdometryMotionModel2d model;
  Pose2d p_in;
  Pose2d p_out;
  Pose2d::Vec x_in;
  Pose2d::Vec x_out;
  Pose2d::Vec x_expect;
  Odometry2d odo;
  Odometry2d::MeasureType u;
  Odometry2d::MeasureType u2;
  Odometry2d::MeasureUncertaintyType Su;
  double t;
};

////////// Test Cases //////////

TEST_F(ProcessModelTest, OdometryMotionModel2dTest){
  
  
  x_in << 1, 2, 0;
  p_in.set(x_in);
  
  u << 0, 0, 0;
  Su << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  odo.set(u, Su, 0.3456);
  
  model.step(p_out, p_in, odo);
  
  p_out.get(x_out);
  EXPECT_NEAR(0, (x_out - x_in).norm(), 1e-15) << "2d odometry model with 0 input failed";
  

  x_in << 2, 2, 0; 
  p_in.set(x_in);
  
  u << 3, 3, 0.678;
  Su << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  odo.set(u, Su, 0.3456);
  
  model.step(p_out, p_in, odo);
  
  p_out.get(x_out);
  x_expect = x_in + u;

  EXPECT_NEAR(0, (x_out - x_expect).norm(), 1e-15) << "2d odometry model with translation input failed";



  x_in << 2, 2, 0.5; 
  p_in.set(x_in);
  
  u << 0, 0, 0.678;
  Su << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  odo.set(u, Su, 0.3456);
  
  model.step(p_out, p_in, odo);
  
  p_out.get(x_out);
  x_expect = x_in + u;

  EXPECT_NEAR(0, (x_out - x_expect).norm(), 1e-15) << "2d odometry model with rotation input failed";

}

TEST_F(ProcessModelTest, sampleTest){

  Pose2d::Mat ProcessNoise;
  ProcessNoise << 3, 2, 1, 2, 4, -1, 1, -1, 5;

  model = OdometryMotionModel2d( ProcessNoise );
  
  u << 1, 2, 0;
  Su << 0, 0, 0, 0, 0, 0, 0, 0, 0;
  odo.set(u, Su, 0.3456);

  x_in << 0, 0, 0;
  p_in.set(x_in);

  int nSamples = 200000;
  std::vector < Pose2d::Vec > v;
  for(int i = 0; i < nSamples; i++){
    model.sample(p_out, p_in, odo);  
    p_out.get(x_out);
    v.push_back(x_out);
  }

  Pose2d::Vec sum;
  sum = Pose2d::Vec::Zero();
  for(int i = 0; i < nSamples; i++){
    sum = sum + v[i];
  }  

  Pose2d::Vec mean = sum / nSamples;
  Pose2d::Vec e = mean - (x_in + u);
  EXPECT_NEAR(0.0, e.norm(), 1e-2) << "Got:\n" << mean << "\nExpected:\n" << (x_in + u) << "\n";


  Pose2d::Mat sum2;
  sum2 = Pose2d::Mat::Zero();
  for(int i = 0; i < nSamples; i++){
    e = v[i] - mean;
    sum2 = sum2 + e * e.transpose();
  }  
  Pose2d::Mat cov = sum2 / nSamples;  

  for( int i = 0; i < mean.size(); i++){
    for( int j = 0; j < mean.size(); j++){
      double e = cov(i,j) - ProcessNoise(i,j);
      EXPECT_NEAR(0.0, e, 1e-2) << "Failure may be due to randomness\nExpected:" 
				<< ProcessNoise(i,j)
				<< "\nGot:" << cov(i,j) << "\n\n";
    }
  }

}
