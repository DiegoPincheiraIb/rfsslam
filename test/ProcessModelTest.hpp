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
