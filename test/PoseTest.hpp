// Test for Pose classes
// Keith Leung 2013

#include "Pose.hpp" 

class PoseTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  PoseTest(){}
  
  /** Destructor */
  virtual ~PoseTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //

};

////////// Test Cases //////////

// Constructor test for Pose2d
TEST_F(PoseTest, constructorTestPose2d){
  
  Pose2d::Vec x, y;

  Pose2d p1;
  p1.get(x);
  EXPECT_EQ(Pose2d::Vec::Zero() ,x);

  y << 1.2, 3.4, 5.6;
  Pose2d p2(y);
  p2.get(x);
  EXPECT_EQ(y, x);

  Pose2d p3(1.2, 3.4, 5.6);
  p3.get(x);
  EXPECT_EQ(y, x);
 
}

// Constructor test for Pose1d
TEST_F(PoseTest, constructorTestPose1d){

  Eigen::Matrix<double, 1, 1> x, y;
  
  y << 0;

  Pose1d p1;
  p1.get(x);
  EXPECT_EQ(y, x);
  
  y << 1.2;
  Pose1d p2(y);
  p2.get(x);
  EXPECT_EQ(y, x);

}
