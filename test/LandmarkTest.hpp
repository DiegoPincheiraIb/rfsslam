// Test for classes derived from Landmark
// Keith Leung 2013

#include "gtest/gtest.h"
#include "Landmark.hpp"

/**
 * \class LandmarkTest
 * \brief Unit testing fixture for Landmark derived classes
 * \author Keith Leung 
 */

class LandmarkTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  LandmarkTest(){}
  
  /** Destructor */
  virtual ~LandmarkTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //

};

////////// Test Cases //////////

// Constructor test for Landmark2d
TEST_F(LandmarkTest, constructorTestLandmark2d){

  Landmark2d::Vec x_in, x_out;
  Landmark2d::Mat Sx_in, Sx_out;
  x_in << 0, 0;
  Sx_in << 0, 0, 0, 0;
  Landmark2d m;
  m.get(x_out);
  EXPECT_EQ(x_in, x_out);
  m.get(x_out, Sx_out);
  EXPECT_EQ(x_in, x_out);
  EXPECT_EQ(Sx_in, Sx_out);
}

// Constructor test for Landmark1d
TEST_F(LandmarkTest, constructorTestLandmark1d){

  Landmark1d::Vec x_in, x_out;
  Landmark1d::Mat Sx_in, Sx_out;
  x_in << 0;
  Sx_in << 0;
  Landmark1d m;
  m.get(x_out);
  EXPECT_EQ(x_in, x_out);
  m.get(x_out, Sx_out);
  EXPECT_EQ(x_in, x_out);
  EXPECT_EQ(Sx_in, Sx_out);
}

// Information setting / getting test using Landmark2d
TEST_F(LandmarkTest, dataSetGetTest){

  Landmark2d::Vec x_in, x_out;
  Landmark2d::Mat Sx_in, Sx_out;

  x_in << 0, 0;
  Sx_in << 0, 0, 0, 0;
  Landmark2d m;
  m.get(x_out);
  EXPECT_EQ(x_in, x_out);
  m.get(x_out, Sx_out);
  EXPECT_EQ(x_in, x_out);
  EXPECT_EQ(Sx_in, Sx_out);

  x_in << 2.1, 1.2;
  Sx_in << 2.2, -1.1, -1.1, 2.2;
  m.set(x_in);
  m.get(x_out);
  EXPECT_EQ(x_in, x_out);

  m.set(x_in, Sx_in);
  m.get(x_out, Sx_out);
  EXPECT_EQ(x_in, x_out);
  EXPECT_EQ(Sx_in, Sx_out);
}

// Reference count test
TEST_F(LandmarkTest, referenceCountTest){

  Landmark2d m;
  unsigned int n = 0;
  
  for( int i = 0; i < 10; i++){
    EXPECT_EQ(i, m.getNRef() );
    EXPECT_EQ(i+1, m.incNRef() );
  }

  for( int i = 10; i > 0; i--){
    EXPECT_EQ(i, m.getNRef() );
    EXPECT_EQ(i-1, m.decNRef() );
  }
  EXPECT_EQ(0, m.getNRef() );

}

// Mahalanobis distance calculation test for Landmark2d
TEST_F(LandmarkTest, mahalanobisDistanceTest2d){

  Landmark2d::Vec x_in, x;
  Landmark2d::Mat Sx_in;
  
  Landmark2d m;

  x_in << 1, 2;
  Sx_in << 2, 0.0, 0.0, 2;
  
  x = x_in;
  m.set(x_in, Sx_in);
  EXPECT_EQ(0, RandomVecMathTools<Landmark2d>::mahalanobisDist(m, x) ); 

  x << 0, 2;
  EXPECT_EQ(sqrt(0.5), RandomVecMathTools<Landmark2d>::mahalanobisDist(m, x) ); 

  x << 0, 1;
  EXPECT_EQ(1, RandomVecMathTools<Landmark2d>::mahalanobisDist(m, x) ); 

  x << 4.4, 7.2;
  Sx_in << 1, -0.5, -0.5, 2;
  m.set(x_in, Sx_in);
  double d2 = (x - x_in).transpose() * Sx_in.inverse() * (x - x_in);
  EXPECT_EQ( sqrt(d2), RandomVecMathTools<Landmark2d>::mahalanobisDist(m, x) ); 
}

// Mahalanobis distance calculation test for Landmark1d
TEST_F(LandmarkTest, mahalanobisDistanceTest1d){
  
  Landmark1d::Vec x_in, x;
  Landmark1d::Mat Sx_in;
  
  Landmark1d m;
  EXPECT_EQ(1, m.getNDim());

  x_in << 1;
  Sx_in << 2;
  
  x << x_in;
  m.set(x_in, Sx_in);
  EXPECT_EQ(0, RandomVecMathTools<Landmark1d>::mahalanobisDist(m, x) ); 

  x << 0;
  EXPECT_EQ(sqrt(0.5), RandomVecMathTools<Landmark1d>::mahalanobisDist(m, x) ); 

  x << 3;
  EXPECT_EQ(sqrt(2), RandomVecMathTools<Landmark1d>::mahalanobisDist(m, x) ); 

  x << 4.4;
  Sx_in << 0.75;
  m.set(x_in, Sx_in);
  double d2 = (x - x_in).transpose() * Sx_in.inverse() * (x - x_in);
  EXPECT_EQ( sqrt(d2), RandomVecMathTools<Landmark1d>::mahalanobisDist(m, x) ); 
}

