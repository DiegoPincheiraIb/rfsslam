// Test for the RB-PHD Filter Class
// Keith Leung 2013

#include "gtest/gtest.h"
#include "RBPHDFilter.hpp"

/**
 * \class RBPHDFilterTest
 * \brief Unit testing fixture for the ParticleFilter class
 * \author Keith Leung 
 */

class RBPHDFilterTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  RBPHDFilterTest(){
    motionModel = OdometryMotionModel2d();
    measurementModel = RangeBearingModel();
  }
  
  /** Destructor */
  virtual ~RBPHDFilterTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //
  OdometryMotionModel2d motionModel;
  RangeBearingModel measurementModel;


};

////////// Test Cases //////////

// Test Particle constructor
TEST_F(RBPHDFilterTest, RBPHDFilterDefaultConstructorTest){

  const unsigned int n = 100;
  RBPHDFilter::TPose::Vec x0 = RBPHDFilter::TPose::Vec::Zeros();
  RBPHDFilter::TPose::Vec Sx0 = RBPHDFilter::TPose::Vec::Identity() * 0.1;
  RBPHDFilter::TPose p0( x0, Sx0, 0.1);
  RBPHDFilter filter(n, p0, motionModel, measurementModel );

}
