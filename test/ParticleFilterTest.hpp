// Test for the Particle Filter Class
// Keith Leung 2013

#include "gtest/gtest.h"
#include "MeasurementModel.hpp"
#include "ProcessModel.hpp"

/**
 * \class ParticleFilterTest
 * \brief Unit testing fixture for the ParticleFilter class
 * \author Keith Leung 
 */

class ParticleFilterTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  ParticleFilterTest(){
    motionModel = OdometryMotionModel2d();
    measurementModel = RangeBearingModel();
  }
  
  /** Destructor */
  virtual ~ParticleFilterTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //
  OdometryMotionModel2d motionModel;
  RangeBearingModel measurementModel;
};

////////// Test Cases //////////

// TODO incomplete

