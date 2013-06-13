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

};

////////// Test Cases //////////
