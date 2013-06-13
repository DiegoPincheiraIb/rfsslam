// Unit testing program for the PHDFilter Library using Google Test
// Keith Leung 2013

#include <gtest/gtest.h>

#include "LandmarkTest.hpp"
#include "MeasurementTest.hpp"
#include "MeasurementModelTest.hpp"
#include "PoseTest.hpp"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
