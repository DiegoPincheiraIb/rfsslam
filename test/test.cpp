// Unit testing program for the PHDFilter Library using Google Test
// Keith Leung 2013

#include <gtest/gtest.h>

#include "PoseTest.hpp"
#include "MeasurementTest.hpp"
#include "ProcessModelTest.hpp"
#include "MeasurementModelTest.hpp"
#include "GaussianMixtureTest.hpp"
#include "ParticleFilterTest.hpp"
#include "KalmanFilterTest.hpp"
#include "LinearAssignmentTest.hpp"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
