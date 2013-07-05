// Test for classes derived from Measurement
// Keith Leung 2013

#include "gtest/gtest.h"
#include "Measurement.hpp"
#include "RandomVecMathTools.hpp"

/**
 * \class MeasurementTest
 * \brief Unit testing fixture for Measurement derived classes
 * \author Keith Leung 
 */
class MeasurementTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  MeasurementTest(){}
  
  /** Destructor */
  virtual ~MeasurementTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //

};

////////// Test Cases //////////

// Test constructor for Odometry2d
TEST_F(MeasurementTest, testConstructorOdometry2d){

  Eigen::Vector3d u_in, u_out;
  Eigen::Matrix3d Su_in, Su_out;
  double t_in, t_out;

  double dx = 1.1;
  double dy = 2.2;
  double dz = 3.3;
  double vdx = 1.0;
  double vdy = 0.5;
  double vdz = 0.1;
  double t = 0.345;
  
  u_in << dx, dy, dz;
  Su_in << vdx, 0, 0, 0, vdy, 0, 0, 0, vdz;
  t_in = t;

  Odometry2d u1(dx, dy, dz, vdx, vdy, vdz, t);
  u1.get(u_out, Su_out, t_out);
  EXPECT_EQ(u_in, u_out);
  EXPECT_EQ(Su_in, Su_out);
  EXPECT_EQ(t_in, t_out);
  u1.get(u_out, t_out);
  EXPECT_EQ(u_in, u_out);
  EXPECT_EQ(t_in, t_out);

  Odometry2d u2;
  u2.set(u_in, Su_in, t_in);
  u2.get(u_out, Su_out, t_out);
  EXPECT_EQ(u_in, u_out);
  EXPECT_EQ(Su_in, Su_out);
  EXPECT_EQ(t_in, t_out);
  u2.get(u_out, t_out);
  EXPECT_EQ(u_in, u_out);
  EXPECT_EQ(t_in, t_out);

  EXPECT_EQ( 0, RandomVecMathTools<Odometry2d>::mahalanobisDist(u2, u1));
}

// Test Measurement2d Class
TEST_F(MeasurementTest, testMeasurement2d){
  
  Eigen::Vector2d z, z_out, z_out2;
  Eigen::Matrix2d Sz, Sz_out;
  double t = 0.345;
  double t_out, t_out2;
  z << 0.6, 9.6;
  Sz << 0.9, 0.5, 0.5, 0.4;
  
  Measurement2d y;
  y.set(z, Sz, t);
  y.get(z_out, Sz_out, t_out);
  EXPECT_EQ(z, z_out);
  EXPECT_EQ(Sz, Sz_out);
  EXPECT_EQ(t, t_out); 
  y.get(z_out2, t_out2);
  EXPECT_EQ(z, z_out2);
  EXPECT_EQ(t, t_out2); 

  EXPECT_EQ(0, RandomVecMathTools<Measurement2d>::mahalanobisDist(y, z) );
  EXPECT_EQ(0, RandomVecMathTools<Measurement2d>::mahalanobisDist(y, y) );
  z << 0.1, 9.5;
  EXPECT_LT(0, RandomVecMathTools<Measurement2d>::mahalanobisDist(y, z) );
  EXPECT_EQ(0, RandomVecMathTools<Measurement2d>::mahalanobisDist(y, y) );
  EXPECT_LT(0, RandomVecMathTools<Measurement2d>::evalGaussianLikelihood(y, z) );
  EXPECT_LT(0, RandomVecMathTools<Measurement2d>::evalGaussianLikelihood(y, y) );
}

// Test Measurement1d Class
TEST_F(MeasurementTest, testMeasurement1d){

  typedef Eigen::Matrix<double, 1, 1> v1d;

  v1d z, z_out, z_out2;
  v1d Sz, Sz_out;
  double t, t_out, t_out2;
  z << 0.6;
  Sz << 0.9;
  t = 0.345;
  
  Measurement1d y;
  EXPECT_EQ(1, y.getNDim());
  y.set(z, Sz, t);
  y.get(z_out, Sz_out, t_out);
  EXPECT_EQ(z, z_out);
  EXPECT_EQ(Sz, Sz_out);
  EXPECT_EQ(t, t_out); 
  y.get(z_out2, t_out2);
  EXPECT_EQ(z, z_out2);
  EXPECT_EQ(t, t_out2); 

  EXPECT_EQ(0, RandomVecMathTools<Measurement1d>::mahalanobisDist(y, z) );
  EXPECT_EQ(0, RandomVecMathTools<Measurement1d>::mahalanobisDist(y, y) );
  z << 9.5;
  EXPECT_LT(0, RandomVecMathTools<Measurement1d>::mahalanobisDist(y, z) );
  EXPECT_EQ(0, RandomVecMathTools<Measurement1d>::mahalanobisDist(y, y) );
  EXPECT_LT(0, RandomVecMathTools<Measurement1d>::evalGaussianLikelihood(y, z) );
  EXPECT_LT(0, RandomVecMathTools<Measurement1d>::evalGaussianLikelihood(y, y) );
}
