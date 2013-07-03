// Test for the Particle Filter Class
// Keith Leung 2013

#include "gtest/gtest.h"
#include "MeasurementModel.hpp"
#include "ProcessModel.hpp"
#include "ParticleFilter.hpp"

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

// Test Particle constructor
TEST_F(ParticleFilterTest, ParticleDefaultConstructorTest){
  Particle<Pose2d> p;
  Particle<Pose2d>::tPose x_in, x_out;
  Particle<Pose2d>::tPose::Vec v_in, v_out;
  EXPECT_EQ(0, p.getId() );
  EXPECT_EQ(0, p.getParentId() );
  EXPECT_EQ(0, p.getWeight() );
  x_in.get(v_in); 
  p.getPose(x_out);
  x_out.get(v_out);
  EXPECT_EQ(v_in, v_out );
}

TEST_F(ParticleFilterTest, ParticleConstructorTest){
  Particle<Pose2d>::tPose x_in( 1, 2, 1);
  Particle<Pose2d>::tPose x_out;
  Particle<Pose2d> p( 123, x_in, 3);
  EXPECT_EQ(123, p.getId() );
  EXPECT_EQ(123, p.getParentId() );
  EXPECT_EQ(3, p.getWeight() );
  Particle<Pose2d>::tPose::Vec v_in, v_out;
  x_in.get(v_in); 
  p.getPose(x_out);
  x_out.get(v_out);
  EXPECT_EQ(v_in, v_out );
}

TEST_F(ParticleFilterTest, ParticleSetGetTest){
  Particle<Pose2d> p;
  Particle<Pose2d>::tPose x_in( 1, 2, 1);
  Particle<Pose2d>::tPose x_out;
  Particle<Pose2d>::tPose::Vec v_in, v_out;
  p.setPose(x_in);
  p.getPose(x_out);
  x_in.get(v_in); 
  x_out.get(v_out);
  EXPECT_EQ(v_in, v_out );

  p.setWeight(2.2);
  EXPECT_EQ(2.2, p.getWeight() );
  
  p.setId( 5 );
  EXPECT_EQ(5, p.getId() );

  p.setParentId( 6);
  EXPECT_EQ(6, p.getParentId() );
  
}

TEST_F(ParticleFilterTest, ParticleCopyTest){
  Particle<Pose2d>::tPose x_in( 1, 2, 1);
  Particle<Pose2d>::tPose x_out;
  Particle<Pose2d>::tPose::Vec v_in, v_out;
  Particle<Pose2d> p1( 3, x_in, 1.1);
  Particle<Pose2d> p2;
  p1.copyStateTo( &p2 );
  p2.getPose(x_out);
  x_in.get(v_in); 
  x_out.get(v_out);
  EXPECT_EQ(v_in, v_out );
  EXPECT_EQ(0, p2.getWeight() );
  EXPECT_EQ(0, p2.getId() );
  EXPECT_EQ(0, p2.getParentId() );
}

// Test PF constructor
TEST_F(ParticleFilterTest, ParticleFilterDefaultConstructorTest){
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel> pf;
  EXPECT_EQ(0, pf.getParticleCount() );
}

// Test PF constructor
TEST_F(ParticleFilterTest, ParticleFilterConstructorTest){
  unsigned int n = 100;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose initState;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose::Vec x0;
  x0 << 1, 2, 3;
  initState.set(x0);
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel> pf(n, &initState);
  EXPECT_EQ(n, pf.getParticleCount() );
  EXPECT_EQ(n, pf.getParticleSet()->size() );
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TParticleSet* pps = pf.getParticleSet();
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose s;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose::Vec x;
  for(int i = 0; i < pps->size(); i++){
    (*pps)[i]->getPose(s);
    s.get(x);
    EXPECT_EQ(x0, x);
  }
}


// Test PF set / get functions
// already tested with constructor test:
// getParticleCount() 
// getParticleSet()
TEST_F(ParticleFilterTest, ParticleFilterSetGetTest){
  unsigned int n = 100;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose initState;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose::Vec x0;
  x0 << 1, 2, 3;
  initState.set(x0);
  typedef ParticleFilter<OdometryMotionModel2d, RangeBearingModel> PF_TYPE;
  PF_TYPE pf(n);
  
  pf.setEffectiveParticleCountThreshold(1.234);
  EXPECT_EQ(1.234, pf.getEffectiveParticleCountThreshold() );

  std::vector<PF_TYPE::TMeasure> z1, z2;
  PF_TYPE::TMeasure m1( (PF_TYPE::TMeasure::Vec() << 1, 0.1).finished() );
  PF_TYPE::TMeasure m2( (PF_TYPE::TMeasure::Vec() << 2, 0.2).finished() );
  z1.push_back(m1);
  z1.push_back(m1);
  z1.push_back(m1);
  z2.push_back(m2);
  z2.push_back(m2);
  z2.push_back(m2);
  z2.push_back(m2);
  z2.push_back(m2);
  z2.push_back(m2);
  pf.setMeasurements( z1 );
  EXPECT_EQ(0, z1.size());
  pf.setMeasurements( z2 );
  EXPECT_EQ(0, z2.size());
}

TEST_F(ParticleFilterTest, ParticleFilterPropagateTest){

  typedef ParticleFilter<OdometryMotionModel2d, RangeBearingModel> PF_TYPE;

  PF_TYPE::TPose::Mat ProcessNoise;
  ProcessNoise << 3, 2, 1, 2, 4, -1, 1, -1, 5;

  unsigned int n = 200000;
  PF_TYPE::TPose initState;
  PF_TYPE::TPose::Vec x0;
  x0 << 0, 0, 0; // no rotation for this test
  initState.set(x0);
  PF_TYPE pf(n, &initState);
  pf.getProcessModel()->setNoise( ProcessNoise );

  PF_TYPE::TParticleSet* pps = pf.getParticleSet();
  PF_TYPE::TPose s;
  PF_TYPE::TPose::Vec x;
  for(int i = 0; i < pps->size(); i++){
    (*pps)[i]->getPose(s);
    s.get(x);
    EXPECT_EQ(x0, x);
  }
  
  Odometry2d::Vec u;
  Odometry2d::Mat Su;
  Odometry2d odo;
  u << 1, 2, 0; // no rotation for this test
  Su = Odometry2d::Mat::Identity();
  odo.set(u, Su, 0.1);

  pf.propagate( odo );
  
  std::vector < PF_TYPE::TPose::Vec > v;
  for(int i = 0; i < pps->size(); i++){
    (*pps)[i]->getPose(s);
    s.get(x);
    v.push_back(x);    
  }

  PF_TYPE::TPose::Vec sum;
  sum = PF_TYPE::TPose::Vec::Zero();
  for(int i = 0; i < n; i++){
    sum = sum + v[i];
  }  
  PF_TYPE::TPose::Vec mean = sum / n;

  PF_TYPE::TPose::Vec e;
  e = mean - (x0 + u);

  for( int i = 0; i < mean.size(); i++){
    EXPECT_NEAR(0.0, e(i), 1e-2);
  }
  
  PF_TYPE::TPose::Mat sum2;
  sum2 = Pose2d::Mat::Zero();
  for(int i = 0; i < n; i++){
    e = v[i] - mean;
    sum2 = sum2 + e * e.transpose();
  }  
  PF_TYPE::TPose::Mat cov = sum2 / n;  

  for( int i = 0; i < mean.size(); i++){
    for( int j = 0; j < mean.size(); j++){
      double e = cov(i,j) - ProcessNoise(i,j);
      EXPECT_NEAR(0.0, e, 1e-2) << "Failure may be due to randomness\nExpected:" 
				<< ProcessNoise(i,j)
				<< "\nGot:" << cov(i,j) << "\n\n";
    }
  }

}

TEST_F(ParticleFilterTest, ParticleFilterResampleTest){

  typedef ParticleFilter<OdometryMotionModel2d, RangeBearingModel> PF_TYPE;

  unsigned int n = 10;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose initState;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::TPose::Vec x0;
  x0 << 1, 2, 3;
  initState.set(x0);
  PF_TYPE pf(n, &initState);
  
  PF_TYPE::TParticleSet* pps = pf.getParticleSet();
  PF_TYPE::TPose s;
  PF_TYPE::TPose::Vec x;
  for(int i = 0; i < pps->size(); i++){
    (*pps)[i]->setWeight(1);
  }

  // All particles have equal weight, Effective_n_of_particles = n / 4 by default
  // We should not resample
  EXPECT_FALSE( pf.resample() );

  // Expect to get n / 2 = 5 samples of particle 5
  pf.setEffectiveParticleCountThreshold( n ); // force resample
  for(int i = 0; i < pps->size(); i++){
    (*pps)[i]->setWeight(1);
  }
  (*pps)[1]->setWeight(2);
  (*pps)[5]->setWeight(10);
  pf.resample();
  int count = 0;
  for(int i = 0; i < pps->size(); i++){
    //std::cout << (*pps)[i]->getParentId() << std::endl;
    if ((*pps)[i]->getParentId() == 5){
      count++;
    }
  }
  EXPECT_EQ(5, count);

  // Expect to get 4 samples of particle 2 and 4 of particle 6
  pf.setEffectiveParticleCountThreshold( n ); // force resample
  for(int i = 0; i < pps->size(); i++){
    (*pps)[i]->setWeight(1);
  }
  (*pps)[0]->setWeight(20);
  (*pps)[6]->setWeight(20);
  (*pps)[7]->setWeight(3);
  pf.resample();
  count = 0;
  for(int i = 0; i < pps->size(); i++){
    if ((*pps)[i]->getParentId() == 0){
      count++;
    }
  }
  EXPECT_EQ(4, count);

  count = 0;
  for(int i = 0; i < pps->size(); i++){
    if ((*pps)[i]->getParentId() == 6){
      count++;
    }
  }
  EXPECT_EQ(4, count);

}
