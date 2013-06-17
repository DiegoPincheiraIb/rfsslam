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
  EXPECT_EQ(NULL, pf.getProcessModel() );
  EXPECT_EQ(NULL, pf.getMeasurementModel() );
}

// Test PF constructor
TEST_F(ParticleFilterTest, ParticleFilterConstructorTest){
  unsigned int n = 100;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::State initState;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::State::Vec x0;
  x0 << 1, 2, 3;
  initState.set(x0);
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel> 
    pf(n, initState, &motionModel, &measurementModel);
  EXPECT_EQ(n, pf.getParticleCount() );
  EXPECT_EQ(&motionModel, pf.getProcessModel() );
  EXPECT_EQ(&measurementModel, pf.getMeasurementModel() );
  EXPECT_EQ(n, pf.getParticleSet()->size() );
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::ParticleSet* pps = pf.getParticleSet();
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::State s;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::State::Vec x;
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
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::State initState;
  ParticleFilter<OdometryMotionModel2d, RangeBearingModel>::State::Vec x0;
  x0 << 1, 2, 3;
  initState.set(x0);
  typedef ParticleFilter<OdometryMotionModel2d, RangeBearingModel> PF_TYPE;
  PF_TYPE pf(n, initState, &motionModel, &measurementModel);
  
  OdometryMotionModel2d motionModel2;
  RangeBearingModel measurementModel2;
  EXPECT_EQ(&motionModel, pf.getProcessModel() );
  EXPECT_EQ(&measurementModel, pf.getMeasurementModel() );
  pf.setProcessModel( &motionModel2 );
  pf.setMeasurementModel( &measurementModel2 );
  EXPECT_EQ(&motionModel2, pf.getProcessModel() );
  EXPECT_EQ(&measurementModel2, pf.getMeasurementModel() );  
  
  pf.setEffectiveParticleCountThreshold(1.234);
  EXPECT_EQ(1.234, pf.getEffectiveParticleCountThreshold() );

  std::vector<PF_TYPE::Measure> z1, z2;
  PF_TYPE::Measure m1( (PF_TYPE::Measure::Vec() << 1, 0.1).finished() );
  PF_TYPE::Measure m2( (PF_TYPE::Measure::Vec() << 2, 0.2).finished() );
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

}

TEST_F(ParticleFilterTest, ParticleFilterSampleTest){

}
