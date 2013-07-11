#include "gtest/gtest.h"
#include "KalmanFilter.hpp"
#include "State.hpp"
#include "Measurement.hpp"
#include <Eigen/Core>


/**
 * \class KalmanFilterTest
 * \brief Unit testing fixture for the Kalman Filter class
 * \author Felipe Inostroza
 */


class KalmanFilterTest : public ::testing::Test{

protected:

  static const double abs_error=1e-14;
  /** Constructor for setting up each test */
  KalmanFilterTest(){}
  
  /** Destructor */
  virtual ~KalmanFilterTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //

};



////////// Test Cases //////////

TEST_F(KalmanFilterTest, TestKalmanFilter2d_example){


  for(double robot_angle=-2*PI ; robot_angle<6*PI ; robot_angle+=PI/300){
  
    Pose2d::Vec robot_x = Pose2d::Vec::Zero();
    Pose2d::Mat robot_Sx = Pose2d::Mat::Zero();
    
    robot_x(2)=robot_angle;
    Landmark2d::Vec l2_x , expected_x;
    Landmark2d::Mat l2_Sx , expected_Sx;
  
    Pose2d robot_pose(robot_x , robot_Sx);
  
    Measurement2d::Vec z;
    Measurement2d::Mat Sz;
    Landmark2d::Vec l_x;
    Landmark2d::Mat l_Sx;
  
    z << 3.1 , PI/2 -robot_angle;
    Sz << 0.1 , 0 , 0 , 0.01;
    l_x << 0 , 3;
    l_Sx << 1 , 0 , 0 , 1;
  
  
    Landmark2d landmark(l_x , l_Sx);
    Landmark2d updatedLandmark;
  
    Measurement2d meas(z , Sz);
  
    StaticProcessModel<Landmark2d> processModel;

    RangeBearingModel measModel( 0.1 , 0.01 );
    measModel.config.probabilityOfDetection_ = 0.7;
    measModel.config.probabilityOfFalseAlarm_ = 0.3;
    measModel.config.rangeLim_=10;
    measModel.config.rangeLimBuffer_=1;
  
    RangeBearingKalmanFilter filter(&processModel, &measModel);
    
    filter.predict(landmark, updatedLandmark);
  
    filter.correct(robot_pose, meas, updatedLandmark, updatedLandmark);
  
    expected_x << 0.0 , 3.09090909090909;
    expected_Sx << 0.0825688073394495 ,	0 ,
                 0 ,	0.0909090909090909 ;

    updatedLandmark.get(l2_x , l2_Sx);
    
    EXPECT_DOUBLE_EQ(l2_x(0), expected_x(0));
    EXPECT_DOUBLE_EQ(l2_x(1), expected_x(1));
    EXPECT_DOUBLE_EQ(l2_Sx(0,0), expected_Sx(0,0));
    EXPECT_DOUBLE_EQ(l2_Sx(0,1), expected_Sx(0,1));
    EXPECT_DOUBLE_EQ(l2_Sx(1,0), expected_Sx(1,0));
    EXPECT_DOUBLE_EQ(l2_Sx(1,1), expected_Sx(1,1));
  }
  
  
}

TEST_F(KalmanFilterTest, TestKalmanFilterExample2){


  
  
  Pose2d::Vec robot_x = Pose2d::Vec::Zero();
  Pose2d::Mat robot_Sx = Pose2d::Mat::Zero();
    
  robot_x(2)=0;
  Eigen::Matrix<double, 4, 1> l2_x , expected_x;
  Eigen::Matrix<double, 4, 4> l2_Sx , expected_Sx;
  
  Pose2d robot_pose(robot_x , robot_Sx);
  
  Eigen::Matrix<double, 6, 1> z;
  Eigen::Matrix<double, 6, 6> Sz;
  Eigen::Matrix<double, 4, 1> l_x;
  Eigen::Matrix<double, 4, 4> l_Sx;
  
  z << 1 , 1 , 1 , 1 , 1 , 1 ;
  Sz << 
    1 , 0 , 0 , 0 , 0 , 0 ,
    0 , 1 , 0 , 0 , 0 , 0 ,
    0 , 0 , 1 , 0 , 0 , 0 ,
    0 , 0 , 0 , 1 , 0 , 0 ,
    0 , 0 , 0 , 0 , 1 , 0 ,
    0 , 0 , 0 , 0 , 0 , 1 ;

  l_x << 0 , 0 , 0 , 0 ;
  l_Sx << 
    1 , 0 , 0 , 0 ,
    0 , 1 , 0 , 0 ,
    0 , 0 , 1 , 0 ,
    0 , 0 , 0 , 1 ;

  
  
  Landmark< Eigen::Matrix<double, 4, 1> , Eigen::Matrix<double, 4, 4> >  
    landmark(l_x , l_Sx); 
  Landmark< Eigen::Matrix<double, 4, 1> , Eigen::Matrix<double, 4, 4> > 
    updatedLandmark;
  
  RandomVec< Eigen::Matrix<double, 6, 1> , Eigen::Matrix<double, 6, 6> > 
    meas(z , Sz);
  
  StaticProcessModel< RandomVec<Eigen::Matrix<double, 4, 1> , Eigen::Matrix<double, 4, 4> > > processModel;
  
  
  Eigen::Matrix<double, 6, 4> H;
  Eigen::Matrix<double, 6, 6> R;
  R=Sz;
  H <<
    1 , 0 , 0 , 0 ,
    0 , 1 , 0 , 0 ,
    0 , 0 , 1 , 0 ,
    0 , 0 , 0 , 1 ,
    1 , 0 , 0 , 0 ,
    0 , 1 , 0 , 0 ;


  LinearModel< Landmark< Eigen::Matrix<double, 4, 1> , Eigen::Matrix<double, 4, 4> > , RandomVec< Eigen::Matrix<double, 6, 1> , Eigen::Matrix<double, 6, 6> > >
    measModel(R , H);
  
  measModel.config.probabilityOfDetection_ = 0.7;
  measModel.config.probabilityOfFalseAlarm_ = 0.3;
  
  KalmanFilter<
    StaticProcessModel< RandomVec<Eigen::Matrix<double, 4, 1> , Eigen::Matrix<double, 4, 4> > > ,
    LinearModel< Landmark< Eigen::Matrix<double, 4, 1> , Eigen::Matrix<double, 4, 4> > , RandomVec< Eigen::Matrix<double, 6, 1> , Eigen::Matrix<double, 6, 6> > > >
    filter(&processModel, &measModel);
    
  filter.predict(landmark, updatedLandmark);
  
  filter.correct(robot_pose, meas, updatedLandmark, updatedLandmark);
  
  expected_x << 
    0.666666666666667,
    0.666666666666667,
    0.500000000000000,
    0.500000000000000;
  expected_Sx <<
    0.333333333333333, 0,                  0,                0,
    0,                 0.333333333333333,  0,                0,
    0,	               0,                  0.500000000000000,0,
    0,	               0,                  0,                0.500000000000000;


  
  updatedLandmark.get(l2_x , l2_Sx);
    
  for(int i=0;i<4;i++){
    EXPECT_NEAR(l2_x(i), expected_x(i), abs_error);
    for(int j=0;j<4;j++)
      EXPECT_NEAR(l2_Sx(i,j), expected_Sx(i,j), abs_error);

  }
  
  
  
}


// Constructor test for Landmark2d
