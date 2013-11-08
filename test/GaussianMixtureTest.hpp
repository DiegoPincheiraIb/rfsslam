#include "GaussianMixture.hpp" 

class GaussianMixtureTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  GaussianMixtureTest(){}
  
  /** Destructor */
  virtual ~GaussianMixtureTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

  // Additional objects to declare //

};

TEST_F(GaussianMixtureTest, GaussianMixtureConstructorTest){

  // Using Landmark2d for template 
  GaussianMixture<Landmark2d> gm;
  EXPECT_EQ(0, gm.getGaussianCount() ); 
}


TEST_F(GaussianMixtureTest, GaussianMixtureAddRemoveGetTest){

  // Using Landmark2d for template 
  GaussianMixture<Landmark2d> gm;
  
  int nLandmarks = 10;
  std::vector<Landmark2d*> landmarks(nLandmarks);
  for(int i = 0; i < nLandmarks; i++){
    Landmark2d::Vec x;
    x << i, i;
    Landmark2d::Mat Sx = Landmark2d::Mat::Identity() * 2;
    landmarks[i] = new Landmark2d;
    landmarks[i]->set( x, Sx );
    gm.addGaussian( landmarks[i], i );
 
  }
  
  // get Gaussian count
  EXPECT_EQ( nLandmarks, gm.getGaussianCount() );

  // get test for state and covariance
  for(int i = 0; i < nLandmarks; i++){
     
    Landmark2d* lm;
    Landmark2d::Vec x, x_expect;
    Landmark2d::Mat Sx;
    double w;
    gm.getGaussian( i, lm, w );
    EXPECT_EQ(i, w);
    lm->get(x, Sx);
    x_expect << i, i;
    EXPECT_EQ(x_expect, x);
    EXPECT_EQ(Landmark2d::Mat::Identity() * 2, Sx);
 
  }

  // add Gaussian with memory allocation
  for( int i = 0; i < 5; i++ ){
    Landmark2d lm;
    Landmark2d::Vec x;
    x << i, i;
    Landmark2d::Mat Sx = Landmark2d::Mat::Identity();
    lm.set(x, Sx);
    gm.addGaussian( &lm, 1.23, true ); // memory should be allocated
  }

  EXPECT_EQ( nLandmarks + 5, gm.getGaussianCount() );

  for( int i = 0; i < 5; i++ ){
    Landmark2d* lm;
    Landmark2d::Vec x, x_expect;
    Landmark2d::Mat Sx;
    double w;
    gm.getGaussian( 10 + i, lm, w );
    EXPECT_EQ(1.23, w);
    lm->get(x, Sx);
    x_expect << i, i;
    EXPECT_EQ(x_expect, x);
    EXPECT_EQ(Landmark2d::Mat::Identity(), Sx) << "\n\n" << Sx << "\n\n";
  }

  // set weight
  for(int i = 0; i < nLandmarks; i++){
    gm.setWeight(i, 15 - i*0.25);
  }

  // get weight
  for(int i = 0; i < nLandmarks; i++){
    Landmark2d* lm;
    double w;
    gm.getGaussian(i, lm, w);
    EXPECT_EQ(15 - i*0.25, w);
  }

  // remove some Gaussians
  EXPECT_EQ(14, gm.removeGaussian(2) );
  EXPECT_EQ(13, gm.removeGaussian(6) );
  EXPECT_EQ(12, gm.removeGaussian(7) );
  EXPECT_EQ(12, gm.removeGaussian(7) );
  EXPECT_EQ(12, gm.removeGaussian(7) );
  EXPECT_EQ(11, gm.removeGaussian(10) );
  EXPECT_EQ(11, gm.removeGaussian(10) );
  EXPECT_EQ(10, gm.removeGaussian(13) );

  EXPECT_EQ(nLandmarks, gm.getGaussianCount());
  
  Landmark2d* l_temp;
  double w;
  gm.getGaussian(2, l_temp, w); 
  EXPECT_EQ(NULL, l_temp);
  EXPECT_EQ(0, w);
  gm.getGaussian(6, l_temp, w); 
  EXPECT_EQ(NULL, l_temp);
  EXPECT_EQ(0, w);
  gm.getGaussian(7, l_temp, w); 
  EXPECT_EQ(NULL, l_temp);
  EXPECT_EQ(0, w);
  gm.getGaussian(10, l_temp, w); 
  EXPECT_EQ(NULL, l_temp);
  EXPECT_EQ(0, w);
  gm.getGaussian(13, l_temp, w); 
  EXPECT_EQ(NULL, l_temp);
  EXPECT_EQ(0, w);

}

// test Gaussian update

TEST_F(GaussianMixtureTest, GaussianMixtureUpdateTest){

  // Using Landmark2d for template 
  GaussianMixture<Landmark2d> gm;

  Landmark2d* lm1 = new Landmark2d; // gm will take care of deallocating mem
  Landmark2d* lm2 = new Landmark2d; // gm will take care of deallocating mem
  Landmark2d* lm3 = new Landmark2d;
  Landmark2d::Vec x;
  Landmark2d::Mat Sx = Landmark2d::Mat::Identity();
  x << 0, 0;
  lm1->set(x, Sx);
  lm2->set(x, Sx);  
  x << 1, 1;
  Sx = Landmark2d::Mat::Identity() * 2;
  lm3->set(x, Sx);  

  gm.addGaussian( lm1, 1 );
  gm.addGaussian( lm2, 2 );

  gm.updateGaussian(0, *lm3);

  Landmark2d::Vec x1, x2;
  Landmark2d::Mat Sx1, Sx2;
  double w1, w2;

  gm.getGaussian(0, lm1, w1);  
  
  lm1->get(x1, Sx1);
  EXPECT_EQ(1, w1);
  EXPECT_EQ(x, x1);
  EXPECT_EQ(Sx, Sx1);
  EXPECT_NE(x, x2);
  EXPECT_NE(Sx, Sx2);

  gm.updateGaussian(1, *lm3, 10);
  gm.getGaussian(1, lm2, w2);
  lm2->get(x2, Sx2);
  EXPECT_EQ(10, w2);
  EXPECT_EQ(x, x2);
  EXPECT_EQ(Sx, Sx2);  

  delete lm3;

}

// test merging

TEST_F(GaussianMixtureTest, GaussianMixtureMergeTest){

  // Using Landmark2d for template 
  GaussianMixture<Landmark2d> gm;
  
  int nLandmarks = 10;
  std::vector<Landmark2d*> landmarks(nLandmarks);
  for(int i = 0; i < nLandmarks; i++){
    Landmark2d::Vec x;
    x << i, i;
    Landmark2d::Mat Sx = Landmark2d::Mat::Identity() * 2;
    landmarks[i] = new Landmark2d;
    landmarks[i]->set( x, Sx );
    double w = i;
    gm.addGaussian( landmarks[i], w );
 
  }
  
  Landmark2d::Vec x;
  x << 0.1, 0.1;
  Landmark2d::Mat Sx = Landmark2d::Mat::Identity() * 2;
  Landmark2d lm(x, Sx);
  gm.updateGaussian(5, lm);

  // Expect 0 and 5 to merge
  EXPECT_EQ(1, gm.merge());
  EXPECT_EQ(nLandmarks - 1, gm.getGaussianCount() );
  Landmark2d* plm;
  double w;
  gm.getGaussian(5, plm, w);
  EXPECT_EQ(NULL, plm);
  EXPECT_EQ(0, w);

  // Expect Gaussian 0 to become Gaussian 5 since
  // weight of Gaussian 0 is 0
  gm.getGaussian(0, plm, w);
  EXPECT_EQ(5, w);
  Landmark2d::Vec x1;
  Landmark2d::Mat Sx1;
  plm->get(x1, Sx1);  
  EXPECT_EQ(x, x1);
  EXPECT_EQ(Sx, Sx1);

  x << 1.2, 1.2;
  Sx = Landmark2d::Mat::Identity() * 2;
  lm.set(x, Sx);
  gm.updateGaussian(6, lm, 1);

  x << 3.2, 3.2;
  Sx = Landmark2d::Mat::Identity() * 2;
  lm.set(x, Sx);
  gm.updateGaussian(7, lm, 3);
  
  // Expect 1 and 6 to merge
  // Expect 3 and 7 to merge
  EXPECT_EQ(2, gm.merge());
  EXPECT_EQ(nLandmarks - 3, gm.getGaussianCount() );

  gm.getGaussian(1, plm, w);
  EXPECT_EQ(2, w);
  plm->get(x1, Sx1);  
  x << 1.1 , 1.1;
  Sx << 2.01, 0.01, 0.01, 2.01;
  EXPECT_EQ(x, x1);
  for(int i = 0; i < Sx.rows(); i++){
    for(int j = 0; i < Sx.rows(); i++){
      EXPECT_NEAR(0, Sx(i,j) - Sx1(i,j), 1e-15);
    }
  }

  gm.getGaussian(3, plm, w);
  EXPECT_EQ(6, w);
  plm->get(x1, Sx1);  
  x << 3.1 , 3.1;
  Sx << 2.01, 0.01, 0.01, 2.01;
  EXPECT_EQ(x, x1);
  for(int i = 0; i < Sx.rows(); i++){
    for(int j = 0; i < Sx.rows(); i++){
      EXPECT_NEAR(0, Sx(i,j) - Sx1(i,j), 1e-15);
    }
  }
}

// Test pruning 
TEST_F(GaussianMixtureTest, GaussianMixturePruningTest){
  
  // Using Landmark2d for template 
  GaussianMixture<Landmark2d> gm;
  
  int nLandmarks = 10;
  std::vector<Landmark2d*> landmarks(nLandmarks);
  for(int i = 0; i < nLandmarks; i++){
    Landmark2d::Vec x;
    x << i, i;
    Landmark2d::Mat Sx = Landmark2d::Mat::Identity() * 2;
    landmarks[i] = new Landmark2d;
    landmarks[i]->set( x, Sx );
    double w = i + 1;
    gm.addGaussian( landmarks[i], w );
  }
  EXPECT_EQ(0, gm.prune(1));
  EXPECT_EQ(10, gm.getGaussianCount());
  
  Landmark2d::Vec x;
  x << 0.1, 0.1;
  Landmark2d::Mat Sx = Landmark2d::Mat::Identity() * 2;
  Landmark2d lm(x, Sx);
  gm.updateGaussian(5, lm, 0.1); // expect to get pruned
  gm.updateGaussian(6, lm, 0.3); // expect to get pruned
  gm.updateGaussian(8, lm, 0.4); // expect to get pruned
  gm.updateGaussian(2, lm, 0.7);
  EXPECT_EQ(3, gm.prune(0.5));
  EXPECT_EQ(7, gm.getGaussianCount());

  for(int i = 0; i < nLandmarks; i++){ // ask to update 10
    bool result = gm.updateGaussian(i, lm, 0.1);
    if( i >= 7 )
      EXPECT_FALSE(result); // gm should only contain 7 landmark
    else
      EXPECT_TRUE(result); 
  }
  // add 3 Gaussians back
  gm.addGaussian( &lm, 0.1, true);
  gm.addGaussian( &lm, 0.1, true);
  gm.addGaussian( &lm, 0.1, true);
  EXPECT_EQ(10, gm.prune(0.2)); // everything should get pruned
  EXPECT_EQ(0, gm.getGaussianCount());

}


// Test with shared landmarks from 2 mixtures
TEST_F(GaussianMixtureTest, GaussianMixtureCopyTest){

  // Using Landmark2d for template 
  GaussianMixture<Landmark2d> gm1;

  int nLandmarks = 10;
  std::vector<Landmark2d*> landmarks(nLandmarks);
  for(int i = 0; i < nLandmarks; i++){
    Landmark2d::Vec x;
    x << i, i;
    Landmark2d::Mat Sx = Landmark2d::Mat::Identity() * 2;
    landmarks[i] = new Landmark2d;
    landmarks[i]->set( x, Sx );
    double w = i + 1;
    gm1.addGaussian( landmarks[i], w );
  }

  GaussianMixture<Landmark2d> gm2( gm1 );

  EXPECT_EQ(gm1.getGaussianCount(), gm2.getGaussianCount());

  Landmark2d* lm1;
  Landmark2d* lm2;
  double w1, w2;
  gm1.getGaussian(3, lm1, w1);
  gm2.getGaussian(3, lm2, w2);
  EXPECT_NE(lm1, lm2);

  gm1.updateGaussian(1, *lm1);
  gm2.getGaussian(1, lm2, w2);
  EXPECT_NE(lm1, lm2);

  for(int i = 0; i < gm1.getGaussianCount(); i++){
    gm1.removeGaussian(i);
  }

  for(int i = 0; i < gm1.getGaussianCount(); i++){
    gm2.getGaussian(i, lm2); // should not crash
  }

}
