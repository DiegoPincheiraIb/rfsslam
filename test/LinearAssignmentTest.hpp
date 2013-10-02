// Test for Linear Assignment
// Keith Leung 2013

#include "LinearAssignment.hpp"

class LinearAssignmentTest : public ::testing::Test{

protected:

  /** Constructor for setting up each test */
  LinearAssignmentTest(){}
  
  /** Destructor */
  virtual ~LinearAssignmentTest(){}

  /** Setup -- called after constructor before each test */
  virtual void SetUp(){}

  /** Teardown -- called after each test but before destructor */
  virtual void TearDown(){}

};

TEST_F(LinearAssignmentTest, hungarianMethodTest){

  HungarianMethod hm;
  
  int n = 5;
  double** C = new double*[n];
  for(int i = 0; i < n; i++){
    C[i] = new double[n];
  }
  /*
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){  
      C[i][j] = 1;
    }
  }
  for(int i = 0; i < n; i++){
    C[i][i] = 5;
    }*/

  /*
  C[0][0] = 7;
  C[0][1] = 4;
  C[0][2] = 3;
  C[1][0] = 3;
  C[1][1] = 1;
  C[1][2] = 2;
  C[2][0] = 3;
  C[2][1] = 0;
  C[2][2] = 0;
  */

  C[0][0] = 10; C[0][1] = 0; C[0][2] =  8; C[0][3] = 15; C[0][4] = 0;
  C[1][0] = 10; C[1][1] = 18; C[1][2] =  7; C[1][3] = 17; C[1][4] = 24;
  C[2][0] = 13; C[2][1] = 16; C[2][2] =  9; C[2][3] = 14; C[2][4] = 15;
  C[3][0] = 12; C[3][1] = 12; C[3][2] =  8; C[3][3] = 18; C[3][4] = 6;
  C[4][0] = 14; C[4][1] = 17; C[4][2] = 10; C[4][3] = 19; C[4][4] = 8;

  int* s = new int[n];
  double c;
  hm.run(C, n, s, &c, true); 
  
  printf("Best Solution:\n");
  for(int i = 0; i < n; i++){
    printf("x[%d] ----- y[%d]\n", i, s[i]);
  }
  printf("Score: %f\n", c);
  

  for(int i = 0; i < n; i++){
    delete[] C[i];
  }
  delete[] C;
  delete[] s;

}


TEST_F(LinearAssignmentTest, MurtyAlgorithmTest){

  int n = 5;
  double** C = new double*[n];
  for(int i = 0; i < n; i++){
    C[i] = new double[n];
  }

  C[0][0] = 10; C[0][1] = 21; C[0][2] = 11; C[0][3] = 15; C[0][4] = 21;
  C[1][0] = 10; C[1][1] = 18; C[1][2] =  7; C[1][3] = 16; C[1][4] = 5;
  C[2][0] = 17; C[2][1] = 15; C[2][2] = 20; C[2][3] = 14; C[2][4] = 15;
  C[3][0] = 12; C[3][1] = 12; C[3][2] =  8; C[3][3] =  9; C[3][4] = 6;
  C[4][0] = 14; C[4][1] = 17; C[4][2] = 10; C[4][3] = 19; C[4][4] = 8;

  /*  
  C[0][0] = 10; C[0][1] = 19; C[0][2] =  8; C[0][3] = 15; C[0][4] = 21;
  C[1][0] = 10; C[1][1] = 18; C[1][2] =  7; C[1][3] = 17; C[1][4] = 24;
  C[2][0] = 13; C[2][1] = 16; C[2][2] =  9; C[2][3] = 14; C[2][4] = 15;
  C[3][0] = 12; C[3][1] = 12; C[3][2] =  8; C[3][3] = 18; C[3][4] = 6;
  C[4][0] = 14; C[4][1] = 17; C[4][2] = 10; C[4][3] = 19; C[4][4] = 8;
  */

  BruteForceLinearAssignment bf;
  int** bfa;
  double* bfs;
  int nbfa = bf.run(C, n, bfa, bfs);
  printf("Brute force approach looked through %d assignments\n", nbfa);
  printf("Best assignment:\n");
  for(int j = 0; j < n; j++){
    printf("x[%d] ----- y[%d]\n", j, bfa[0][j]);
  }
  printf("Score: %f\n", bfs[0]);
  printf("Worst assignment:\n");
  for(int j = 0; j < n; j++){
    printf("x[%d] ----- y[%d]\n", j, bfa[nbfa-1][j]);
  }
  printf("Score: %f\n\n\n", bfs[nbfa-1]);

  int* a;
  double score;
  int k;
  Murty murty(C, n);

  do
  {
    k = murty.findNextBest(a, &score); 
    printf("\nThe %d-best solution:\n", k);
    for(int i = 0; i < n; i++){
      printf("x[%d] ----- y[%d]\n", i, a[i]);
    }
    printf("Score: %f\n", score);
    ASSERT_EQ(score, bfs[k-1]);
  }while(k < nbfa);

  for(int i = 0; i < 2; i++){
    k = murty.findNextBest(a, &score); 
    EXPECT_EQ(score, 0);
    EXPECT_TRUE(a == NULL);
    EXPECT_EQ(k, -1);
  }
}
