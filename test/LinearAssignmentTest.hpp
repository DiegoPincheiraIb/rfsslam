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

  HungarianMethod hm();
  
  n = 3;
  double** C = new double*[n];
  for(int i = 0; i < n; i++){
    C[n] = new double[n];
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){  
      C[i][j] = 1;
    }
  }
  for(int i = 0; i < n; i++){
    C[i][i] = 5;
  }

  int* s = new int[n];
  double c;
  hm.run(C, n, s, &c); 
  
  printf("Solution:\n");
  for(int i = 0; i < n; i++){
    printf("x[%d] ----- y[%d]\n", i, s[i]);
  }
  printf("Score: %f\n", c);

  for(int i = 0; i < n; i++){
    delete[] C[n];
  }
  delete[] C;
  delete[] s;

}
