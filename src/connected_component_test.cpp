#include <boost/config.hpp>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <utility>
#include "LinearAssignment.hpp"

#include <boost/timer/timer.hpp>
#include <boost/format.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

int main(int argc, char *argv[])
{

  boost::mt19937 rng;
  boost::uniform_01<> ud;
  boost::variate_generator<boost::mt19937, boost::uniform_01<> > gen(rng, ud);

  int const C_size = 5;
  double** C = new double*[C_size];
  for(int i = 0; i < C_size; i++){
    C[i] = new double[C_size];
    for(int j = 0; j < C_size; j++){
      C[i][j] = 0;
    }
  }
  C[0][3] = gen();
  C[1][0] = gen();
  C[2][2] = gen();
  C[2][4] = gen();
  C[4][1] = gen();
  C[4][2] = gen();
  C[4][4] = gen();

  for(int i = 0; i < C_size; i++){
    for(int j = 0; j < C_size; j++){
      printf("%f  ", C[i][j]);
    }
    printf("\n");
  }


  boost::timer::cpu_timer timer;

  CostMatrix likelihoodMat(C, C_size);
  int nP = likelihoodMat.partition();
  
  boost::timer::cpu_times t = timer.elapsed();
  std::cout << boost::format("Elapsed time: %1% [ns]\n") % t.wall;

  for(int n = 0; n < nP; n++){
    
    size_t nCols, nRows;
    likelihoodMat.getPartitionSize(n, nCols, nRows);
    double** Cp = new double* [nCols + nRows];
    for(int i = 0; i < nCols + nRows; i++){
      Cp[i] = new double[nCols + nRows];
    }
    likelihoodMat.getPartition(n, Cp);
    for(int i = 0; i < nRows; i++){
      for(int j = 0; j < nCols; j++){
	printf("%f  ", Cp[i][j]);
      }
      printf("\n");
    }
    
    for(int i = 0; i < nCols + nRows; i++){
      delete[] Cp[i];
    }
    delete[] Cp;

  }

  for(int i = 0; i < C_size; i++){
    delete[] C[i];
  }
  delete[] C;

  return 0;
}
