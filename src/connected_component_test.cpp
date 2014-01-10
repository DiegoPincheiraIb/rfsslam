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

  int const C_size = 7;
  double** C = new double*[C_size];
  for(int i = 0; i < C_size; i++){
    C[i] = new double[C_size];
    for(int j = 0; j < C_size; j++){
      C[i][j] = 0;
    }
  }
  C[1][4] = gen();
  C[2][1] = gen();
  C[3][3] = gen();
  C[3][5] = gen();
  C[5][2] = gen();
  C[5][3] = gen();
  C[5][5] = gen();

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
    
    printf("Partition %d\n", n);

    size_t nCols, nRows;
    likelihoodMat.getPartitionSize(n, nRows, nCols);
    double** Cp = new double* [nCols + nRows];
    for(int i = 0; i < nCols + nRows; i++){
      Cp[i] = new double[nCols + nRows];
    }

    bool isZeroPartition;
    int* rowIdx = new int[nRows];
    int* colIdx = new int[nCols];
    likelihoodMat.getPartition(n, Cp, &isZeroPartition, rowIdx, colIdx);
 
    if(isZeroPartition)
      printf("Zero partition\n");
    printf("Original rows: ");
    for(int i = 0; i < nRows; i++){
      printf("%d  ", rowIdx[i]);
    }
    printf("\n");
    printf("Original cols: ");
    for(int j = 0; j < nCols; j++){
      printf("%d  ", colIdx[j]);
    }
    printf("\n");

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

    delete[] rowIdx;
    delete[] colIdx;

  }

  for(int i = 0; i < C_size; i++){
    delete[] C[i];
  }
  delete[] C;


  int* ordering = new int[5];
  for(int i = 0; i < 5; i++){
    ordering[i] = i;
  }
  do{
    printf("%d %d %d %d %d\n", ordering[0], ordering[1], ordering[2], ordering[3], ordering[4]);
  }while(std::next_permutation(ordering, ordering + 5));


  return 0;
}
