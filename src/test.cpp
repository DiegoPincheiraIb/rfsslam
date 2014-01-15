/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung, Felipe Inostroza
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Advanced Mining Technology Center (AMTC), the
 *       Universidad de Chile, nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AMTC, UNIVERSIDAD DE CHILE, OR THE COPYRIGHT 
 * HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
 * THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "BruteForceAssignment.hpp"
#include "HungarianMethod.hpp"
#include "MurtyAlgorithm.hpp"
#include <stdio.h>
#include <iostream>

int main(int argc, char* argv[]){
 
  int n = 5;
  double** C = new double* [n];
  for(int i = 0; i < n; i++ ){
    C[i] = new double[n];
  }
  C[0][0] =  1; C[0][1] = 1; C[0][2] =  9; C[0][3] =  14; C[0][4] =12;
  C[1][0] =  1; C[1][1] = 1; C[1][2] =  7; C[1][3] =  12; C[1][4] = 15;
  C[2][0] =  2; C[2][1] = 15; C[2][2] = 2; C[2][3] =  2; C[2][4] = 2;
  C[3][0] =  1; C[3][1] = 1; C[3][2] = 18; C[3][3] =  9; C[3][4] = 6;
  C[4][0] = 14; C[4][1] = 1; C[4][2] =  0; C[4][3] =  9; C[4][4] = 8;

  printf("\n");
  for(int i = 0; i < n; i++ ){
    for(int j = 0; j < n; j++ ){
      printf("%f   ", C[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  double** C_1;
  double** C_2;
  int* a1 = new int[n];
  double s1;
  int* i_remap;
  int* j_remap;
  CostMatrix CMat(C, n);
  CMat.reduce(10);
  int n1 = CMat.getCostMatrix(C_1);
  int n2 = CMat.getCostMatrixReduced(C_2, a1, &s1, i_remap, j_remap);

  printf("\n");
  for(int i = 0; i < n1; i++ ){
    for(int j = 0; j < n1; j++ ){
      printf("%f   ", C_1[i][j]);
    }
    printf("\n");
  }
  printf("\n");

  printf("\n");
  for(int i = 0; i < n2; i++ ){
    for(int j = 0; j < n2; j++ ){
      printf("%f   ", C_2[i][j]);
    }
    printf("\n");
  }
  printf("\n");
  
  printf("i_remap:\n");
  for(int i = 0; i < n2; i++ ){
    printf("%d   ", i_remap[i]);
  }
  printf("\n");
 
  printf("j_remap:\n");
  for(int i = 0; i < n2; i++ ){
    printf("%d   ", j_remap[i]);
  }
  printf("\n");

  printf("Fixed assignments:\n");
  for(int i = 0; i < n1; i++ ){
    printf("%d   ", a1[i]);
  }
  printf("\n");

  delete a1;




  HungarianMethod hm;
  int hm_soln[n];
  double hm_score = 0;
  hm.run(C, n, hm_soln, &hm_score);
  printf("Hungarian Method\n");
  printf("Best assignment:\n");
  for(int j = 0; j < n; j++){
    printf("x[%d] ----- y[%d]\n", j, hm_soln[j]);
  }
  printf("Score: %f\n\n", hm_score);


  hm.run(CMat, hm_soln, &hm_score);
  printf("Hungarian Method with CostMatrix\n");
  printf("Best assignment:\n");
  for(int j = 0; j < n; j++){
    printf("x[%d] ----- y[%d]\n", j, hm_soln[j]);
  }
  printf("Score: %f\n\n", hm_score);


  BruteForceLinearAssignment bf;
  unsigned int** bfa;
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
  }while(k < nbfa);
  

}
