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

#include "LinearAssignment.hpp"
#include <iostream>


int main(int argc, char* argv[]){
  
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
  }while(k < nbfa);

}
