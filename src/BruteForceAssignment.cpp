/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2013, Keith Leung
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

#include <cmath>
#include <stdio.h>

BruteForceLinearAssignment::BruteForceLinearAssignment() : a_(NULL), s_(NULL), nAssignments(0){}

BruteForceLinearAssignment::~BruteForceLinearAssignment(){
  if(a_ != NULL ){
    for( unsigned int m = 0; m < nAssignments; m++ ){
      delete[] a_[m];
    }
    delete[] a_;
  }
  if(s_ != NULL ){
    delete[] s_;
  }
  s_ = new double [nAssignments];
}

unsigned int BruteForceLinearAssignment::run(double** C, int n, int** &a, double* &s, bool maxToMin){

  double offset = 0;
  if(!maxToMin){
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	C[x][y] *= -1;
	if(C[x][y] < offset)
	  offset = C[x][y];
      }
    }
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	C[x][y] -= offset;
      }
    }
  }


  unsigned int nAssignments = 0;
  bool lastPermutationSequence = false;
  std::vector<int> currentPermutation(n);
  for(int m = 0; m < n; m++){
    currentPermutation[m] = m;
  }

  while( !lastPermutationSequence ){
    
    assignment as;
    as.score = 0;
    as.a = new int[n];
    for(int z = 0; z < n; z++){
      int m = currentPermutation[ z ];
      as.a[z] = m;
      if( maxToMin )
	as.score += C[m][z];
      else
	as.score -= ( C[m][z] + offset );
    }
    pq.push(as);
    nAssignments++;

    // Generate the next permutation sequence
    for(int m = n - 2; m >= -1; m--){

      if( m == -1){
	lastPermutationSequence = true;
	break;
      }
      
      // Find the highest index m such that currentPermutation[m] < currentPermutation[m+1]
      if(currentPermutation[m] < currentPermutation[ m + 1 ]){

	// Find highest index i such that currentPermutation[i] > currentPermutation[m] 
	// then swap the elements
	for(int i = n - 1; i >= 0; i--){
	  if( currentPermutation[i] > currentPermutation[m] ){
	    int temp = currentPermutation[i];
	    currentPermutation[i] = currentPermutation[m];
	    currentPermutation[m] = temp;
	    break;
	  }
	}

	// reverse order of elements after currentPermutation[m]
	int nElementsToSwap = n - (m + 1);
	int elementsToSwapMidPt = nElementsToSwap / 2;
	int idx1 = m + 1;
	int idx2 = n - 1;
	for(int i = 1; i <= elementsToSwapMidPt; i++){
	  int temp = currentPermutation[idx1];
	  currentPermutation[idx1] = currentPermutation[idx2];
	  currentPermutation[idx2] = temp;
	  idx1++;
	  idx2--;
	}

	break;
      }

    }

    // now we should have the next permutation sequence
  }


  if(a_ != NULL ){
    for( unsigned int m = 0; m < nAssignments; m++ ){
      delete[] a_[m];
    }
    delete[] a_;
  }
  a_ = new int* [nAssignments];
  if(s_ != NULL ){
    delete[] s_;
  }
  s_ = new double [nAssignments];

  for(int m = 0; m < nAssignments; m++){
    assignment as = pq.top();
    pq.pop();
    a_[m] = as.a;
    s_[m] = as.score;
  }
     
  a = a_;
  s = s_;
  return nAssignments;
}
