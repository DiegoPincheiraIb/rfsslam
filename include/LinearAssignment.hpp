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
 *     * Neither the name of the owner nor the names of its contributors may be 
 *       used to endorse or promote products derived from this software without 
 *       specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef LINEAR_ASSIGNMENT
#define LINEAR_ASSIGNMENT

#include <queue> // std::queue
#include <cfloat>

/** 
 * \class HungarianMethod
 * This class is an implementation of the Hungarian method for a linear assignment
 * problem specified by a NxN cost matrix. It has complexity O(N^3).
 * This code is referenced from the Top Coder tutorial on the Hungarian method: 
 * http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=hungarianAlgorithm 
 * The following pdf is also a good reference for the method:
 * www.cse.ust.hk/~golin/COMP572/Notes/Matching.pdf
 * \brief The Hungarian Method for linear assignmnet
 * \author Keith Leung
 */
class HungarianMethod
{
public:

  /** Default constructor */
  HungarianMethod();

  /** Default destructor */
  ~HungarianMethod();

  /**
   * Run the Hungarian method
   * \param[in] C square cost matrix
   * \param[in] n size of cost matrix
   * \param[out] soln assignment solution
   * \param[out] cost assignment solution cost
   */
  run(double** C, int n, int* soln, double* cost );

private:

  int* lx; // label for set X
  int* ly; // label for set Y
  int* xy; // vertex of Y matched to x
  int* yx; // vertex of X matched to y
  bool* S[N]; // true if X[i] is in set S
  bool* T[N]; // true if Y[i] is in set T
  double* slack[N]; // slack variable
  double* slackx[N]; // slack variable
  int* prev[N]; // alternating path

};

HungarianMethod::HungarianMethod(){}

HungarianMethod::~HungarianMethod(){}

HungarianMethod::run(double** C, int n, int* soln, double* cost){

  // Initialize
  double* lx = new double[n]; // label for X
  double* ly = new double[n]; // label for Y
  int* xy = new int[n]; // match to x = X[n] 
  int* yx = new int[n]; // match to y = Y[n]
  bool* S = new bool[n];
  bool* T = new bool[n];
  bool* NS = new bool[n]; // neighbourhood of S, a subset or equal to Y
  double* slack = new double[n];
  
  std::queue<int> q; // queue for breadth-first-search (bfs)
  bool* x_q = new bool[n]; // flag indicating x has already been queued in bfs
  bool* y_q = new bool[n]; // flag indicating y has already been queued in bfs
  int* p = new int[2*n]; // parent of vertices in bfs

  bool** E = new bool*[n]; // Equality graph (may not need)
  for(int x = 0; x < n; x++){
    E[x] = new bool[n];
  }
  
  int x, x_t, y; // indices
  int root; // root index
  bool pickFreeVertex = true; // flag to pick a new root, step 2
  bool updateLabel; // flag to perform update label, step 3

  for(int x = 0; x < n; x++){
    xy[x] = -1;
    S[x] = false;
  }
  for(int y = 0; y < n; y++){
    yx[y] = -1;
    T[y] = false;
  }

  // Step 1 - initial labeling and matching M in equality graph El
  for( int x = 0; x < n; x++){

    lx[x] = 0;
    ly[x] = 0;
    for( int y = 0; y < n; y++){
      if( C[x][y] > lx[x] ){
	lx[x] = C[x][y];
	xy[x] = y; // temporary matching of x to y, but there may be conflicts
      }
    }
        
    // check of y = xy[x] is matched to another x_t = yx[y] = yx[xy[x]]
    y = xy[x];
    x_t = yx[y];
    if(yx[y] != -1){ // matched with another x, conflict! Keep the one with the best score
      if( C[x][y] > C[x_t][y] ){
	xy[x_t] = -1; // current match is better than other match, remove other match
	yx[y] = x;
      }else{
	xy[x] = -1; // current match is worse than other match, keep the other match
      }
    }    

  }

  while(true){

    // Step 2 - Check if we have perfect matching
    if(pickFreeVertex){
    
      for(x = 0; x < n; x++){
	S[x] = false;
      }
      for(y = 0; y < n; y++){
	T[y] = false;
	NS[y] = false;
      }
      for(x = 0; x < n; x++){
	if(xy[x] == -1)
	  break; // x is not matched and therefore a free vertex
      }
      if(x == n){ // no more free / unmatched x, we have found a solution
	cost = 0;
	for(x = 0; x < n; x++){
	  soln[x] == xy[x];
	  cost += C[x][xy[x]];
	}
	return;
      }
      root = x; // root of alternating tree
      S[x] = true; // put the free vertex in S
      for(y = 0; y < n; y++ ){
	if( E[x][y] == true )
	  NS[y] = true;
      }

    }

    // Step 3 - Update labels (if neighbourhood N(S) == T)
    updateLabel = true;
    for(y = 0; y < n; y++ ){
      if( NS[y] != T[y] ){
	updateLabel = false;
	break;
      }
    }
    if(updateLabel){ // After updating label, N(S) != T

      double a = DBL_MAX;
      //calculate delta using slack    
      for (y = 0; y < n; y++){   	     
	if (!T[y]){
	  a = min(a, slack[y]);
	}
      }
    
      //update X labels
      for (x = 0; x < n; x++){            
	if (S[x]){ 
	  lx[x] -= a;
	}
      }
    
      //update Y labels
      for (y = 0; y < n; y++){            
	if (T[y]){ 
	  ly[y] += a;
	}
      }

      //update slack
      for (y = 0; y < n; y++){            
	if (!T[y]){
	  slack[y] -= delta;
	}
      }

    }//end updateLabel

    // Step 4 - N(S) != T, select y in N(S) but not in T
    for(y = 0; y < n; y++ ){
      if( NS[y] && !Y[y] )
	break;
    }
    x_t = yx[y]; // match to y
    if(x_t == -1){ // if y is free, path root -> y is an augmenting path, goto step 2

      // Breadth first search for path root -> y
      // We need to make distinct indices for vertices in X and in Y
      // So in this search only, indices for X are 0, 1, 2, 3 ... n-1
      // and indices for Y are n+0, n+1, n+2 ... 2n-1 when storing indices in the queue
      int target = y + n;
      q.clear();
      q.push(root);
      for(x = 0; x < n; x++){
	x_q[x] = false;
	y_q[x] = false;
      }
      x_q[root] = true;
      for(x = 0; x < 2*n; x++){
	p[x] = -1;
      }
      while(!q.empty()){
      
	int t = q.front();
	if (t == target){
	  // path found, perform augmentation in inverse the matching along this path
	  while(t != root){
	    if( t >= n){ // t is in Y, p[t]--t should be a match
	      x_t = p[t - n];
	      xy[x_t] = t-n;
	      yx[t-n] = x_t;
	    }
	    t = p[t];
	  }
	  break;
	}
	q.pop();

	if( t < n){ // t is a vertex in X 
	  for(y = 0; y < n; y++){
	    // check if y is in equality graph and in set T and not already queued
	    if(E[t][y] && T[y] && !y_q[y]){
	      y_q[y] = true;
	      y_p[y] = t;
	      q.push( y+n );
	    }
	  }
	}else{ // t >= n and is a vertex in Y
	  t -= n; // fix the index back to original
	  for(x = 0; x < n; x++){
	    // check if x is in equality graph and in set S and not already queued
	    if(E[x][t] && S[x] && !x_q[x]){
	      x_q[x] = true;
	      x_p[x] = t+n;
	      q.push( x );
	    }
	  }
	}

      } // end while

      pickFreeVertex = true; // go back to step 2

    }else{ // y is matched to x_t, add to alternating tree, go to step 3
      S[x_t] = true;
      for(y = 0; y < n; y++ ){
	if( E[x_t][y] == true )
	  NS[y] = true;
      }
      T[y] = true;

      // It is necessary to update the slack variables here since S now includes x_t
      for(y = 0; y < n; y++ ){
	double slack_x_t = lx[x_t] + ly[y] - C[x_t][y]; 
	if( slack_x_t < slack[y]){
	  slack[y] = slack_x_t;
	}
      }

      pickFreeVertex = false; // go back to step 3
    }

  }

  // Cleanup
  delete[] lx;
  delete[] ly;
  delete[] xy;
  delete[] yx;
  delete[] S;
  delete[] T;
  delete[] NS;
  delete[] slack;
  delete[] slackx;
  delete[] prev;  
  for(int x = 0; x < n; x++)
    delete[] E[n];
  delete[] E;
  delete[] x_q;
  delete[] y_q;
  delete[] p
}

#endif
