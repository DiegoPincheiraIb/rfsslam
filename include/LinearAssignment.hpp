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

#ifndef LINEAR_ASSIGNMENT
#define LINEAR_ASSIGNMENT

#include <cfloat>
#include <queue> // std::queue
#include <cmath>
#include <stdio.h>
#include "Tree.hpp"

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
   * \param[in] C costMat square score / cost matrix.
   * \param[in] n size of cost matrix
   * \param[out] soln assignment solution
   * \param[out] cost assignment solution cost
   * \param[in] maximize true if we want to find maximum score, false for minimum score
   * \param[in] debug creates debug printouts if true
   * \return whether a solution has been found
   */
  bool run(double** C, int n, int* soln, double* cost, bool maximize = true, bool debug = false );


};

HungarianMethod::HungarianMethod(){}

HungarianMethod::~HungarianMethod(){}

bool HungarianMethod::run(double** C, int n, int* soln, double* cost, bool maximize, bool debug){

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

  double offset = 0;
  if(!maximize){
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
  }else{
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	if(C[x][y] < offset){
	  offset = C[x][y];
	}
      }
    }
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	C[x][y] -= offset;
      }
    }
  }
  
  if(debug){
    printf("Score / Cost Matrix:\n");
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	printf("%f   ", C[x][y]);
      }
      printf("\n");
    }
  }
  

  // Step 1 - initial labeling and matching M in equality graph El
  for( int x = 0; x < n; x++){

    lx[x] = 0.0;
    ly[x] = 0.0;
    for( int y = 0; y < n; y++){
      if( C[x][y] >= lx[x] ){
	lx[x] = C[x][y];
	xy[x] = y; // temporary matching of x to y, but there may be conflicts
      }
    }

    // check if y = xy[x] is matched to another x_t = yx[y] = yx[xy[x]]
    int y = xy[x];
    x_t = yx[y];
    //printf("Trying to match x[%d] to y[%d]\n", x, xy[x]); 
    if(yx[y] != -1){ // matched with another x, conflict! Keep the one with the best score
      //printf("  y[%d] is already matched with x[%d]\n", y, yx[y]);
      if( C[x][y] > C[x_t][y] ){
	xy[x_t] = -1; // current match is better than other match, remove other match
	yx[y] = x;
	//printf("  x[%d] is no longer matched with y[%d]\n", x_t, xy[x]);
	//printf("  x[%d] is now matched with y[%d]\n", x, xy[x]);
      }else{
	xy[x] = -1; // current match is worse than other match, keep the other match
	//printf("  x[%d] cannot be matched\n", x);
      }
    }else{
      yx[y] = x;
      //printf("  x[%d] is now matched with y[%d]\n", x, xy[x]);
    }

  }

  if(debug){
    printf("Initial labeling:\n");
    for( int x = 0, y = 0; x < n; x++, y++){
      printf("  lx[%d] = %f     ly[%d] = %f\n", x, lx[x], y, ly[y]);
    }

    printf("Initial matching:\n");
    for( int x = 0; x < n; x++ )
      if( xy[x] != -1)
	printf("  x[%d] ----- y[%d]\n", x, xy[x]);
      else
	printf("  x[%d] -----      \n", x);
    for( int y = 0; y < n; y++ )
      if( yx[y] != -1)
	printf("  y[%d] ----- x[%d]\n", y, yx[y]);
      else
	printf("  y[%d] -----      \n", y);

    printf("Initial Slack Table:\n");
    for(int x = 0; x < n; x++){
      for(int y = 0; y < n; y++){
	printf("%e   ", lx[x] + ly[y] - C[x][y]);
      }
      printf("\n");
    }
  
  }

  while(true){

    // Step 2 - Check if we have perfect matching, if not, pick a new root
    
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

	// undo offset if necessary
	if(offset != 0){
	  for(x = 0; x < n; x++){
	    for(y = 0; y < n; y++){
	      C[x][y] = C[x][y] + offset;
	    }
	  }
	}

	// undo negative if necessary
	if(!maximize){
	  for(x = 0; x < n; x++){
	    for(y = 0; y < n; y++){
	      C[x][y] *= -1;
	    }
	  }
	}

	*cost = 0;
	for(x = 0; x < n; x++){
	  soln[x] = xy[x];
	  *cost += C[x][ xy[x] ];
	}
	/*printf("Solution found with score = %f\n", *cost);
        for( int x = 0; x < n; x++ )
	printf("  x[%d] ----- y[%d]\n", x, xy[x]);*/

	// Cleanup
	delete[] lx;
	delete[] ly;
	delete[] xy;
	delete[] yx;
	delete[] S;
	delete[] T;
	delete[] NS;
	delete[] slack; 
	delete[] x_q;
	delete[] y_q;
	delete[] p;

	return true;
      }
      root = x; // root of alternating tree
      S[x] = true; // put the free vertex in S
      if(debug)
	printf("Slack:\n");
      for(y = 0; y < n; y++ ){
	slack[y] = lx[x] + ly[y] - C[x][y];
	if(debug)
	  printf("    s[%d] = %f\n", y, slack[y]);
	if(  fabs(slack[y]) < 1e-14 ){
	  slack[y] = 0;
	  NS[y] = true;
	}
      }

      if(debug){
	printf("New root selected\n");  
	printf("S = { ");
	for( int x = 0; x < n; x++ ){
	  if(S[x])
	    printf("x[%d] ", x);
	}
	printf(" }\n");
	printf("T = { ");
	for( int y = 0; y < n; y++ ){
	  if(T[y])
	    printf("y[%d] ", y);
	}
	printf(" }\n");
	printf("NS = { ");
	for( int y = 0; y < n; y++ ){
	  if(NS[y])
	    printf("y[%d] ", y);
	}
	printf(" }\n");
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


      if(debug){
	printf("Slack Table:\n");
	for(int x = 0; x < n; x++){
	  for(int y = 0; y < n; y++){
	    printf("%e   ", lx[x] + ly[y] - C[x][y]);
	  }
	  printf("\n");
	}
      }

      double a = DBL_MAX;
      //calculate delta using slack    
      for (y = 0; y < n; y++){   	     
	if (!T[y]){
	  a = fmin(a, slack[y]);
	}
      }

      if(debug){
	printf("Update labels:\n");
	printf("  minimum slack of !T = %f\n",a);
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
	  slack[y] -= a;
	}
	if(slack[y] == 0){
	  NS[y] = true;
	}
      }

      if(debug){
	for( int x = 0, y = 0; x < n; x++, y++)
	printf("  lx[%d] = %f     ly[%d] = %f\n", x, lx[x], y, ly[y]);
	printf("  Update slack variables:\n");
	for( int y = 0; y < n; y++ )
	  printf("    s[%d] = %f\n", y, slack[y]);

	printf("Slack Table:\n");
	for(int x = 0; x < n; x++){
	  for(int y = 0; y < n; y++){
	    printf("%e   ", lx[x] + ly[y] - C[x][y]);
	  }
	  printf("\n");
	}
      }

    }//end updateLabel

    // Step 4 - N(S) != T, select y in N(S) but not in T
    //printf("Finding y in NS and not T\n");
    for(y = 0; y < n; y++ ){
      if( NS[y] && !T[y] )
	break;
    }
    x_t = yx[y]; // match to y
    if(x_t == -1){ // if y is free, path root -> y is an augmenting path, goto step 2

      bool augmentingPathFound = false;
      
      if(debug){
	printf("  Found free y[%d]\n", y);
      }

      // Breadth first search for path root -> y
      // We need to make distinct indices for vertices in X and in Y
      // So in this search only, indices for X are 0, 1, 2, 3 ... n-1
      // and indices for Y are n+0, n+1, n+2 ... 2n-1 when storing indices in the queue
      int target = y + n;
      std::queue<int> q_empty;
      std::swap(q, q_empty); // clear the queue
      q.push(root);
      for(x = 0; x < n; x++){
	x_q[x] = false;
	y_q[x] = false;
      }
      x_q[root] = true;
      for(x = 0; x < 2*n; x++){
	p[x] = -1;
      }
      int t = q.front();
      while(!q.empty()){
 
	t = q.front();
	if(debug){
	  if( t < n){
	    printf("Breadth first search queue front: x[%d]\n",t);
	  }else{
	    printf("Breadth first search queue front: y[%d]\n",t-n);	  
	  }
	}

	if (t == target){
	  // path found, perform augmentation in inverse the matching along this path
	  if(debug){
	    printf("  Augmentation path found:\n    y[%d]", t-n);
	  }
	  while(t != root){
	    if( t >= n){ // t is in Y, p[t]--t should be a match
	      x_t = p[t];
	      xy[x_t] = t-n;
	      yx[t-n] = x_t;
	    }
	    t = p[t];

	    if(debug){
	      if( t < n ){
		printf(" --- x[%d]", t);
	      }else{
		printf(" --- y[%d]", t-n);
	      }
	    }
	  }

	  if(debug){
	    printf("\nCurrent pairing\n");
	  
	    for( int x = 0; x < n; x++ )
	      if( xy[x] != -1){
		printf("  x[%d] ----- y[%d]\n", x, xy[x]);
	      }else{
		printf("  x[%d] -----      \n", x);
	      }
	    for( int y = 0; y < n; y++ )
	      if( yx[y] != -1){
		printf("  y[%d] ----- x[%d]\n", y, yx[y]);
	      }else{
		printf("  y[%d] -----      \n", y);
	      }
	  }
	  augmentingPathFound = true;
	  break;
	}
	q.pop();

	if( t < n){ // t is a vertex in X 
	  for(y = 0; y < n; y++){ // look through children of t
	    // check if y is in equality graph and not already queued and unmatched
	    if(fabs(lx[t] + ly[y] - C[t][y]) < 1e-12 && !y_q[y] && xy[t] != y){
	      if(debug){
		printf("Breadth first search inserting y[%d] in queue with parent x[%d]\n", y, t); 
	      }
	      y_q[y] = true;
	      p[y+n] = t;
	      q.push( y+n );
	    }
	  }
	}else{ // t >= n and is a vertex in Y
	  t -= n; // fix the index back to original
	  for(x = 0; x < n; x++){ // look through children of t
	    // check if x is in equality graph and in set S and not already queued and matched
	    if(fabs(lx[x] + ly[t] - C[x][t]) < 1e-12 && S[x] && !x_q[x] && yx[t] == x){
	      if(debug){
		printf("Breadth first search inserting x[%d] in queue with parent y[%d]\n", x, t);
	      }
	      x_q[x] = true;
	      p[x] = t+n;
	      q.push( x );
	    }
	  }
	}

      } // end while

      if(!augmentingPathFound){
	printf("Cannot find alternating path\n");
	/*printf("Slack table:\n");
	for(int x = 0; x < n; x++){
	  for(int y = 0; y < n; y++){
	    printf("%e   ", lx[x] + ly[y] - C[x][y]);
	  }
	  printf("\n");
	  }*/
	return false;
      }

      pickFreeVertex = true; // go back to step 2

    }else{ // y is matched to x_t, add to alternating tree, go to step 3

      if(debug){
	printf("  Found y[%d] matched with x[%d]\n", y, x_t);
	printf("  Examine slack for x[%d]:\n", x_t);
      }
      S[x_t] = true;
      T[y] = true;
      for(int y_t = 0; y_t < n; y_t++ ){
	if(debug)
	  printf("    s[%d -- %d] = %f\n", x_t, y_t, lx[x_t] + ly[y_t] - C[x_t][y_t]);
	if( fabs(lx[x_t] + ly[y_t] - C[x_t][y_t]) < 1e-14 )
	  NS[y_t] = true;
      }
      
      if(debug){
	printf("  S = { ");
	for( int x = 0; x < n; x++ ){
	  if(S[x]){
	    printf("x[%d] ", x);
	  }
	}
	printf(" }\n");
	printf("  T = { ");
	for( int y = 0; y < n; y++ ){
	  if(T[y]){
	    printf("y[%d] ", y);
	  }
	}
	printf(" }\n");
	printf("  NS = { ");
	for( int y = 0; y < n; y++ ){
	  if(NS[y]){
	    printf("y[%d] ", y);
	  }
	}
	printf(" }\n");
      }

      // It is necessary to update the slack variables here since S now includes x_t
      for(int y = 0; y < n; y++ ){
	double slack_x_t = lx[x_t] + ly[y] - C[x_t][y]; 
	if( slack_x_t < slack[y]){
	  slack[y] = slack_x_t;
	}
      }

      if(debug){
	printf("  Update slack variables for !T:\n");
	for( int y = 0; y < n; y++ )
	  if(!T[y])
	    printf("    s[%d] = %f\n", y, slack[y]);
      }      

      pickFreeVertex = false; // go back to step 3

    }

  }

}

#endif

/**
 * \class MurtyNode
 * A node used in the Murty class for k-best linear assignment
 * \brief A node used in the Murty class for k-best linear assignment
 */
class MurtyNode : public Node
{
public:

  /** Constructor 
   * \param[in] id The partition id number for this node. Default is -1
   * \param[in] n_children_exp the number of children this node is expected to have,
   * used for more efficient memory allocation. Not mandatory.
   */
  MurtyNode(int id = -1, int n_children_exp = -1) : Node(n_children_exp),
						    id_(id),
						    assignment_(NULL), 
						    score_(0){}

  /** Destructor */
  ~MurtyNode(){
    if( assignment_ != NULL )
      delete[] assignment_;
    assignment_ = NULL;
  }

  /** Comparison operator for sorting */
  bool operator< (const MurtyNode& other){
    if( score_ < other.score_ )
      return true;
    return false;
  }

  /**
   * Set the assignment and score presented by this node 
   * \param[in] a assignment, where b = a[i] implies b is assigned to a_i. This
   * will be deallocated when this node is destroyed.
   * \param[in] s total score from all assignments;
   */
  void setAssignment(int* a, double s){
    if(assignment_ != NULL)
      delete[] assignment_;
    assignment_ = a;
    score_ = s;
  }

  /**
   * Get the score of this node
   * \return score
   */
  double getScore(){
    return score_;
  }

  /**
   * Get the assignment represented by this node 
   * \return pointer to assignment array 
   */
  int* getAssignment(){
    return assignment_;
  }

  /** 
   * Get the ID number of the node 
   */
  int getPartitionID(){
    return id_;
  }


private:

  int id_; 
  int* assignment_;
  double score_;

};

/**
 * \class MurtyNodeCompare
 * A comparison class used for the priorty queue in Murty
 * \brief A comparison class used for the priorty queue in Murty
 */
class MurtyNodeCompare{
public:
  bool operator()(MurtyNode* p1, MurtyNode* p2){
    if( p1->getScore() < p2->getScore() )
      return true;
    return false;
  }
};

/**
 * \class Murty
 * Murty's k-best linear assignment algorithm
 * Given the partitioning of a linear assignment problem, and the
 * current best possible assignment,  Murty's algorithm will find 
 * the next best assignment
 * \brief Murty's k-best linear assignment algorithm
 */ 
class Murty{

public:
   
  /** constructor 
   *  \param[in] C square score matrix from which we will find the best assignments
   *  \param[in] n dimension of C
   *  \param[in] bigNum A number whose negative is used to replace entries in the cost matrix that cannot be used for assignment
   */ 
  Murty(double** C, int n, double bigNum = 10000) : k_(0), n_(n), bigNumber_(bigNum){
    C_ = C;
    C_t_ = new double* [n_];
    for( int i = 0; i < n_; i++){
      C_t_[i] = new double [n_];
    }
    root_ = new MurtyNode(0, n_);
    bestScore_ = 0;
  }
  
  /** destructor */
  ~Murty(){
    for(int i = 0; i < n_; i++)
      delete[] C_t_[i];
    delete[] C_t_;
    delete root_;
  }

  /**
   * Get the best score.
   * \note The best score = 0 when object is first initiated.
   * \return the best score
   */
  double getBestScore(){
    return bestScore_;
  }

  /** Find the next best assignment 
   *  \param[out] assignment pointer to array of assignment (no need to allocate/deallocate )
   *  \param[out] score the assignment score 
   *  \return k, the k-best solution index, or -1 if no more solutions are available
   */
  int findNextBest( int* &assignment, double* score){
    
    if( k_ == 0 ){
      int* a = new int[n_];
      double s; 
      lam_.run(C_, n_, a, &s);
      root_->setAssignment(a, s);
      k_++;
      pq.push(root_);
      assignment = a;
      *score = s;
      bestScore_ = s;
      return k_;
    }

    // pop the highest-score node
    if(pq.empty()){
      //printf("No more solutions available\n");
      assignment = NULL;
      *score = 0;
      return -1;
    }
    MurtyNode* parent = pq.top();
    int parent_partition = parent->getPartitionID();
    pq.pop();

    // partition this node (into n_ - 1 - parent_partition parts) -- make children nodes
    //printf("\n************************\n");
    for(int n = parent_partition; n < n_-1; n++ ){ 
      
      //printf("\nPartition %d\n", n);
      MurtyNode* p = new MurtyNode(n, n_);
      parent->addChild(p);
      int* a = new int[n_];
      bool freeCol[n_];
      for(int i = 0; i < n_; i++){
	freeCol[i] = true;
      }
      double assignmentFixedScore = 0;

      // Copy all the fixed assignments from the parent node
      //printf("Must have: ");
      int* a_parent = parent->getAssignment();
      for(int i = 0; i < parent_partition; i++){
	a[i] = a_parent[i];
	freeCol[a[i]] = false;
	assignmentFixedScore += C_[i][a[i]];
	//printf("(%d,%d)[%f] ", i, a[i], C_[i][a[i]]);
      }

      // assignments that have to be included in this partition
      for(int i = parent_partition; i < n; i++){
	a[i] = a_parent[i];
	freeCol[a[i]] = false;
	assignmentFixedScore += C_[i][a[i]];
	//printf("(%d,%d)[%f] ", i, a[i], C_[i][a[i]]);
      } 
      //printf("\n");

      // build up the cost matrix to solve, which may be smaller than the original cost matrix due to fixed assignments
      int nFreeAssignments = n_ - n;
      int rowRemap[nFreeAssignments]; // find entry in original matrix from reduced matrix
      int rowRemapR[n_]; // find extry in reduced matrix from original matrix
      int colRemap[nFreeAssignments];
      int colRemapR[n_];
      int nFreeCols = 0;

      for(int i = 0; i < nFreeAssignments; i++){
	rowRemap[i] = n + i;
	rowRemapR[n+i] = i;
	//printf("Remap row %d to %d\n", rowRemap[i], i);
      }
      for(int j = 0; j < n_; j++){
	if(freeCol[j]){
	  colRemap[nFreeCols] = j;
	  colRemapR[j] = nFreeCols;
	  //printf("Remap col %d to %d\n", colRemap[nFreeCols], nFreeCols);
	  nFreeCols++; 
	}
      } // at this point nFreeCol should equal nFreeAssignments
      for(int i = 0; i < nFreeAssignments; i++ ){
	
	for(int j = 0; j < nFreeAssignments; j++ ){
	  //printf("Remap [%d,%d] to [%d,%d]\n", original_i, original_j, i, j);
	  C_t_[i][j] = C_[rowRemap[i]][colRemap[j]];
	}
      }

      // Set all the assignments that cannot occur
      //printf("Must not have: ");
      MurtyNode* current = p;
      MurtyNode* next;
      do{
	int currentPart = current->getPartitionID();
	if( currentPart < p->getPartitionID() ){
	  break;
	}
	next = static_cast< MurtyNode* > ( current->getParent() );
	int* nextAssignment = next->getAssignment();
	int doNotAssign_i = rowRemapR[currentPart];
	int doNotAssign_j = colRemapR[nextAssignment[currentPart]];
	C_t_[doNotAssign_i][doNotAssign_j] = -bigNumber_;	
	//printf("(%d,%d)=>(%d,%d) ", currentPart, nextAssignment[currentPart],doNotAssign_i,doNotAssign_j);
	current = next;
	
      }while (current != root_ );
      //printf("\n");
      bool solutionPossible = false;
      int constraintRow = rowRemapR[p->getPartitionID()];
      for(int j = 0; j < nFreeAssignments; j++){
	if( C_t_[constraintRow][j] != -bigNumber_){
	  solutionPossible = true;
	  break;
	}
      }

      // Solve best-linear assignment for the n partitions
      if(solutionPossible){
	int aTmp[nFreeAssignments];
	double s = 0;
	/*printf("Solving best linear assignment for cost matrix:\n");
	for(int i = 0; i < nFreeAssignments; i++ ){
	  for(int j = 0; j < nFreeAssignments; j++ ){
	    printf("%f   ", C_t_[i][j]);
	  }
	  printf("\n");
	  }*/
	lam_.run(C_t_, nFreeAssignments, aTmp, &s);
      
	double s_more_accurate = 0;
	for(int i = 0; i < nFreeAssignments; i++){
	  int i_actual = rowRemap[i];
	  int j_actual = colRemap[aTmp[i]];
	  a[i_actual] = j_actual;
	  s_more_accurate += C_[i_actual][a[i_actual]];
	}
	/*printf("Assignment:\n");
	for(int i = 0; i < n_; i++){
	  printf("x[%d] ----- y[%d]\n", i, a[i]);
	  }*/
	s_more_accurate += assignmentFixedScore;
	// printf("Score: %f\n", s);
      
	// Insert the solution in the priority queue and set it as a child of the highest-score node
	p->setAssignment(a, s_more_accurate);
	pq.push(p);
      }else{
	//printf("No solution possible with this partition\n");
	for(int i = 0; i < n_; i++){
	  a[i] = -1;
	}
	p->setAssignment(a, -1);
      }
      
    } // end for all partitions

    // All new partitions have been inserted, now look for the highest score
    if(pq.empty()){
      //printf("No more solutions available\n");
      assignment = NULL;
      *score = 0;
      return -1;
    }
    MurtyNode* node_hi = pq.top();
    assignment = node_hi->getAssignment();
    *score = node_hi->getScore();
    
    k_++;

    return k_;
  }

private:
  
  MurtyNode* root_; /**< The best linear assignment */
  int k_; /**< number of solutions found */
  double** C_; /**< the nxn score matrix */
  double** C_t_; /**< temporary score matrix */
  int n_; /**< dimension of C_ and C_t_*/
  HungarianMethod lam_; /**< Linear assignment method */
  std::priority_queue<MurtyNode*, std::vector<MurtyNode*>, MurtyNodeCompare> pq; /**< Priority queue for easy retrieval of node with highest score */
  double bestScore_;
  double bigNumber_;

};

/**
 * \class BruteForceLinearAssignment
 * This function finds all the linear assignments and orders them from best to worst.
 * It is intended for testing and checking the Hungarian method and Murty's k-best algorithm.
 * \brief Brute-force linear assignment 
 */
class BruteForceLinearAssignment
{

public:

  BruteForceLinearAssignment() : a_(NULL), s_(NULL), nAssignments(0){}
  ~BruteForceLinearAssignment(){
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

  /**
   *  Run brute-force linear assignment by finding the cost of all possible linear assignments
   *  in a lexicographical order, and then returning the results in an ordered format with 
   *  respect to the cost / score.
   *  \param[in] C square cost / score matrix
   *  \param[in] n dimension of C
   *  \param[out] a ordered assignments a[k][n], where k is the assignment order
   *  \param[out] s ordered assignment scores
   *  \param[in] maxToMin ordering of results
   *  \return number of assignments
   */
   unsigned int run(double** C, int n, int** &a, double* &s, bool maxToMin = true){

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

private:

  /** \brief A linear assignment */
  struct assignment{
    int* a;
    double score;
    
    bool operator<(const assignment& rhs) const{
      if(score < rhs.score)
	return true;
      return false;
    }

  };

  std::priority_queue<assignment> pq;
  unsigned int nAssignments;
  int** a_;
  double* s_;
};
