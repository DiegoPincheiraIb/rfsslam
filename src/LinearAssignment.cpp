#include <assert.h>
#include "LinearAssignment.hpp"
#include <stdio.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>


/// Cost Matrix Implementation

CostMatrix::CostMatrix(double** C, int nDim) : C_(C), 
					       n_(nDim), 
					       C_reduced_(NULL), 
					       n_reduced_(0),
					       reducedCostMatAvailable_(false),
					       i_remap_(NULL),
					       j_remap_(NULL),
					       components_i(NULL),
					       components_j(NULL),
					       combinedZeroPartition_(-1)
					       
{

  a_fixed_.resize(nDim, -1);
  a_fixed_reverse_.resize(nDim, -1);
}

CostMatrix::~CostMatrix(){

  if( n_reduced_ > 0){
    for(int i = 0; i < n_reduced_; i++){
      delete[] C_reduced_[i];
    }
    delete[] C_reduced_;
  }
  if(i_remap_ != NULL)
    delete[] i_remap_;
  if(j_remap_ != NULL)
    delete[] j_remap_;

  if(components_i != NULL){
    delete[] components_i;
  }
  if(components_j != NULL){
    delete[] components_j;
  }

}

void CostMatrix::reduce(double lim, bool minVal){

  int nMatch_i[n_]; // number of possible matches
  int nMatch_j[n_];
  int best_score[n_];
  for(int n = 0; n < n_; n++){
    nMatch_i[n] = 0;
    nMatch_j[n] = 0;
    best_score[n] = lim;
  }

  for(int i = 0; i < n_; i++){
    for(int j = 0; j < n_; j++){

      if( (C_[i][j] >= lim && !minVal) || (C_[i][j] <= lim && minVal) ){ 
	C_[i][j] = lim;
      }else{
	nMatch_i[i]++;
	nMatch_j[j]++;

	if( nMatch_i[i] == 1 && nMatch_j[j] == 1){
	  a_fixed_[i] = j;
	  a_fixed_reverse_[j] = i;
	}

	if( nMatch_i[i] > 1){
	  a_fixed_[i] = -1;
	}
	if( nMatch_j[j] > 1){
	  a_fixed_reverse_[j] = -1;
	}
      }
      
    }
  }

  for(int n = 0; n < n_; n++){
    if( nMatch_j[ a_fixed_[n] ] != 1 ){
      a_fixed_[n] = -1;
    }
    if( nMatch_i[ a_fixed_reverse_[n] ] != 1){
      a_fixed_reverse_[n] = -1;
    }
    if( a_fixed_[n] == -1){
      i_reduced_.push_back(n);
    }
    if( a_fixed_reverse_[n] == -1){
      j_reduced_.push_back(n);
    }
  }

  assert( i_reduced_.size() == j_reduced_.size());
  n_reduced_ = i_reduced_.size();

  if(n_reduced_ != 0){

    C_reduced_ = new double* [n_reduced_];
    for(int i = 0; i < n_reduced_; i++ ){
      C_reduced_[i] = new double [n_reduced_];
      for(int j = 0; j < n_reduced_; j++ ){
	C_reduced_[i][j] = C_[ i_reduced_[i] ][ j_reduced_[j] ];
      } 
    }

    if(n_reduced_ == 1){

      a_fixed_[ i_reduced_[0] ] = j_reduced_[0];
      n_reduced_ = 0;

    }
  }

  reducedCostMatAvailable_ = true;
}

size_t CostMatrix::partition(){

  // partitioning using connected-component analysis
  boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> G;
  for(int i = 0; i < n_; i++){
    for(int j = 0; j < n_; j++){
      if (C_[i][j] != 0 ){
	boost::add_edge(i, j+n_, G);
      }
    }
  }
  boost::add_edge(n_*2-1, n_*2-1, G); // add edge to self so the graph has the desired number of vetices

  std::vector<int> cc_results( n_*2, -1 );
  int ncc = boost::connected_components(G, &cc_results[0]);

  if(components_i != NULL){
    delete[] components_i;
  }
  components_i = new std::vector<unsigned int> [ncc];
  for(int i = 0; i <  n_; i++){
    components_i[cc_results[i]].push_back(i);
  }
  if(components_j != NULL){
    delete[] components_j;
  }
  components_j = new std::vector<unsigned int> [ncc];
  for(int j = n_; j < n_*2; j++){
    components_j[cc_results[j]].push_back(j - n_);
  }

  combinedZeroPartition_ = -1;
  int nMergedZeroPartitions = 0; 
  for(int n = 0; n < ncc; n++){
    if(components_i[n].size() == 0 || components_j[n].size() == 0){
      
      if(combinedZeroPartition_ == -1){
	combinedZeroPartition_ = n;
      }else if(components_i[n].size() != 0){
	components_i[combinedZeroPartition_].push_back(components_i[n][0]);
	nMergedZeroPartitions++;
      }else{
	components_j[combinedZeroPartition_].push_back(components_j[n][0]);
	nMergedZeroPartitions++;
      }

    }
  }
  return ncc - nMergedZeroPartitions;

}

void CostMatrix::getPartitionSize(int p, size_t& nRows, size_t& nCols){
  nRows = components_i[p].size();
  nCols = components_j[p].size();
}

void CostMatrix::getPartition(int p, double** Cp, bool* isZeroPartition, int* row_indices, int* col_indices){

  if(Cp != NULL){
    for(int i = 0; i < components_i[p].size(); i++){
      for(int j = 0; j < components_j[p].size(); j++){
	Cp[i][j] = C_[ components_i[p][i] ][ components_j[p][j]];
      }  
    }
  }

  if(isZeroPartition != NULL){
    *isZeroPartition = false;
    if(p == combinedZeroPartition_)
      *isZeroPartition = true;
  }

  if(row_indices != NULL){
    for(int i = 0; i < components_i[p].size(); i++){
      row_indices[i] = components_i[p][i];
    }
  }

  if(col_indices != NULL){
     for(int j = 0; j < components_j[p].size(); j++){
       col_indices[j] = components_j[p][j];
     }
  }

}

int CostMatrix::getCostMatrix(double** &C){
  C = C_;
  return n_;
}

int CostMatrix::getCostMatrixReduced(double** &C, int* &fixedAssignments, double* fixedScore, int* &i_remap, int* &j_remap ){

  if(!reducedCostMatAvailable_)
    return -1;

  if(i_remap_ == NULL)
    i_remap_ = new int[n_];
  if(j_remap_ == NULL) 
    j_remap_ = new int[n_];

  C = C_reduced_;
  *fixedScore = 0;
  for(int n = 0; n < n_; n++){
    fixedAssignments[n] = a_fixed_[n];
    if(fixedAssignments[n] != -1){
      *fixedScore += C_[n][fixedAssignments[n]];
    }
  }

  for(int n = 0; n < n_reduced_; n++){
    i_remap_[n] = i_reduced_[n];
    j_remap_[n] = j_reduced_[n];
  }
  i_remap = i_remap_;
  j_remap = j_remap_;

  return n_reduced_;
}

