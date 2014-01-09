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

#include <cstddef>
#include <vector>

/**
 * \class CostMatrix
 * \brief Cost / Reward matrix for linear assignment
 * \author Keith Leung
 */
class CostMatrix{

public:
  
  /**
   * Constructor
   * \param[in] C A square cost matrix. Caller must take care of deallocation.
   * \param[in] nDim size of C
   */
  CostMatrix(double** C, int nDim);

  /**
   * Destructor 
   */
  ~CostMatrix();

  /**
   * Partition the cost matrix such that the rows and columns can be re-ordered to get a block diagonal matrix,
   * where each block is a partition.
   * \return the number of partitions
   */
  size_t partition();

  /**
   * Reduce the cost matrix based on assignments that are almost certain. C will be modified.
   * \param[in] lim the limit for elements in C
   * \param[in] minVal true if the limit is a lower limit for elements in C as opposed to an upper limit
   */
  void reduce(double lim, bool minVal = true);

  /**
   * Get the cost matrix
   * \param[out] C cost matrix pointer
   * \return size of C
   */
  int getCostMatrix(double** &C);

  /**
   * Get the size of a partition
   * \param[in] p partition number 
   * \param[out] nRows number of rows
   * \param[out] nCols number of cols
   */
  void getPartitionSize(int p, size_t& nRows, size_t& nCols);
  
  /**
   * Get a partition of the cost matrix
   * \param[in] p partition number 
   * \param[out] Cp partition p of the cost matrix C. Memory needs to be allocated by caller. Use getPartitionSize() to determine the size of the partition for memory allocation.
   * \param[out] isZeroPartition if not equal to NULL, indicates if partition p contains all zeros
   * \param[out] if not NULL, row_indices the original indices of the rows in partition p
   * \param[out] if not NULL, col_indices the original indices of the columns in partition p
   */
  void getPartition(int p, double** Cp, bool* isZeroPartition = NULL, int* row_indices = NULL, int* col_indices = NULL);
  

  /**
   * Get the reduced cost matrix
   * \param[out] C reduced cost matrix pointer
   * \param[out] fixedAssignments pointer to an array of size n = getCostMatrix(C_full) assignments, where -1 indicates unassgined.
   * The caller must allocate memory for this array before calling the function.
   * \param[out] score pointer to the score of the fixed assignments. Memory needs to be allocated by caller.
   * \param[out] i_remap pointer to an array of the mapping of row indices from the reduced cost matrix to the full cost matrix 
   * \param[out] j_remap pointer to an array of the mapping of col indices from the reduced cost matrix to the full cost matrix
   * \return size of C
   */
  int getCostMatrixReduced(double** &C, int* &fixedAssignments, double* score, int* &i_remap, int* &j_remap);  

private:

  double** C_; /**< cost matrix */
  double** C_reduced_; /**< reduced cost matrix with fixed assignments */
  std::vector<int> a_fixed_; /**< fixed assignments */
  std::vector<int> a_fixed_reverse_; /**< reverse fixed assignments */
  int n_; /**< size of C */
  int n_reduced_; /**< size of C_reduced */
  std::vector<int> i_reduced_; /**< Index remapping for the reduced matrix */
  std::vector<int> j_reduced_; /**< Index remapping for the reduced matrix */
  int* i_remap_; /**< Index remapping for the reduced matrix */
  int* j_remap_; /**< Index remapping for the reduced matrix */
  bool reducedCostMatAvailable_; /**< Flag indicating if the reduced cost matrix has been calculated */

  
  std::vector<unsigned int> *components_i;
  std::vector<unsigned int> *components_j;
  int combinedZeroPartition_;

};

#endif
