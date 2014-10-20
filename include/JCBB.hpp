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

#ifndef JCBB_HPP
#define JCBB_HPP

#include <boost/foreach.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/shared_ptr.hpp>
#include <cstddef>
#include <float.h>
#include <queue>
#include <vector>

#include "MeasurementModel.hpp"
#include "Tree.hpp"

namespace rfs{

  /**
   * \class JCBB_TreeNode
   * \brief Node for the JCBB interpretation tree
   * Node for the JCBB interpretation tree
   */
  class JCBB_TreeNode: public TreeNode
  {
  public:

    /**
     * Default constructor
     */
    JCBB_TreeNode();

    /**
     * Destructor
     */
    ~JCBB_TreeNode();

    /** 
     * Less than operator for use with priorty queue
     */
    bool operator < (JCBB_TreeNode &other) const;

    /**
     * Set the data association represented by this node
     */
    void setAssociation(size_t z_idx, size_t m_idx, double dist_mahalanobist2);

    /**
     * Get the association represented by this node
     */ 
    void getAssociation(size_t &z_idx, size_t &m_idx);

    /**
     * Get the joint Mahalanobis distance squared of all associations up to the root
     * \return Mahalanobis distance squared
     */
    double getMahalanobisDist2();

    /**
     * Get the measurement innovation of all associations up to the root
     * \return pointer to innovation vector
     */
    ::Eigen::VectorXd *getInnovation();

    /**
     * Get the measurement innovation covariance of all assoications up to the root
     * \return pointer to innovation covariance
     */
    ::Eigen::MatrixXd* getInnovationCovInv();

  private:

    size_t z_idx_;
    size_t m_idx_;
    double dist_mahalanobist2_; // joint Mahalanobis distance for all associations up to the root
    ::Eigen::MatrixXd C_inverse_; // innovation covariance for all associations up to the root
    ::Eigen::VectorXd innov_; // innovation for all associations up to the root

  };

  /**
   * \class JCBB
   * \brief Joint compatibility branch and bound data association algorithm
   * An implementation of the Joint compatibility branch and bound data association algorithm
   * Reference:  Neira, J. Tardos; Tardos, J.D. (2001), "Data association in stochastic mapping 
   * using the joint compatibility test", Robotics and Automation, IEEE Transactions on 17 (6): 
   * pp 890–897, doi:10.1109/70.976019.
   */
  template<class MeasurementModelType>
  class JCBB
  {

  public:

    typedef typename MeasurementModelType::TPose TPose;
    typedef typename MeasurementModelType::TLandmark TLandmark;
    typedef typename MeasurementModelType::TMeasurement TMeasurement;
    typedef typename MeasurementModelType::TJacobianPose TJacobianPose;
    typedef typename MeasurementModelType::TJacobianLmk TJacobianLmk;

    /**
     * Constructor
     */
    JCBB();

    /**
     * Constructor
     */
    JCBB(double confidenceItvl,
	 MeasurementModelType *measurementModel,
	 std::vector< TMeasurement > *measurements,
	 TPose const * robot,
	 std::vector< TLandmark const * > *landmarks,
	 ::Eigen::MatrixXd const *estCovDense = NULL);

    /**
     * Destructor
     */
    ~JCBB();


  private:

    double confidence_interval_;

    JCBB_TreeNode iTreeRoot_;
    
    MeasurementModelType *measurementModel_;
    std::vector<TMeasurement> *Z_actual_ptr_;
    std::vector<TMeasurement> Z_predict_;
    std::vector< boost::shared_ptr< TJacobianPose > > jacobians_wrt_pose_; 
    std::vector< boost::shared_ptr< TJacobianLmk  > > jacobians_wrt_lmk_; 

    TPose const *robot_;
    std::vector<TLandmark*> *landmarks_;
    bool isEstCovBlockDiag_;
    ::Eigen::MatrixXd estCovDense_;

    std::priority_queue<JCBB_TreeNode*> branchList_; // priority queue of nodes to branch

    int nAssociations_best;
    double dist_mahalanobist2_best;
    JCBB_TreeNode* node_best;

    /**
     * calculate expected measurements and covariance 
     */
    void calcExpectedMeasurements();
    
    /**
     * Branch a node of the interpretation tree
     */
    void branch(JCBB_TreeNode *node);

  };

} // namespace rfs


// implementation JCBB_TreeNode

namespace rfs{

  JCBB_TreeNode::JCBB_TreeNode(){

  }

  JCBB_TreeNode::~JCBB_TreeNode(){

  }

  bool JCBB_TreeNode::operator<(JCBB_TreeNode &other) const{
    if(this->dist_mahalanobist2_ < other.dist_mahalanobist2_)
      return true;
    else
      return false;
  }

  void JCBB_TreeNode::setAssociation(size_t z_idx, size_t m_idx, double dist_mahalanobist2){
    z_idx_ = z_idx;
    m_idx_ = m_idx;
    dist_mahalanobist2_ = dist_mahalanobist2;
  }
  
  void JCBB_TreeNode::getAssociation(size_t &z_idx, size_t &m_idx){
    z_idx = z_idx_;
    m_idx = m_idx_;
  }

  double JCBB_TreeNode::getMahalanobisDist2(){
    return dist_mahalanobist2_;
  }

  ::Eigen::VectorXd* JCBB_TreeNode::getInnovation(){
    return &innov_;
  }

  ::Eigen::MatrixXd* JCBB_TreeNode::getInnovationCovInv(){
    return &C_inverse_;
  }

}


// implementation JCBB

namespace rfs{

  template<class MeasurementModelType>
  JCBB<MeasurementModelType>::JCBB(){

  }

  template<class MeasurementModelType>
  JCBB<MeasurementModelType>::JCBB(double confidenceItvl,
				   MeasurementModelType *measurementModel,
				   std::vector< TMeasurement > *measurements,
				   TPose const *robot,
				   std::vector< TLandmark const * > *landmarks,
				   ::Eigen::MatrixXd const *estCovDense){
    
    confidence_interval_ = confidenceItvl;
    measurementModel_ = measurementModel;
    Z_actual_ptr_ = measurements;
    robot_ = robot;
    landmarks_ = landmarks;
    isEstCovBlockDiag_ = true;
    if(estCovDense != NULL){
      isEstCovBlockDiag_ = false;
      estCovDense_ = estCovDense;
    }
    nAssociations_best = 0;
    dist_mahalanobist2_best = DBL_MAX;
    node_best = NULL;

    calcExpectedMeasurements();

    branchList_.push(&iTreeRoot_);
    while(!branchList_.empty()){
      branch( branchList_.top() );
      branchList_.pop();
    }

  }

  template<class MeasurementModelType>
  JCBB<MeasurementModelType>::~JCBB(){

  }

  template<class MeasurementModelType>
  void JCBB<MeasurementModelType>::calcExpectedMeasurements(){

    size_t nLandmarks = landmarks_.size();
    TMeasurement z_predict; 
    Z_predict_.clear();
    Z_predict_.reserve( nLandmarks );
    jacobians_wrt_pose_.reserve( nLandmarks );
    jacobians_wrt_lmk_.reserve( nLandmarks );

    BOOST_FOREACH( TLandmark* lmkPtr, landmarks_){
      
      TJacobianPose *H_pose_ptr = new TJacobianPose;
      TJacobianLmk  *H_lmk_ptr = new TJacobianLmk;

      measurementModel_->measure(*robot_, *lmkPtr, z_predict, H_lmk_ptr, H_pose_ptr);
      Z_predict_.push_back(z_predict);
      jacobians_wrt_pose_.push_back( H_pose_ptr );
      jacobians_wrt_lmk_.push_back( H_lmk_ptr );

    }
  }

  template<class MeasurementModelType>
  void JCBB<MeasurementModelType>::branch(JCBB_TreeNode *node){

    // Check previous associations
    size_t z_idx = -1;
    size_t m_idx;
    JCBB_TreeNode *node_current = node;
    std::vector<int> isUsed_lmk_idx(0, landmarks_.size());
    std::vector<int> used_lmk_indices;
    while(node_current->getParent() != NULL){
      node_current->getAssociation(z_idx, m_idx);
      isUsed_lmk_idx[ m_idx ] = 1;
      used_lmk_indices.push_back( m_idx );
      node_current = node_current->getParent();
    }
    size_t nAssociations = used_lmk_indices.size();
    ::Eigen::VectorXd *innov_m = node->getInnovation();

    // Branch (for unused landmark indices), if not rejected by gating
    double d2_m = node->getMahalanobisDist2();
    node->getAssociation(z_idx, m_idx);
    z_idx++;
    m_idx = 0;
    size_t nInsideGate = 0;

    while(nInsideGate == 0 && z_idx < Z_actual_ptr_.size()){

      BOOST_FOREACH(int isUsed, isUsed_lmk_idx){
	if( isUsed == 0 ){
	
	  // innovation
	  typename TMeasurement::Vec z_actual = Z_actual_ptr_->at(z_idx).get();
	  typename TMeasurement::Vec z_predict = Z_predict_[m_idx].get();
	  typename TMeasurement::Vec innov = z_actual - z_predict; 
	  typename TMeasurement::Mat Q = measurementModel_->getNoise();

	  int robotDim = robot_->getNDim();
	  int lmkDim = landmarks_->at(0)->getNDim();
	  int measureDim = z_actual->size();
	  int nPreviousMeasurements = used_lmk_indices.size();
	  double dist_mahalanobist2 = -1;

	  ::Eigen::MatrixXd W;
	  W.resize(measureDim, measureDim * nPreviousMeasurements);

	  bool isRobotPoseHasUncertainty = true;

	  if(!isEstCovBlockDiag_){ // Estimate covariance is not block diagonal
	
	    ::Eigen::MatrixXd P_r_rows = estCovDense_.topRows(robotDim); // first row-block (robotDim rows) of estimate cov	
	    ::Eigen::MatrixXd P_lmk_rows = estCovDense_.block(robotDim + m_idx*lmkDim, 0, lmkDim, estCovDense_.cols() ); // the row-block of landmark m_idx 
	    ::Eigen::MatrixXd HP = *(jacobians_wrt_pose_[m_idx]) * P_r_rows + *(jacobians_wrt_lmk_[m_idx]) * P_lmk_rows;
	    ::Eigen::MatrixXd HP_cols_r = HP.leftCols(robotDim); // first col-block (robotDim cols) of HP

	    int W_block_idx = 0;
	    BOOST_REVERSE_FOREACH(int m_idx_used, used_lmk_indices){	 
	      ::Eigen::MatrixXd HP_cols_lmk = HP.block(0, robotDim + m_idx_used*lmkDim, HP.rows(), lmkDim);  // Pick off the col-block of HP for landmark m_idx_used
	      W.block(0, measureDim * W_block_idx, measureDim, measureDim) = HP_cols_r * jacobians_wrt_pose_[m_idx_used]->transpose() + 
		HP_cols_lmk * jacobians_wrt_lmk_[m_idx_used]->transpose() + Q;
	      W_block_idx++;
	    }
	  }else{ // Estimate covariance is block diagonal

	    ::Eigen::Matrix<double, TPose::Vec::RowsAtCompileTime, TPose::Vec::RowsAtCompileTime> P_rr = robot_->getCov();
	    if(!P_rr.isZero()){
	      for(int W_block_idx = 0; W_block_idx < nPreviousMeasurements; W_block_idx++){		
		W.block(0, measureDim * W_block_idx, measureDim, measureDim) = *(jacobians_wrt_pose_[m_idx]) * P_rr * jacobians_wrt_pose_[m_idx]->transpose() + Q;
	      }
	    }else{
	      W.setZero();
	      isRobotPoseHasUncertainty = false;
	    }

	  }

	  ::Eigen::Matrix<double, TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> C_current = *(jacobians_wrt_pose_[m_idx]) * estCovDense_.block(0, 0, robotDim, robotDim) * (jacobians_wrt_pose_[m_idx])->transpose() + *(jacobians_wrt_lmk_[m_idx]) * estCovDense_.block(robotDim + m_idx*lmkDim, robotDim + m_idx*lmkDim, lmkDim, lmkDim) * (jacobians_wrt_lmk_[m_idx])->transpose();

	  ::Eigen::Matrix<double, TMeasurement::Vec::RowsAtCompileTime, TMeasurement::Vec::RowsAtCompileTime> N;
	  ::Eigen::MatrixXd L;
	  ::Eigen::MatrixXd K;

	  if(isRobotPoseHasUncertainty){
	    
	    ::Eigen::MatrixXd *C_inverse_m = node->getInnovationCovInv(); 
	    N = ( C_current - W * *C_inverse_m * W.transpose() ).inverse();	
	    L = (-N * W) * *C_inverse_m;	  
	    K = *C_inverse_m + L.transpose() * N.inverse() * L;
	    dist_mahalanobist2 = d2_m + innov_m->transpose() * L.transpose() * N * L * *innov_m 
	      + 2 * innov.transpose() * L * *innov_m + innov * N * innov; 

	  }else{ // robot pose has no uncertainty

	    dist_mahalanobist2 = d2_m + innov.transpose() * C_current.inverse() * innov;
	    isRobotPoseHasUncertainty = false;

	  }

	  // Chi-square gating
	  int nDof = measureDim * (1 + nPreviousMeasurements);
	  boost::math::chi_squared_distribution<> chi2Dist(nDof);
	  double threshold = quantile(chi2Dist, confidence_interval_);
	  if(dist_mahalanobist2 < threshold * threshold){

	    // gating passed - new node in interpretation tree
	    JCBB_TreeNode* new_node = node->addChild();
	    new_node->setAssociation(z_idx, m_idx, dist_mahalanobist2);
	    *(new_node->getInnovation()) << *innov_m, innov;
	    if(isRobotPoseHasUncertainty){
	      *(new_node->getInnovationCovInv()) << K, L.transpose(), L, N;
	    }
	    
	    // add to priority queue
	    branchList_.push(new_node);
	    nInsideGate++;

	    // Keep track of best data association hypothesis
	    if( (nAssociations_best < nAssociations + 1) ||  
		( nAssociations_best == nAssociations + 1 && dist_mahalanobist2_best > dist_mahalanobist2) ){
	      nAssociations_best = nAssociations + 1;
	      dist_mahalanobist2_best = dist_mahalanobist2;
	      node_best = new_node;
	    }
	  }

	}

	m_idx++;

      }// BOOST_FOREACH
    
      z_idx++;

    }

  }

} // namespace rfs


#endif
