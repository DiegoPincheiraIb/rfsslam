/*
 * Software License Agreement (New BSD License)
 *
 * Copyright (c) 2014, Keith Leung, Felipe Inostroza
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
#ifndef VECTORGLMBSLAM_HPP
#define VECTORGLMBSLAM_HPP

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Timer.hpp"
#include <Eigen/Core>
#include "MurtyAlgorithm.hpp"
#include "PermutationLexicographic.hpp"
#include "RFSPSOParticle.hpp"
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include "OSPA.hpp"
#include "RandomVecMathTools.hpp"
#include "ceres/ceres.h"
#include <unordered_map>

#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include "g2o/core/robust_kernel_impl.h"

#include "g2o/types/slam2d/vertex_point_xy.h"
#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/types/slam2d/edge_pointxy.h"
#include "g2o/types/slam2d/edge_se2_pointxy.h"
#include "g2o/types/slam2d/edge_se2.h"

#ifdef _PERFTOOLS_CPU
#include <gperftools/profiler.h>
#endif
#ifdef _PERFTOOLS_HEAP
#include <gperftools/heap-profiler.h>
#endif

namespace rfs {

/**
 *  The weight and index of a landmarks from which a measurement can come from
 */
struct AssociationProbability{
    int i; /**< index of a landmark */
    double l; /**< log probability of association*/
};


/**
 * Struct to store a single component of a VGLMB , with its own g2o optimizer
 */
struct VectorGLMBComponent2D{

    typedef g2o::VertexPointXY  PointType;
    typedef g2o::VertexSE2  PoseType;
    typedef g2o::EdgeSE2PointXY MeasurementEdge;
    typedef g2o::EdgePointXY PointAnchorEdge;





     typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
     typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
     g2o::SparseOptimizer * optimizer_;
     g2o::OptimizationAlgorithmLevenberg *solverLevenberg_;
     SlamLinearSolver * linearSolver_;



     std::vector<std::vector<int> > DA_; /**< vector containing data association hypothesis,
     -1 meaning unknown data association (should not be used but kept for consistency with the known data associations), -2 known to be false alarm */
     std::vector<std::vector<MeasurementEdge*> > Z_; /**< Measurement edges stored, in order to set data association and add to graph later */
     std::vector<std::vector<AssociationProbability>> DAProbs_; /**< The association probability of each measurement, used for switching using gibbs sampling*/
     std::vector<PoseType*> poses_;
     double logweight_;
     int numPoses_,numPoints_;
};

/**
 *  \class VectorGLMBSLAM2D
 *  \brief Random Finite Set  optimization using ceres solver for  feature based SLAM
 *
 *
 *  \author  Felipe Inostroza
 */
    class VectorGLMBSLAM2D  {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW

      typedef g2o::VertexPointXY  PointType;
      typedef g2o::VertexSE2  PoseType;
      typedef g2o::EdgeSE2 OdometryEdge;
      typedef g2o::EdgeSE2PointXY MeasurementEdge;
      typedef g2o::EdgePointXY PointAnchorEdge;
      typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  SlamBlockSolver;
      typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> SlamLinearSolver;
      /**
       * \brief Configurations for this RFSBatchPSO optimizer
       */
      struct Config {


        /** The threshold used to determine if a possible meaurement-landmark
         *  pairing is significant to worth considering
         */
        double MeasurementLikelihoodThreshold_;

        double mapFromMeasurementProb_; /**< probability that each measurement will initialize a landmark on map initialization*/

        int numLandmarks_;

        int lmExistenceProb_;

        int numComponents_;



      } config;

      /**
       * Constructor
       */
      VectorGLMBSLAM2D ();

      /** Destructor */
      ~VectorGLMBSLAM2D ();



      /**
     *  Load a g2o style file , store groundtruth data association.
     * @param s Input stream
     */
    void
    load(std::ifstream s);

    /**
     * initialize the components , set the initial data associations to all false alarms
     */
    void initComponents();

    /**
     * Use the data association stored in DA_ to create the graph.
     * @param c the GLMB component
     */
    void constructGraph(VectorGLMBComponent2D &c);
    /**
     *
     * Initialize a VGLMB component , setting the data association to all false alarms
     * @param c the GLMB component
     */
    void init(VectorGLMBComponent2D &c);


    /**
     * Calculate the probability of each measurement being associated with a specific landmark
     * @param c the GLMB component
     */
    void updateDAProbs(VectorGLMBComponent2D &c);



      int nThreads_; /**< Number of threads  */

      VectorGLMBComponent2D gt_graph;

      std::vector<VectorGLMBComponent2D> components_; /**< VGLMB components */


    };


  //////////////////////////////// Implementation ////////////////////////



  VectorGLMBSLAM2D::VectorGLMBSLAM2D () {
      nThreads_ = 1;

#ifdef _OPENMP
      nThreads_ = omp_get_max_threads();
#endif

    }

  VectorGLMBSLAM2D::~VectorGLMBSLAM2D () {


    }




  void
  VectorGLMBSLAM2D::load(std::ifstream s){

      gt_graph.optimizer_->load(s);
  }





inline void VectorGLMBSLAM2D::initComponents() {
    components_.resize(config.numComponents_);

    for(auto &c:components_){
        init(c);
        constructGraph(c);
    }

}

inline void VectorGLMBSLAM2D::updateDAProbs(VectorGLMBComponent2D& c){

    for(PoseType p:c.poses_){

    }

}

inline void VectorGLMBSLAM2D::constructGraph(VectorGLMBComponent2D& c) {


    c.numPoses_  = 0;
    c.numPoints_ = 0;
    //Copy Vertices from optimizer with data association
    for(g2o::OptimizableGraph::Vertex*  v:gt_graph.optimizer_->vertices()){
        PoseType* pose = dynamic_cast<PoseType>(v);
        if (pose != NULL) {
            PoseType*  poseCopy= new PoseType();
            double poseData[3];
            pose->getEstimateData(poseData);
            poseCopy->setEstimateData(poseData);
            poseCopy->setId(pose->id());
            c.optimizer_->addVertex(poseCopy);
            c.poses_.push_back(poseCopy);
            c.numPoses_++;
        }
        PointType* point = dynamic_cast<PoseType>(v);
        if (point != NULL) {
            PointType*  pointCopy= new PointType();
            double pointData[2];
            point->getEstimateData(pointData);
            pointCopy->setEstimateData(pointData);
            pointCopy->setId(point->id());
            c.optimizer_->addVertex(pointCopy);
            c.numPoints_++;
        }

    }
    //Copy odometry measurements, Copy and save landmark measurements


    c.DA_.resize(c.numPoses_);
    c.Z_.resize(c.numPoses_);
    c.DAProbs_(c.numPoses_);

    for(g2o::OptimizableGraph::Edge * e:gt_graph.optimizer_->edges()){
        OdometryEdge* odo = dynamic_cast<OdometryEdge>(e);
        if(odo!=NULL){
            OdometryEdge* odocopy = new OdometryEdge();
            int firstvertex = odo->vertex(0)->id();
            odocopy->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(c.optimizer_->vertices().find(firstvertex)->second));
            int secondvertex = odo->vertex(1)->id();
            odocopy->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(c.optimizer_->vertices().find(secondvertex)->second));
            double measurementData[3];
            odo->getMeasurementData(measurementData);
            odocopy->setMeasurementData(measurementData);
            odocopy->setInformation(odo->information());
            odocopy->setParameterId(0, 0);
            c.optimizer_->addEdge(odocopy);
        }

        MeasurementEdge* z = dynamic_cast<MeasurementEdge>(e);
        if(z!=NULL){
            MeasurementEdge* zcopy = new MeasurementEdge();
            int firstvertex = z->vertex(0)->id();
            zcopy->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(c.optimizer_->vertices().find(firstvertex)->second));
            int secondvertex = z->vertex(1)->id();
            // zcopy->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(c.optimizer_->vertices().find(secondvertex)->second));
            double measurementData[2];
            z->getMeasurementData(measurementData);
            zcopy->setMeasurementData(measurementData);
            zcopy->setInformation(odo->information());
            zcopy->setParameterId(0, 0);

            c.DA_[firstvertex].push_back(-2);
            c.Z_[firstvertex].push_back(zcopy);
        }

    }




}



inline void VectorGLMBSLAM2D::init(VectorGLMBComponent2D& c) {
    auto linearSolver = g2o::make_unique<SlamLinearSolver>();
    linearSolver->setBlockOrdering(false);
    c.linearSolver_   = linearSolver.get();
    c.solverLevenberg_ =  new g2o::OptimizationAlgorithmLevenberg(
            g2o::make_unique<SlamBlockSolver>(std::move(linearSolver))
          );
    c.optimizer_->setAlgorithm(c.solverLevenberg_);
}

}
#endif
