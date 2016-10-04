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

#define BOOST_NO_CXX11_SCOPED_ENUMS // required for boost/filesystem to work with C++11
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "ProcessModel_Odometry6D.hpp"
#include "RBPHDFilter.hpp"
#include "MeasurementModel_6D.hpp"
#include "Visualizer6D.hpp"
#include <stdio.h>
#include <string>
#include <sys/ioctl.h>

#ifdef _PERFTOOLS_CPU
#include <gperftools/profiler.h>
#endif
#ifdef _PERFTOOLS_HEAP
#include <gperftools/heap-profiler.h>
#endif



using namespace rfs;

/**
 * \class Simulator_RBPHDSLAM_6d
 * \brief A 6d SLAM Simulator using the RB-PHD Filter
 * \author Felipe Inostroza
 */
class Simulator_RBPHDSLAM_6d {

public:

	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	;

	Simulator_RBPHDSLAM_6d() {
		pFilter_ = NULL;
	}

	~Simulator_RBPHDSLAM_6d() {

		if (pFilter_ != NULL) {
			delete pFilter_;
		}

	}

	/** Read the simulator configuration file */
	bool readConfigFile(const char* fileName) {

		cfgFileName_ = fileName;

		boost::property_tree::ptree pt;
		boost::property_tree::xml_parser::read_xml(fileName, pt);

		logResultsToFile_ = false;
		if (pt.get("config.logging.logResultsToFile", 0) == 1)
			logResultsToFile_ = true;
		logTimingToFile_ = false;
		if (pt.get("config.logging.logTimingToFile", 0) == 1)
			logTimingToFile_ = true;
		logDirPrefix_ = pt.get<std::string>("config.logging.logDirPrefix",
				"./");
		if (*logDirPrefix_.rbegin() != '/')
			logDirPrefix_ += '/';

		use_gui_= false;
		if( pt.get<bool>("config.use_gui") != 0){
			use_gui_=true;
		}
		kMax_ = pt.get<int>("config.timesteps");
		dT_ = pt.get<double>("config.sec_per_timestep");
		dTimeStamp_ = TimeStamp(dT_);

		nSegments_ = pt.get<int>("config.trajectory.nSegments");
		max_dx_ = pt.get<double>("config.trajectory.max_dx_per_sec");
		max_dy_ = pt.get<double>("config.trajectory.max_dy_per_sec");
		max_dz_ = pt.get<double>("config.trajectory.max_dz_per_sec");
		max_dqx_ = pt.get<double>("config.trajectory.max_dqx_per_sec");
		max_dqy_ = pt.get<double>("config.trajectory.max_dqy_per_sec");
		max_dqz_ = pt.get<double>("config.trajectory.max_dqz_per_sec");
		max_dqw_ = pt.get<double>("config.trajectory.max_dqw_per_sec");
		min_dx_ = pt.get<double>("config.trajectory.min_dx_per_sec");
		vardx_ = pt.get<double>("config.trajectory.vardx");
		vardy_ = pt.get<double>("config.trajectory.vardy");
		vardz_ = pt.get<double>("config.trajectory.vardz");
		vardqx_ = pt.get<double>("config.trajectory.vardqx");
		vardqy_ = pt.get<double>("config.trajectory.vardqy");
		vardqz_ = pt.get<double>("config.trajectory.vardqz");
		vardqw_ = pt.get<double>("config.trajectory.vardqw");

		nLandmarks_ = pt.get<int>("config.landmarks.nLandmarks");
		varlmx_ = pt.get<double>("config.landmarks.varlmx");
		varlmy_ = pt.get<double>("config.landmarks.varlmy");
		varlmz_ = pt.get<double>("config.landmarks.varlmz");

		rangeLimitMax_ = pt.get<double>("config.measurements.rangeLimitMax");
		rangeLimitMin_ = pt.get<double>("config.measurements.rangeLimitMin");
		rangeLimitBuffer_ = pt.get<double>(
				"config.measurements.rangeLimitBuffer");
		Pd_ = pt.get<double>("config.measurements.probDetection");
		c_ = pt.get<double>("config.measurements.clutterIntensity");
		varzx_ = pt.get<double>("config.measurements.varzx");
		varzy_ = pt.get<double>("config.measurements.varzy");
		varzz_ = pt.get<double>("config.measurements.varzz");

		nParticles_ = pt.get("config.filter.nParticles", 200);

		pNoiseInflation_ = pt.get(
				"config.filter.predict.processNoiseInflationFactor", 1.0);
		birthGaussianWeight_ = pt.get(
				"config.filter.predict.birthGaussianWeight", 0.01);

		zNoiseInflation_ = pt.get(
				"config.filter.update.measurementNoiseInflationFactor", 1.0);
		innovationRangeThreshold_ = pt.get<double>(
				"config.filter.update.KalmanFilter.innovationThreshold.range");
		innovationBearingThreshold_ =
				pt.get<double>(
						"config.filter.update.KalmanFilter.innovationThreshold.bearing");
		newGaussianCreateInnovMDThreshold_ = pt.get<double>(
				"config.filter.update.GaussianCreateInnovMDThreshold");

		importanceWeightingEvalPointCount_ = pt.get(
				"config.filter.weighting.nEvalPt", 15);
		importanceWeightingEvalPointGuassianWeight_ = pt.get(
				"config.filter.weighting.minWeight", 0.75);
		importanceWeightingMeasurementLikelihoodMDThreshold_ = pt.get(
				"config.filter.weighting.threshold", 3.0);
		useClusterProcess_ = false;
		if (pt.get("config.filter.weighting.useClusterProcess", 0) == 1)
			useClusterProcess_ = true;

		effNParticleThreshold_ = pt.get("config.filter.resampling.effNParticle",
				nParticles_);
		minUpdatesBeforeResample_ = pt.get(
				"config.filter.resampling.minTimesteps", 1);

		gaussianMergingThreshold_ = pt.get<double>(
				"config.filter.merge.threshold");
		gaussianMergingCovarianceInflationFactor_ = pt.get(
				"config.filter.merge.covInflationFactor", 1.0);

		gaussianPruningThreshold_ = pt.get("config.filter.prune.threshold",
				birthGaussianWeight_);

		return true;
	}

	/** Generate a 6d trajectory in 3d space
	 *  \param[in] randSeed random seed for generating trajectory
	 */
	void generateTrajectory(int randSeed = 0) {

		srand48(randSeed);

		TimeStamp t;
		int seg = 0;
		MotionModel_Odometry6d::TState::Mat Q;
		Q.setZero();
		Q(0, 0) = vardx_;
		Q(1, 1) = vardy_;
		Q(2, 2) = vardz_;
		Q(3, 3) = vardqw_;
		Q(4, 4) = vardqx_;
		Q(5, 5) = vardqy_;
		Q(6, 6) = vardqz_;
		MotionModel_Odometry6d motionModel(Q);
		MotionModel_Odometry6d::TInput input_k(t);
		MotionModel_Odometry6d::TState pose_k(t);
		pose_k[6]=1.0;
		MotionModel_Odometry6d::TState pose_km(t);
		pose_km[6]=1.0;
		groundtruth_displacement_.reserve(kMax_);
		groundtruth_pose_.reserve(kMax_);
		groundtruth_displacement_.push_back(input_k);
		groundtruth_pose_.push_back(pose_k);

		for (int k = 1; k < kMax_; k++) {

			t += dTimeStamp_;

			if (k <= 50) {
				double dx = 0;
				double dy = 0;
				double dz = 0;
				double dqx = 0;
				double dqy = 0;
				double dqz = 0;
				double dqw = 1;
				MotionModel_Odometry6d::TInput::Vec d;
				MotionModel_Odometry6d::TInput::Vec dCovDiag;
				d << dx, dy, dz, dqx, dqy, dqz, dqw;
				dCovDiag << 0, 0, 0, 0, 0, 0, 0;
				input_k = MotionModel_Odometry6d::TInput(d,
						dCovDiag.asDiagonal(), k);
			} else if (k >= kMax_ / nSegments_ * seg) {
				seg++;
				double dx = drand48() * (max_dx_ - min_dx_) * dT_
						+ min_dx_ * dT_;

				double dy = (drand48() * max_dy_ * 2 - max_dy_) * dT_;
				double dz = (drand48() * max_dz_ * 2 - max_dz_) * dT_;
				double dqx = (drand48() * max_dqx_ * 2 - max_dqx_) * dT_;
				double dqy = (drand48() * max_dqy_ * 2 - max_dqy_) * dT_;
				double dqz = (drand48() * max_dqz_ * 2 - max_dqz_) * dT_;
				double dqw = 1+(drand48() * max_dqw_ * 2 - max_dqw_) * dT_;

				MotionModel_Odometry6d::TInput::Vec d;
				MotionModel_Odometry6d::TInput::Vec dCovDiag;
				d << dx, dy, dz, dqx, dqy, dqz, dqw;

				dCovDiag << Q(0, 0), Q(1, 1), Q(2, 2), Q(3, 3), Q(4, 4), Q(5,
						5), Q(6, 6);
				input_k = MotionModel_Odometry6d::TInput(d,
						dCovDiag.asDiagonal(), k);
			}

			groundtruth_displacement_.push_back(input_k);
			groundtruth_displacement_.back().setTime(t);

			MotionModel_Odometry6d::TState x_k;
			motionModel.step(x_k, groundtruth_pose_[k - 1], input_k,
					dTimeStamp_);
			groundtruth_pose_.push_back(x_k);
			groundtruth_pose_.back().setTime(t);

		}

	}

	/** Generate odometry measurements */
	void generateOdometry() {

		odometry_.reserve(kMax_);
		MotionModel_Odometry6d::TInput zero;
		MotionModel_Odometry6d::TInput::Vec u0;
		u0.setZero();
		zero.set(u0, 0);
		odometry_.push_back(zero);

		MotionModel_Odometry6d::TState::Mat Q;
		Q.setZero();
		Q(0, 0) = vardx_;
		Q(1, 1) = vardy_;
		Q(2, 2) = vardz_;
		Q(3, 3) = vardqw_;
		Q(4, 4) = vardqx_;
		Q(5, 5) = vardqy_;
		Q(6, 6) = vardqz_;
		MotionModel_Odometry6d motionModel(Q);
		deadReckoning_pose_.reserve(kMax_);
		deadReckoning_pose_.push_back(groundtruth_pose_[0]);

		TimeStamp t;

		for (int k = 1; k < kMax_; k++) {

			t += dTimeStamp_;
			double dt = dTimeStamp_.getTimeAsDouble();

			MotionModel_Odometry6d::TInput in = groundtruth_displacement_[k];
			MotionModel_Odometry6d::TState::Mat Qk = Q * dt * dt;
			in.setCov(Qk);
			MotionModel_Odometry6d::TInput out;
			in.sample(out);

			odometry_.push_back(out);

			MotionModel_Odometry6d::TState p;
			motionModel.step(p, deadReckoning_pose_[k - 1], odometry_[k],
					dTimeStamp_);
			p.setTime(t);
			deadReckoning_pose_.push_back(p);
		}

	}

	/** Generate landmarks */
	void generateLandmarks() {

		MeasurementModel_6D measurementModel(varzx_, varzy_, varzz_);
		MeasurementModel_6D::TPose pose;

		groundtruth_landmark_.reserve(nLandmarks_);

		int nLandmarksCreated = 0;
		for (int k = 1; k < kMax_; k++) {

			if (k >= kMax_ / nLandmarks_ * nLandmarksCreated) {

				MeasurementModel_6D::TPose pose;
				MeasurementModel_6D::TMeasurement measurementToCreateLandmark;
				MeasurementModel_6D::TMeasurement::Vec z;
				do {
					for (int i = 0; i < z.size(); i++)
						z(i) = drand48() * 2 * rangeLimitMax_ - rangeLimitMax_;
				} while (z.norm() > rangeLimitMax_);

				measurementToCreateLandmark.set(z);
				MeasurementModel_6D::TLandmark lm;

				measurementModel.inverseMeasure(groundtruth_pose_[k],
						measurementToCreateLandmark, lm);

				groundtruth_landmark_.push_back(lm);

				nLandmarksCreated++;

			}

		}

	}

	/** Generate landmark measurements */
	void generateMeasurements() {

		MeasurementModel_6D measurementModel(varzx_, varzy_, varzz_);
		MeasurementModel_6D::TMeasurement::Mat R;
		measurementModel.getNoise(R);
		measurementModel.config.rangeLimMax_ = rangeLimitMax_;
		measurementModel.config.rangeLimMin_ = rangeLimitMin_;
		measurementModel.config.probabilityOfDetection_ = Pd_;
		measurementModel.config.uniformClutterIntensity_ = c_;
		double meanClutter = measurementModel.clutterIntensityIntegral();

		double expNegMeanClutter = exp(-meanClutter);
		double poissonPmf[100];
		double poissonCmf[100];
		double mean_pow_i = 1;
		double i_factorial = 1;
		poissonPmf[0] = expNegMeanClutter;
		poissonCmf[0] = poissonPmf[0];
		for (int i = 1; i < 100; i++) {
			mean_pow_i *= meanClutter;
			i_factorial *= i;
			poissonPmf[i] = mean_pow_i / i_factorial * expNegMeanClutter;
			poissonCmf[i] = poissonCmf[i - 1] + poissonPmf[i];
		}

		lmkFirstObsTime_.resize(groundtruth_landmark_.size());
		for (int m = 0; m < lmkFirstObsTime_.size(); m++) {
			lmkFirstObsTime_[m] = -1;
		}

		TimeStamp t;

		for (int k = 1; k < kMax_; k++) {

			t += dTimeStamp_;

			groundtruth_pose_[k];

			// Real detections
			for (int m = 0; m < groundtruth_landmark_.size(); m++) {

				bool success;
				MeasurementModel_6D::TMeasurement z_m_k;
				success = measurementModel.sample(groundtruth_pose_[k],
						groundtruth_landmark_[m], z_m_k);
				if (success) {

					if ( drand48() <= Pd_) {
						z_m_k.setTime(t);
						// z_m_k.setCov(R);
						measurements_.push_back(z_m_k);
					}

					if (lmkFirstObsTime_[m] == -1) {
						lmkFirstObsTime_[m] = t.getTimeAsDouble();
					}
				}

			}

			// False alarms
			double randomNum = drand48();
			int nClutterToGen = 0;
			while (randomNum > poissonCmf[nClutterToGen]) {
				nClutterToGen++;
			}
			for (int i = 0; i < nClutterToGen; i++) {

				MeasurementModel_6D::TMeasurement z_clutter;
				MeasurementModel_6D::TMeasurement::Vec z;
				do {
					for (int i = 0; i < z.size(); i++)
						z(i) = drand48() * 2 * rangeLimitMax_ - rangeLimitMax_;
				} while (z.norm() < rangeLimitMax_ && z.norm() > rangeLimitMin_);

				z_clutter.set(z, t);
				measurements_.push_back(z_clutter);

			}

		}

	}

	/** Data Logging */
	void exportSimData() {

		if (logResultsToFile_ || logTimingToFile_) {
			boost::filesystem::path dir(logDirPrefix_);
			boost::filesystem::create_directories(dir);
			boost::filesystem::path cfgFilePathSrc(cfgFileName_);
			std::string cfgFileDst(logDirPrefix_);
			cfgFileDst += "simSettings.xml";
			boost::filesystem::path cfgFilePathDst(cfgFileDst.data());
			boost::filesystem::copy_file(cfgFilePathSrc, cfgFilePathDst,
					boost::filesystem::copy_option::overwrite_if_exists);
		}

		if (!logResultsToFile_)
			return;

		TimeStamp t;

		FILE* pGTPoseFile;
		std::string filenameGTPose(logDirPrefix_);
		filenameGTPose += "gtPose.dat";
		pGTPoseFile = fopen(filenameGTPose.data(), "w");
		MotionModel_Odometry6d::TState::Vec x;
		for (int i = 0; i < groundtruth_pose_.size(); i++) {
			groundtruth_pose_[i].get(x, t);
			fprintf(pGTPoseFile, "%f   %f   %f   %f   %f   %f   %f   %f\n",
					t.getTimeAsDouble(), x(0), x(1), x(2), x(3), x(4), x(5),
					x(6));
		}
		fclose(pGTPoseFile);

		FILE* pGTLandmarkFile;
		std::string filenameGTLandmark(logDirPrefix_);
		filenameGTLandmark += "gtLandmark.dat";
		pGTLandmarkFile = fopen(filenameGTLandmark.data(), "w");
		MeasurementModel_6D::TLandmark::Vec m;
		for (int i = 0; i < groundtruth_landmark_.size(); i++) {
			groundtruth_landmark_[i].get(m);
			fprintf(pGTLandmarkFile, "%f   %f   %f   %f\n", m(0), m(1), m(2),
					lmkFirstObsTime_[i]);
		}
		fclose(pGTLandmarkFile);

		FILE* pOdomFile;
		std::string filenameOdom(logDirPrefix_);
		filenameOdom += "odometry.dat";
		pOdomFile = fopen(filenameOdom.data(), "w");
		MotionModel_Odometry6d::TInput::Vec u;
		for (int i = 0; i < odometry_.size(); i++) {
			odometry_[i].get(u, t);
			fprintf(pOdomFile, "%f   %f   %f   %f   %f   %f   %f   %f\n",
					t.getTimeAsDouble(), u(0), u(1), u(2), u(3), u(4), u(5),
					u(6));
		}
		fclose(pOdomFile);

		FILE* pMeasurementFile;
		std::string filenameMeasurement(logDirPrefix_);
		filenameMeasurement += "measurement.dat";
		pMeasurementFile = fopen(filenameMeasurement.data(), "w");
		MeasurementModel_6D::TMeasurement::Vec z;
		for (int i = 0; i < measurements_.size(); i++) {
			measurements_[i].get(z, t);
			fprintf(pMeasurementFile, "%f   %f   %f   %f\n",
					t.getTimeAsDouble(), z(0), z(1), z(2));
		}
		fclose(pMeasurementFile);

		FILE* pDeadReckoningFile;
		std::string filenameDeadReckoning(logDirPrefix_);
		filenameDeadReckoning += "deadReckoning.dat";
		pDeadReckoningFile = fopen(filenameDeadReckoning.data(), "w");
		MotionModel_Odometry6d::TState::Vec odo;
		for (int i = 0; i < deadReckoning_pose_.size(); i++) {
			deadReckoning_pose_[i].get(odo, t);
			fprintf(pDeadReckoningFile, "%f   %f   %f   %f   %f   %f   %f   %f\n",
					t.getTimeAsDouble(), odo(0), odo(1), odo(2), odo(3), odo(4), odo(5), odo(6));
		}
		fclose(pDeadReckoningFile);

	}

	/** RB-PHD Filter Setup */
	void setupRBPHDFilter() {

		pFilter_ =
				new RBPHDFilter<MotionModel_Odometry6d,
						StaticProcessModel<Landmark3d>, MeasurementModel_6D,
						KalmanFilter<StaticProcessModel<Landmark3d>,
								MeasurementModel_6D> >(nParticles_);

		double dt = dTimeStamp_.getTimeAsDouble();

		// configure robot motion model (only need to set once since timesteps are constant)
		MotionModel_Odometry6d::TState::Mat Q;
		Q.setZero();
		Q(0, 0) = vardx_;
		Q(1, 1) = vardy_;
		Q(2, 2) = vardz_;
		Q(3, 3) = vardqx_;
		Q(4, 4) = vardqy_;
		Q(5, 5) = vardqz_;
		Q(6, 6) = vardqw_;
		Q *= (pNoiseInflation_ * dt * dt);
		pFilter_->getProcessModel()->setNoise(Q);

		// configure landmark process model (only need to set once since timesteps are constant)
		Landmark3d::Mat Q_lm;
		Q_lm.setZero();
		Q_lm(0, 0) = varlmx_;
		Q_lm(1, 1) = varlmy_;
		Q_lm(2, 2) = varlmz_;
		Q_lm = Q_lm * dt * dt;
		pFilter_->getLmkProcessModel()->setNoise(Q_lm);

		// configure measurement model
		MeasurementModel_6D::TMeasurement::Mat R;
		R << varzx_, 0, 0, 0, varzy_, 0, 0, 0, varzz_;
		R *= zNoiseInflation_;
		pFilter_->getMeasurementModel()->setNoise(R);
		pFilter_->getMeasurementModel()->config.probabilityOfDetection_ = Pd_;
		pFilter_->getMeasurementModel()->config.uniformClutterIntensity_ = c_;
		pFilter_->getMeasurementModel()->config.rangeLimMax_ = rangeLimitMax_;
		pFilter_->getMeasurementModel()->config.rangeLimMin_ = rangeLimitMin_;
		pFilter_->getMeasurementModel()->config.rangeLimBuffer_ =
				rangeLimitBuffer_;

		// configure the Kalman filter for landmark updates
		/* Thresholds not implemented!!
		 pFilter_->getKalmanFilter()->config.rangeInnovationThreshold_ = innovationRangeThreshold_;
		 pFilter_->getKalmanFilter()->config.bearingInnovationThreshold_ = innovationBearingThreshold_;
		 */
		// configure the filter
		pFilter_->config.birthGaussianWeight_ = birthGaussianWeight_;
		pFilter_->setEffectiveParticleCountThreshold(effNParticleThreshold_);
		pFilter_->config.minUpdatesBeforeResample_ = minUpdatesBeforeResample_;
		pFilter_->config.newGaussianCreateInnovMDThreshold_ =
				newGaussianCreateInnovMDThreshold_;
		pFilter_->config.importanceWeightingMeasurementLikelihoodMDThreshold_ =
				importanceWeightingMeasurementLikelihoodMDThreshold_;
		pFilter_->config.importanceWeightingEvalPointCount_ =
				importanceWeightingEvalPointCount_;
		pFilter_->config.importanceWeightingEvalPointGuassianWeight_ =
				importanceWeightingEvalPointGuassianWeight_;
		pFilter_->config.gaussianMergingThreshold_ = gaussianMergingThreshold_;
		pFilter_->config.gaussianMergingCovarianceInflationFactor_ =
				gaussianMergingCovarianceInflationFactor_;
		pFilter_->config.gaussianPruningThreshold_ = gaussianPruningThreshold_;
		pFilter_->config.useClusterProcess_ = useClusterProcess_;

		// Visualization
		if(use_gui_){
			visualizer = new Visualizer6D();
			visualizer->setup(groundtruth_landmark_,groundtruth_pose_,deadReckoning_pose_);
			visualizer->start();
		}
	}

	/** Run the simulator */
	void run() {

		printf("Running simulation\n\n");

#ifdef _PERFTOOLS_CPU
		std::string perfCPU_file = logDirPrefix_ + "rbphdslam2dSim_cpu.prof";
		ProfilerStart(perfCPU_file.data());
#endif
#ifdef _PERFTOOLS_HEAP
		std::string perfHEAP_file = logDirPrefix_ + "rbphdslam2dSim_heap.prof";
		HeapProfilerStart(perfHEAP_file.data());
#endif

		//////// Initialization at first timestep //////////

		if (!logResultsToFile_) {
			std::cout
					<< "Note: results are NOT being logged to file (see config xml file)\n";
		}
		FILE* pParticlePoseFile;
		if (logResultsToFile_) {
			std::string filenameParticlePoseFile(logDirPrefix_);
			filenameParticlePoseFile += "particlePose.dat";
			pParticlePoseFile = fopen(filenameParticlePoseFile.data(), "w");
		}
		FILE* pLandmarkEstFile;
		if (logResultsToFile_) {
			std::string filenameLandmarkEstFile(logDirPrefix_);
			filenameLandmarkEstFile += "landmarkEst.dat";
			pLandmarkEstFile = fopen(filenameLandmarkEstFile.data(), "w");
		}
		MotionModel_Odometry6d::TState x_i;
		int zIdx = 0;

		if (logResultsToFile_) {
			for (int i = 0; i < pFilter_->getParticleCount(); i++) {
				x_i = *(pFilter_->getParticleSet()->at(i));
				fprintf(pParticlePoseFile,
						"%f   %d   %f   %f   %f   %f   %f   %f   %f   1.0\n",
						0.0, i, x_i.get(0), x_i.get(1), x_i.get(2), x_i.get(3),
						x_i.get(4), x_i.get(5), x_i.get(6));
			}
		}

		/////////// Run simulator from k = 1 to kMax_ /////////

		TimeStamp time;

		for (int k = 1; k < kMax_; k++) {

			time += dTimeStamp_;

			if (k % 100 == 0 || k == kMax_ - 1) {
				float progressPercent = float(k + 1) / float(kMax_);
				int progressBarW = 50;
				struct winsize ws;
				if (ioctl(1, TIOCGWINSZ, &ws) >= 0)
					progressBarW = ws.ws_col - 30;
				int progressPos = progressPercent * progressBarW;
				if (progressBarW >= 50) {
					std::cout << "[";
					for (int i = 0; i < progressBarW; i++) {
						if (i < progressPos)
							std::cout << "=";
						else if (i == progressPos)
							std::cout << ">";
						else
							std::cout << " ";
					}
					std::cout << "] ";
				}
				std::cout << "k = " << k << " (" << int(progressPercent * 100.0)
						<< " %)\r";
				std::cout.flush();
			}
			if (k == kMax_ - 1)
				std::cout << std::endl << std::endl;

#ifdef _PERFTOOLS_HEAP
			if( k % 20 == 0)
			HeapProfilerDump("Timestep interval dump");
#endif

			////////// Prediction Step //////////

			// configure robot motion model ( not necessary since in simulation, timesteps are constant)
			// MotionModel_Odometry2d::TState::Mat Q;
			// Q << vardx_, 0, 0, 0, vardy_, 0, 0, 0, vardz_;
			// Q *= (pNoiseInflation_ * dt * dt);
			// pFilter_->getProcessModel()->setNoise(Q);

			// configure landmark process model ( not necessary since in simulation, timesteps are constant)
			// Landmark2d::Mat Q_lm;
			// Q_lm << varlmx_, 0, 0, varlmy_;
			// Q_lm = Q_lm * dt * dt;
			// pFilter_->getLmkProcessModel()->setNoise(Q_lm);

			pFilter_->predict(odometry_[k], dTimeStamp_);

			if (k <= 100) {
				for (int i = 0; i < nParticles_; i++)
					pFilter_->setParticlePose(i, groundtruth_pose_[k]);
			}

			// Prepare measurement vector for update
			std::vector<MeasurementModel_6D::TMeasurement> Z;
			TimeStamp kz = measurements_[zIdx].getTime();
			while (kz == time) {
				Z.push_back(measurements_[zIdx]);
				zIdx++;
				if (zIdx >= measurements_.size())
					break;
				kz = measurements_[zIdx].getTime();
			}

			////////// Update Step //////////
			pFilter_->update(Z);

			// Log particle poses
			int i_w_max = 0;
			double w_max = 0;
			if (logResultsToFile_) {
				for (int i = 0; i < pFilter_->getParticleCount(); i++) {
					x_i = *(pFilter_->getParticleSet()->at(i));
					double w = pFilter_->getParticleSet()->at(i)->getWeight();
					if (w > w_max) {
						i_w_max = i;
						w_max = w;
					}
					fprintf(pParticlePoseFile, "%f   %d   %f   %f   %f   %f   %f   %f   %f   %f\n",
							time.getTimeAsDouble(), i, x_i.get(0), x_i.get(1),
							x_i.get(2), x_i.get(3), x_i.get(4), x_i.get(5), x_i.get(6), w);
				}
				fprintf(pParticlePoseFile, "\n");
			}

			// Log landmark estimates
			if (logResultsToFile_) {

				int mapSize = pFilter_->getGMSize(i_w_max);
				for (int m = 0; m < mapSize; m++) {
					MeasurementModel_6D::TLandmark::Vec u;
					MeasurementModel_6D::TLandmark::Mat S;
					double w;
					pFilter_->getLandmark(i_w_max, m, u, S, w);

					fprintf(pLandmarkEstFile, "%f   %d   ",
							time.getTimeAsDouble(), i_w_max);
					fprintf(pLandmarkEstFile, "%f   %f   %f      ", u(0), u(1), u(2));
					fprintf(pLandmarkEstFile, "%f   %f   %f   %f   %f   %f", S(0, 0), S(0, 1), S(0, 2),
							S(1, 1), S(1, 2), S(2, 2));
					fprintf(pLandmarkEstFile, "   %f\n", w);

				}
			}

			// Visualization
			if(use_gui_){
				// diplay particle poses

				visualizer->update(pFilter_);


			}
		}

#ifdef _PERFTOOLS_HEAP
		HeapProfilerStop();
#endif
#ifdef _PERFTOOLS_CPU
		ProfilerStop();
#endif

		std::cout << "Elapsed Timing Information [nsec]\n";
		std::cout << std::setw(15) << std::left << "Prediction" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->predict_wall << std::setw(6)
				<< std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->predict_cpu << std::endl;
		std::cout << std::setw(15) << std::left << "Map Update" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapUpdate_wall << std::setw(6)
				<< std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapUpdate_cpu << std::endl;
		std::cout << std::setw(15) << std::left << "Map Update (KF)"
				<< std::setw(15) << std::setw(6) << std::right << "wall:"
				<< std::setw(15) << pFilter_->getTimingInfo()->mapUpdate_kf_wall
				<< std::setw(6) << std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapUpdate_kf_cpu << std::endl;
		std::cout << std::setw(15) << std::left << "Weighting" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->particleWeighting_wall
				<< std::setw(6) << std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->particleWeighting_cpu
				<< std::endl;
		std::cout << std::setw(15) << std::left << "Map Merge" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapMerge_wall << std::setw(6)
				<< std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapMerge_cpu << std::endl;
		std::cout << std::setw(15) << std::left << "Map Prune" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapPrune_wall << std::setw(6)
				<< std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->mapPrune_cpu << std::endl;
		std::cout << std::setw(15) << std::left << "Resampling" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->particleResample_wall
				<< std::setw(6) << std::left << std::right << "cpu:"
				<< std::setw(15)
				<< pFilter_->getTimingInfo()->particleResample_cpu << std::endl;
		std::cout << std::setw(15) << std::left << "Total" << std::setw(15)
				<< std::setw(6) << std::right << "wall:" << std::setw(15)
				<< pFilter_->getTimingInfo()->predict_wall
						+ pFilter_->getTimingInfo()->mapUpdate_wall
						+ pFilter_->getTimingInfo()->particleWeighting_wall
						+ pFilter_->getTimingInfo()->mapMerge_wall
						+ pFilter_->getTimingInfo()->mapPrune_wall
						+ pFilter_->getTimingInfo()->particleResample_wall
				<< std::setw(6) << std::right << "cpu:" << std::setw(15)
				<< pFilter_->getTimingInfo()->predict_cpu
						+ pFilter_->getTimingInfo()->mapUpdate_cpu
						+ pFilter_->getTimingInfo()->particleWeighting_cpu
						+ pFilter_->getTimingInfo()->mapMerge_cpu
						+ pFilter_->getTimingInfo()->mapPrune_cpu
						+ pFilter_->getTimingInfo()->particleResample_cpu
				<< std::endl;

		if (logTimingToFile_) {
			std::ofstream timingFile((logDirPrefix_ + "timing.dat").data());
			timingFile << "Elapsed Timing Information [nsec]\n";
			timingFile << std::setw(15) << std::left << "Prediction"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15) << pFilter_->getTimingInfo()->predict_wall
					<< std::setw(6) << std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->predict_cpu << std::endl;
			timingFile << std::setw(15) << std::left << "Map Update"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15)
					<< pFilter_->getTimingInfo()->mapUpdate_wall << std::setw(6)
					<< std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->mapUpdate_cpu << std::endl;
			timingFile << std::setw(15) << std::left << "Map Update (KF)"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15)
					<< pFilter_->getTimingInfo()->mapUpdate_kf_wall
					<< std::setw(6) << std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->mapUpdate_kf_cpu << std::endl;
			timingFile << std::setw(15) << std::left << "Weighting"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15)
					<< pFilter_->getTimingInfo()->particleWeighting_wall
					<< std::setw(6) << std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->particleWeighting_cpu
					<< std::endl;
			timingFile << std::setw(15) << std::left << "Map Merge"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15) << pFilter_->getTimingInfo()->mapMerge_wall
					<< std::setw(6) << std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->mapMerge_cpu << std::endl;
			timingFile << std::setw(15) << std::left << "Map Prune"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15) << pFilter_->getTimingInfo()->mapPrune_wall
					<< std::setw(6) << std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->mapPrune_cpu << std::endl;
			timingFile << std::setw(15) << std::left << "Resampling"
					<< std::setw(15) << std::setw(6) << std::right << "wall:"
					<< std::setw(15)
					<< pFilter_->getTimingInfo()->particleResample_wall
					<< std::setw(6) << std::left << std::right << "cpu:"
					<< std::setw(15)
					<< pFilter_->getTimingInfo()->particleResample_cpu
					<< std::endl;
			timingFile << std::setw(15) << std::left << "Total" << std::setw(15)
					<< std::setw(6) << std::right << "wall:" << std::setw(15)
					<< pFilter_->getTimingInfo()->predict_wall
							+ pFilter_->getTimingInfo()->mapUpdate_wall
							+ pFilter_->getTimingInfo()->particleWeighting_wall
							+ pFilter_->getTimingInfo()->mapMerge_wall
							+ pFilter_->getTimingInfo()->mapPrune_wall
							+ pFilter_->getTimingInfo()->particleResample_wall
					<< std::setw(6) << std::right << "cpu:" << std::setw(15)
					<< pFilter_->getTimingInfo()->predict_cpu
							+ pFilter_->getTimingInfo()->mapUpdate_cpu
							+ pFilter_->getTimingInfo()->particleWeighting_cpu
							+ pFilter_->getTimingInfo()->mapMerge_cpu
							+ pFilter_->getTimingInfo()->mapPrune_cpu
							+ pFilter_->getTimingInfo()->particleResample_cpu
					<< std::endl;
			timingFile.close();
		}

		if (logResultsToFile_) {
			fclose(pParticlePoseFile);
			fclose(pLandmarkEstFile);
		}
	}

private:

	const char* cfgFileName_;

	int kMax_; /**< number of timesteps */
	double dT_; /**< duration of timestep in seconds */
	TimeStamp dTimeStamp_; /**< duration of timestep in timestamp */

	// Trajectory
	int nSegments_;
	double max_dx_;
	double max_dy_;
	double max_dz_;
	double max_dqx_;
	double max_dqy_;
	double max_dqz_;
	double max_dqw_;
	double min_dx_;
	double vardx_;
	double vardy_;
	double vardz_;
	double vardqx_;
	double vardqy_;
	double vardqz_;
	double vardqw_;
	std::vector<MotionModel_Odometry6d::TInput> groundtruth_displacement_;
	std::vector<MotionModel_Odometry6d::TState> groundtruth_pose_;
	std::vector<MotionModel_Odometry6d::TInput> odometry_;
	std::vector<MotionModel_Odometry6d::TState> deadReckoning_pose_;

	// Landmarks
	int nLandmarks_;
	std::vector<MeasurementModel_6D::TLandmark> groundtruth_landmark_;
	double varlmx_;
	double varlmy_;
	double varlmz_;
	std::vector<double> lmkFirstObsTime_;

	// Range-Bearing Measurements
	double rangeLimitMax_;
	double rangeLimitMin_;
	double rangeLimitBuffer_;
	double Pd_;
	double c_;
	double varzx_;
	double varzy_;
	double varzz_;
	std::vector<MeasurementModel_6D::TMeasurement> measurements_;

	// Filters
	KalmanFilter<StaticProcessModel<Landmark3d>, MeasurementModel_6D> kf_;
	RBPHDFilter<MotionModel_Odometry6d, StaticProcessModel<Landmark3d>,
			MeasurementModel_6D,
			KalmanFilter<StaticProcessModel<Landmark3d>, MeasurementModel_6D> > *pFilter_;
	int nParticles_;
	double pNoiseInflation_;
	double zNoiseInflation_;
	double innovationRangeThreshold_;
	double innovationBearingThreshold_;
	double birthGaussianWeight_;
	double newGaussianCreateInnovMDThreshold_;
	double importanceWeightingMeasurementLikelihoodMDThreshold_;
	double importanceWeightingEvalPointGuassianWeight_;
	double effNParticleThreshold_;
	int minUpdatesBeforeResample_;
	double gaussianMergingThreshold_;
	double gaussianMergingCovarianceInflationFactor_;
	double gaussianPruningThreshold_;
	int importanceWeightingEvalPointCount_;
	bool useClusterProcess_;

	bool logResultsToFile_;
	bool logTimingToFile_;

	// 3D visualization
	bool use_gui_;
	Visualizer6D *visualizer;

public:
	std::string logDirPrefix_;
};

int main(int argc, char* argv[]) {

	Simulator_RBPHDSLAM_6d sim;

	int seed = time(NULL);
	srand(seed);
	int trajNum = rand();
	std::string cfgFileName;
	boost::program_options::options_description desc("Options");
	desc.add_options()("help,h", "produce this help message")("cfg,c",
			boost::program_options::value<std::string>(&cfgFileName)->default_value(
					"cfg/rbphdslam6dSim.xml"), "configuration xml file")(
			"trajectory,t", boost::program_options::value<int>(&trajNum),
			"trajectory number (default: a random integer)")("seed,s",
			boost::program_options::value<int>(&seed),
			"random seed for running the simulation (default: based on current system time)");
	boost::program_options::variables_map vm;
	boost::program_options::store(
			boost::program_options::parse_command_line(argc, argv, desc), vm);
	boost::program_options::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}

	if (vm.count("cfg")) {
		cfgFileName = vm["cfg"].as<std::string>();
	}
	std::cout << "Configuration file: " << cfgFileName << std::endl;
	if (!sim.readConfigFile(cfgFileName.data())) {
		return -1;
	}

	if (vm.count("trajectory")) {
		trajNum = vm["trajectory"].as<int>();
	}
	std::cout << "Trajectory: " << trajNum << std::endl;
	sim.generateTrajectory(trajNum);

	sim.generateLandmarks();
	sim.generateOdometry();
	sim.generateMeasurements();
	sim.exportSimData();
	sim.setupRBPHDFilter();

	if (vm.count("seed")) {
		seed = vm["seed"].as<int>();
		std::cout << "Simulation random seed manually set to: " << seed
				<< std::endl;
	}
	srand48(seed);

	// boost::timer::auto_cpu_timer *timer = new boost::timer::auto_cpu_timer(6, "Simulation run time: %ws\n");

	sim.run();

	// std::cout << "mem use: " << MemProfile::getCurrentRSS() << "(" << MemProfile::getPeakRSS() << ")\n";
	//delete timer;

	return 0;

}
