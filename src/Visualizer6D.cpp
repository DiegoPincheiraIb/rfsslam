

#include "Visualizer6D.hpp"


namespace rfs
{

Visualizer6D::Visualizer6D(){

	display_mutex_ = new std::mutex();
}

void Visualizer6D::start(){

	display_thread_ = new std::thread(&Visualizer6D::run,this);
}
void Visualizer6D::setup(const std::vector<MeasurementModel_6D::TLandmark> &groundtruth_landmark,
		const std::vector<MotionModel_Odometry6d::TState> &groundtruth_pose,
		const std::vector<MotionModel_Odometry6d::TState> &deadreckoning_pose){

	groundtruth_pose_ = &groundtruth_pose;

	vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
	// chart->GetAxis(0)->SetGridVisible(true);
	// chart->GetAxis(0)->SetGridVisible(true);

	vtkSmartPointer<vtkAxesActor> axesActor =
            vtkSmartPointer<vtkAxesActor>::New();
	renderer_ = vtkSmartPointer<vtkRenderer>::New();
	renderer_ -> AddActor(axesActor);
	renderWindow_ = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow_->AddRenderer(renderer_);
	renderWindowInteractor_ = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor_->SetRenderWindow(renderWindow_);

	sphereSource_ = vtkSmartPointer<vtkSphereSource>::New();

	mapPoints_ = vtkSmartPointer<vtkPoints>::New();
	mapColors_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
	mapPolydata_ = vtkSmartPointer<vtkPolyData>::New();
	mapPolydata_->SetPoints(mapPoints_);
	mapPolydata_->GetPointData()->SetScalars(mapColors_);
	mapGlyph3D_ = vtkSmartPointer<vtkGlyph3D>::New();
	mapGlyph3D_->SetSourceConnection(sphereSource_->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
	mapGlyph3D_->SetInput(mapPolydata_);
#else
        mapGlyph3D_->SetInputData(mapPolydata_);
#endif

	mapGlyph3D_->SetColorModeToColorByScalar();
	mapGlyph3D_->ScalingOn();
	mapGlyph3D_->SetScaleModeToDataScalingOff();
	mapGlyph3D_->SetScaleFactor(0.7);
	mapGlyph3D_->Update();
	mapMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapMapper_->SetInputConnection(mapGlyph3D_->GetOutputPort());
	mapActor_ = vtkSmartPointer<vtkActor>::New();
	mapActor_->SetMapper(mapMapper_);



	gtmapPoints_ = vtkSmartPointer<vtkPoints>::New();
	gtmapColors_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
	gtmapPolydata_ = vtkSmartPointer<vtkPolyData>::New();
	gtmapPolydata_->SetPoints(gtmapPoints_);
	gtmapPolydata_->GetPointData()->SetScalars(gtmapColors_);
	gtmapGlyph3D_ = vtkSmartPointer<vtkGlyph3D>::New();
	gtmapGlyph3D_->SetSourceConnection(sphereSource_->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
	gtmapGlyph3D_->SetInput(gtmapPolydata_);
#else
        gtmapGlyph3D_->SetInputData(gtmapPolydata_);
#endif

	gtmapGlyph3D_->SetColorModeToColorByScalar();
	gtmapGlyph3D_->ScalingOff();
	gtmapGlyph3D_->Update();
	gtmapMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
	gtmapMapper_->SetInputConnection(gtmapGlyph3D_->GetOutputPort());
	gtmapActor_ = vtkSmartPointer<vtkActor>::New();
	gtmapActor_->SetMapper(gtmapMapper_);
	gtmapActor_->GetProperty()->SetOpacity(0.5);



	particlePoints_ = vtkSmartPointer<vtkPoints>::New();
	particleColors_ = vtkSmartPointer<vtkUnsignedCharArray>::New();
	particlePolydata_ = vtkSmartPointer<vtkPolyData>::New();
	particlePolydata_->SetPoints(particlePoints_);
	particlePolydata_->GetPointData()->SetScalars(particleColors_);
	particleGlyph3D_ = vtkSmartPointer<vtkGlyph3D>::New();
	particleGlyph3D_->SetSourceConnection(sphereSource_->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
        particleGlyph3D_->SetInput(particlePolydata_);
#else
        particleGlyph3D_->SetInputData (particlePolydata_);
#endif

	particleGlyph3D_->SetColorModeToColorByScalar();
	particleGlyph3D_->ScalingOn();
	particleGlyph3D_->SetScaleModeToDataScalingOff();
	particleGlyph3D_->SetScaleFactor(0.2);
	particleGlyph3D_->Update();
	particleMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
	particleMapper_->SetInputConnection(particleGlyph3D_->GetOutputPort());
	particleActor_ = vtkSmartPointer<vtkActor>::New();
	particleActor_->SetMapper(particleMapper_);



	gtTrajectoryPoints_ = vtkSmartPointer<vtkPoints>::New();
	gtTrajectoryPoints_->SetNumberOfPoints(groundtruth_pose_->size());

	for (int i=0 ;i< groundtruth_pose.size();i++){
		gtTrajectoryPoints_->SetPoint(i, groundtruth_pose[i].get(0), groundtruth_pose[i].get(1), groundtruth_pose[i].get(2));
	}
	//gtTrajectoryPolyLine_ = vtkSmartPointer<vtkPolyLine>::New();
	gtTrajectoryCells_ = vtkSmartPointer<vtkCellArray>::New();
	gtTrajectoryCells_->InsertNextCell(1);
	gtTrajectoryCells_->InsertCellPoint(0);

	gtTrajectoryPolydata_= vtkSmartPointer<vtkPolyData>::New();
	gtTrajectoryPolydata_->SetPoints(gtTrajectoryPoints_);
	gtTrajectoryPolydata_->SetLines(gtTrajectoryCells_);
	gtTrajectoryMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	gtTrajectoryMapper_->SetInput(gtTrajectoryPolydata_);
#else
        gtTrajectoryMapper_->SetInputData(gtTrajectoryPolydata_);
#endif

	gtTrajectoryMapper_->Update();
	gtTrajectoryActor_ =  vtkSmartPointer<vtkActor>::New();
	gtTrajectoryActor_->SetMapper(gtTrajectoryMapper_);
	gtTrajectoryActor_->GetProperty()->SetColor(0,191/255.0,255/255.0);

	drTrajectoryPoints_ = vtkSmartPointer<vtkPoints>::New();
	drTrajectoryPoints_->SetNumberOfPoints(deadreckoning_pose.size());

	for (int i=0 ;i< deadreckoning_pose.size();i++){
		drTrajectoryPoints_->SetPoint(i, deadreckoning_pose[i].get(0), deadreckoning_pose[i].get(1), deadreckoning_pose[i].get(2));
	}
	//gtTrajectoryPolyLine_ = vtkSmartPointer<vtkPolyLine>::New();
	drTrajectoryCells_ = vtkSmartPointer<vtkCellArray>::New();
	drTrajectoryCells_->InsertNextCell(1);
	drTrajectoryCells_->InsertCellPoint(0);

	drTrajectoryPolydata_= vtkSmartPointer<vtkPolyData>::New();
	drTrajectoryPolydata_->SetPoints(drTrajectoryPoints_);
	drTrajectoryPolydata_->SetLines(drTrajectoryCells_);
	drTrajectoryMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();

#if VTK_MAJOR_VERSION <= 5
	drTrajectoryMapper_->SetInput(drTrajectoryPolydata_);
#else
        drTrajectoryMapper_->SetInputData(drTrajectoryPolydata_);
#endif

	drTrajectoryMapper_->Update();
	drTrajectoryActor_ =  vtkSmartPointer<vtkActor>::New();
	drTrajectoryActor_->SetMapper(drTrajectoryMapper_);
	drTrajectoryActor_->GetProperty()->SetColor(0,255,0);

	estTrajectoryPoints_ = vtkSmartPointer<vtkPoints>::New();
	estTrajectoryPoints_->SetNumberOfPoints(groundtruth_pose.size());


	//gtTrajectoryPolyLine_ = vtkSmartPointer<vtkPolyLine>::New();
	estTrajectoryCells_ = vtkSmartPointer<vtkCellArray>::New();
	estTrajectoryCells_->InsertNextCell(1);
	estTrajectoryCells_->InsertCellPoint(0);

	estTrajectoryPolydata_= vtkSmartPointer<vtkPolyData>::New();
	estTrajectoryPolydata_->SetPoints(estTrajectoryPoints_);
	estTrajectoryPolydata_->SetLines(estTrajectoryCells_);
	estTrajectoryMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	estTrajectoryMapper_->SetInput(estTrajectoryPolydata_);
#else
        estTrajectoryMapper_->SetInputData(estTrajectoryPolydata_);
#endif

	estTrajectoryMapper_->Update();
	estTrajectoryActor_ =  vtkSmartPointer<vtkActor>::New();
	estTrajectoryActor_->SetMapper(estTrajectoryMapper_);
	estTrajectoryActor_->GetProperty()->SetColor(0.9,0.0,0.0);



	measurementPoints_ = vtkSmartPointer<vtkPoints>::New();
	measurementPoints_->SetNumberOfPoints(1);
	measurementCells_ = vtkSmartPointer<vtkCellArray>::New();
	measurementCells_->InsertNextCell(1);
	measurementCells_->InsertCellPoint(0);

	measurementPolydata_= vtkSmartPointer<vtkPolyData>::New();
	measurementPolydata_->SetPoints(measurementPoints_);
	measurementPolydata_->SetLines(measurementCells_);
	measurementMapper_ = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	measurementMapper_->SetInput(measurementPolydata_);
#else
        measurementMapper_->SetInputData(measurementPolydata_);
#endif

	measurementMapper_->Update();
	measurementActor_ =  vtkSmartPointer<vtkActor>::New();
	measurementActor_->SetMapper(measurementMapper_);
	measurementActor_->GetProperty()->SetColor(0.9,0.0,0.9);



	int gtmapSize = groundtruth_landmark.size();


	gtmapColors_->SetNumberOfComponents(3);
	gtmapColors_->SetNumberOfTuples(gtmapSize);
	gtmapPoints_->SetNumberOfPoints(gtmapSize);
	for (int m = 0; m < gtmapSize; m++) {
		MeasurementModel_6D::TLandmark::Vec u;

		groundtruth_landmark[m].get(u);

 		gtmapPoints_->SetPoint(m, u(0), u(1), u(2));
		gtmapColors_->SetTuple3(m,00,191,255);
	}




	renderer_->AddActor(mapActor_);
	renderer_->AddActor(gtmapActor_);
	renderer_->AddActor(particleActor_);
	renderer_->AddActor(gtTrajectoryActor_);
	renderer_->AddActor(estTrajectoryActor_);
	renderer_->AddActor(drTrajectoryActor_);
	renderer_->AddActor(measurementActor_);

	double center[3];
	center[0]=0;
	center[1]=0;
	center[2]=0;
	sphereSource_->SetCenter(center);
	sphereSource_->Update();
	sphereSource_->SetRadius(sphere_radius_);

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
	    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	renderWindowInteractor_->SetInteractorStyle(style);






	renderer_->SetBackground(0, 0, 0); // Background color green
	renderWindow_->Render();
	renderWindowInteractor_->Initialize();
	renderWindowInteractor_->AddObserver(vtkCommand::TimerEvent,this, &Visualizer6D::pause);


}
void Visualizer6D::run(){

	display_mutex_->lock();
	int timerId = renderWindowInteractor_->CreateRepeatingTimer(100);
	renderWindowInteractor_->Start();
/*
	while(!stopped){

		display_mutex_->lock();
		renderWindow_->Render();
		renderWindowInteractor_->Render();

		display_mutex_->unlock();
		std::this_thread::sleep_for(std::chrono::milliseconds(10));

	}
*/

}
void Visualizer6D::pause(){
	display_mutex_->unlock();
	std::this_thread::sleep_for(std::chrono::milliseconds(10));

	display_mutex_->lock();
}

void Visualizer6D::update(RBPHDFilter<MotionModel_Odometry6d, StaticProcessModel<Landmark3d>,
		MeasurementModel_6D,
		KalmanFilter<StaticProcessModel<Landmark3d>, MeasurementModel_6D> > *pFilter_){
	display_mutex_->lock();
	// particles
	int i_w_max = 0;
	double w_max = 0;
	MotionModel_Odometry6d::TState x_i;
	particleColors_->SetNumberOfComponents(3);
	particleColors_->SetNumberOfTuples(pFilter_->getParticleCount());
	particlePoints_->SetNumberOfPoints(pFilter_->getParticleCount());
	for (int i = 0; i < pFilter_->getParticleCount(); i++) {
		x_i = *(pFilter_->getParticleSet()->at(i));
		double w = pFilter_->getParticleSet()->at(i)->getWeight();
		if (w > w_max) {
			i_w_max = i;
			w_max = w;
		}
		particlePoints_->SetPoint(i,x_i.get(0),x_i.get(1),x_i.get(2));
		particleColors_->SetTuple3(i,255,0,0);

	}

	// map

        int gmSize = pFilter_->getGMSize(i_w_max);

	int mapSize = 0;


	mapColors_->SetNumberOfComponents(3);
	mapColors_->SetNumberOfTuples(gmSize);
	mapPoints_->SetNumberOfPoints(gmSize);
	for (int m = 0; m < gmSize; m++) {
		MeasurementModel_6D::TLandmark::Vec u;
		MeasurementModel_6D::TLandmark::Mat S;
		double w;
		pFilter_->getLandmark(i_w_max, m, u, S, w);
		if(w > 0.5){
		  mapPoints_->SetPoint(mapSize, u(0), u(1), u(2));
		  mapColors_->SetTuple3(mapSize,200,200,200);
		  mapSize++;
		}
	}

        mapColors_->SetNumberOfTuples(mapSize);
        mapPoints_->SetNumberOfPoints(mapSize);

	// estimated Trajectory
	x_i = *(pFilter_->getParticleSet()->at(i_w_max));
	if (!init_trajectory){
	  estTrajectoryPoints_->SetPoint(estTrajectoryCells_->GetNumberOfCells()-1 , x_i.get(0),x_i.get(1),x_i.get(2));
	  init_trajectory= true;
	}else{
	  //gtTrajectory

	  gtTrajectoryCells_->InsertNextCell(2);
	  gtTrajectoryCells_->InsertCellPoint(gtTrajectoryCells_->GetNumberOfCells()-1);
	  gtTrajectoryCells_->InsertCellPoint(gtTrajectoryCells_->GetNumberOfCells()-2);
	  gtTrajectoryCells_->Modified();
	  gtTrajectoryPoints_->Modified();

	  drTrajectoryCells_->InsertNextCell(2);
	  drTrajectoryCells_->InsertCellPoint(drTrajectoryCells_->GetNumberOfCells()-1);
	  drTrajectoryCells_->InsertCellPoint(drTrajectoryCells_->GetNumberOfCells()-2);
	  drTrajectoryCells_->Modified();
	  drTrajectoryPoints_->Modified();



	  estTrajectoryCells_->InsertNextCell(2);
	  estTrajectoryCells_->InsertCellPoint(estTrajectoryCells_->GetNumberOfCells()-1);
	  estTrajectoryCells_->InsertCellPoint(estTrajectoryCells_->GetNumberOfCells()-2);


	  estTrajectoryPoints_->SetPoint(estTrajectoryCells_->GetNumberOfCells()-1 , x_i.get(0),x_i.get(1),x_i.get(2));

	  estTrajectoryCells_->Modified();
	  estTrajectoryPoints_->Modified();
	 }
	// measurements

	std::vector<MeasurementModel_6D::TMeasurement> Z = pFilter_->getMeasurements();
	MeasurementModel_6D *pMM = pFilter_->getMeasurementModel();
	measurementPoints_->SetNumberOfPoints(Z.size()+1);
	measurementPoints_->SetPoint(0, x_i.get(0),x_i.get(1),x_i.get(2)); // robotPose
	measurementCells_->Reset();

	for(int i =0; i<Z.size();i++){
		MeasurementModel_6D::TLandmark lm;

		pMM->inverseMeasure( x_i,Z[i],lm);
		measurementPoints_->SetPoint(i+1, lm.get(0),lm.get(1),lm.get(2));
		measurementCells_->InsertNextCell(2);
		measurementCells_->InsertCellPoint(0);
		measurementCells_->InsertCellPoint(i+1);

	}
	measurementCells_->Modified();
	measurementPoints_->Modified();

	renderWindowInteractor_->Render();
	display_mutex_->unlock();





}
void Visualizer6D::update(RBLMBFilter<MotionModel_Odometry6d, StaticProcessModel<Landmark3d>,
                MeasurementModel_6D,
                KalmanFilter<StaticProcessModel<Landmark3d>, MeasurementModel_6D> > *pFilter_){
        display_mutex_->lock();
        // particles
        int i_w_max = 0;
        double w_max = 0;
        MotionModel_Odometry6d::TState x_i;
        particleColors_->SetNumberOfComponents(3);
        particleColors_->SetNumberOfTuples(pFilter_->getParticleCount());
        particlePoints_->SetNumberOfPoints(pFilter_->getParticleCount());
        for (int i = 0; i < pFilter_->getParticleCount(); i++) {
                x_i = *(pFilter_->getParticleSet()->at(i));
                double w = pFilter_->getParticleSet()->at(i)->getWeight();
                if (w > w_max) {
                        i_w_max = i;
                        w_max = w;
                }
                particlePoints_->SetPoint(i,x_i.get(0),x_i.get(1),x_i.get(2));
                particleColors_->SetTuple3(i,255,0,0);

        }

        // map

        int trackNum = pFilter_->getTrackNum(i_w_max);

        int mapSize = 0;


        mapColors_->SetNumberOfComponents(3);
        mapColors_->SetNumberOfTuples(trackNum);
        mapPoints_->SetNumberOfPoints(trackNum);
        for (int m = 0; m < trackNum; m++) {
                MeasurementModel_6D::TLandmark::Vec u;
                MeasurementModel_6D::TLandmark::Mat S;
                GMBernoulliComponent<MeasurementModel_6D::TLandmark> track;
                pFilter_->getTrack(i_w_max, m, track);
                if(track.getP() > 0.1){
                  MeasurementModel_6D::TLandmark *landmark;
                  track.getMaxGaussian(landmark);
                  landmark->get(u,S);
                  mapPoints_->SetPoint(mapSize, u(0), u(1), u(2));
                  mapColors_->SetTuple3(mapSize,200,200,200);
                  mapSize++;
                }



        }

        mapColors_->SetNumberOfTuples(mapSize);
        mapPoints_->SetNumberOfPoints(mapSize);

        // estimated Trajectory
        x_i = *(pFilter_->getParticleSet()->at(i_w_max));
        if (!init_trajectory){
          estTrajectoryPoints_->SetPoint(estTrajectoryCells_->GetNumberOfCells()-1 , x_i.get(0),x_i.get(1),x_i.get(2));
          init_trajectory= true;
        }else{
          //gtTrajectory

          gtTrajectoryCells_->InsertNextCell(2);
          gtTrajectoryCells_->InsertCellPoint(gtTrajectoryCells_->GetNumberOfCells()-1);
          gtTrajectoryCells_->InsertCellPoint(gtTrajectoryCells_->GetNumberOfCells()-2);
          gtTrajectoryCells_->Modified();
          gtTrajectoryPoints_->Modified();

          drTrajectoryCells_->InsertNextCell(2);
          drTrajectoryCells_->InsertCellPoint(drTrajectoryCells_->GetNumberOfCells()-1);
          drTrajectoryCells_->InsertCellPoint(drTrajectoryCells_->GetNumberOfCells()-2);
          drTrajectoryCells_->Modified();
          drTrajectoryPoints_->Modified();



          estTrajectoryCells_->InsertNextCell(2);
          estTrajectoryCells_->InsertCellPoint(estTrajectoryCells_->GetNumberOfCells()-1);
          estTrajectoryCells_->InsertCellPoint(estTrajectoryCells_->GetNumberOfCells()-2);


          estTrajectoryPoints_->SetPoint(estTrajectoryCells_->GetNumberOfCells()-1 , x_i.get(0),x_i.get(1),x_i.get(2));

          estTrajectoryCells_->Modified();
          estTrajectoryPoints_->Modified();
         }
        // measurements

        std::vector<MeasurementModel_6D::TMeasurement> Z = pFilter_->getMeasurements();
        MeasurementModel_6D *pMM = pFilter_->getMeasurementModel();
        measurementPoints_->SetNumberOfPoints(Z.size()+1);
        measurementPoints_->SetPoint(0, x_i.get(0),x_i.get(1),x_i.get(2)); // robotPose
        measurementCells_->Reset();

        for(int i =0; i<Z.size();i++){
                MeasurementModel_6D::TLandmark lm;

                pMM->inverseMeasure( x_i,Z[i],lm);
                measurementPoints_->SetPoint(i+1, lm.get(0),lm.get(1),lm.get(2));
                measurementCells_->InsertNextCell(2);
                measurementCells_->InsertCellPoint(0);
                measurementCells_->InsertCellPoint(i+1);

        }
        measurementCells_->Modified();
        measurementPoints_->Modified();

        renderWindowInteractor_->Render();
        display_mutex_->unlock();





}

}
