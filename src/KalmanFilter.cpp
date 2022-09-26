


  // KalmanFilter Specialization for range bearing measurement model

  // 
  
#include "KalmanFilter.hpp"  
namespace rfs{
  template<>
  KalmanFilter<StaticProcessModel<Landmark2d>, MeasurementModel_RngBrg>::KalmanFilter(){
  config.rangeInnovationThreshold_ = -1;
  config.bearingInnovationThreshold_ = -1;
}
  
  template<>
  KalmanFilter<StaticProcessModel<Landmark2d>, MeasurementModel_RngBrg>::KalmanFilter(StaticProcessModel<Landmark2d> *pProcessModel,
						   MeasurementModel_RngBrg *pMeasurementModel){
    
  pMeasurementModel_= pMeasurementModel;
  pProcessModel_= pProcessModel;
  I.setIdentity();
    
  config.rangeInnovationThreshold_ = -1;
  config.bearingInnovationThreshold_ = -1;
}
  template<>
  bool KalmanFilter<StaticProcessModel<Landmark2d>, MeasurementModel_RngBrg>::calculateInnovation(TMeasurement::Vec &z_exp, TMeasurement::Vec &z_act, TMeasurement::Vec &z_innov){
    
  z_innov = z_act - z_exp;
  
  if(config.rangeInnovationThreshold_ > 0 && fabs(z_innov(0)) > config.rangeInnovationThreshold_)
    return false;
  while(z_innov(1) > PI)
    z_innov(1) -= 2 * PI;
  while(z_innov(1) < -PI)
    z_innov(1) += 2 * PI;
  if(config.bearingInnovationThreshold_ > 0 && fabs(z_innov(1)) > config.bearingInnovationThreshold_ )
    return false;
  return true;
}


}
