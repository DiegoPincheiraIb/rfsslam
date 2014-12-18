
#ifndef FRAME_HPP
#define FRAME_HPP

#include <boost/shared_ptr.hpp>
#include "Pose.hpp"

namespace rfs{

  /** 
   * \class Frame2d
   * \brief A 2d reference frame.
   *
   * Given that the relative base frame is frame b, 
   * and the current frame is frame c,
   * this class contains the rotation information R_0_1, and translation 
   * t_1_0_0 (displacement from 1 to 0 in the frame of 0). If F_b points
   * to NULL, the base frame is assumed to be the inertial frame.
   */
  class Frame2d : public Pose<3, 2, 1>
  {

  public:

    /** \brief Smart pointer to Frame2d */
    typedef boost::shared_ptr<Frame2d> FramePtr;
    
    /** \brief 2d Rotation matrix */ 
    typedef Eigen::Matrix<double, 2, 2> RotMat;
    
    /** \brief Eigen vector of size 2 */
    typedef Pose<3, 2, 1>::PosVec PosVec;

    /** \brief RandomVec of size 2 */
    typedef RandomVec<2> PosRandomVec;

    /** \brief Eigen vector of size 3 */
    typedef Pose<3, 2, 1>::Vec PoseVec;

    Frame2d();
    Frame2d(PoseVec const &T_c_b, TimeStamp const &t = TimeStamp(), FramePtr F_b = FramePtr() );
    Frame2d(Pose<3, 2, 1> T_c_b, FramePtr F_b = FramePtr() );
    ~Frame2d();

    Frame2d operator*(Frame2d const &F_c_d);

    void setBaseFrame(FramePtr const &F_b);
    
    void getRotMat( RotMat &R_b_c ) const;

    void getRelToBaseFrame(PosVec &p_b, PosRandomVec::Cov *p_b_cov = NULL) const;
    void getRelToBaseFrame(PosVec const &p_c, PosVec &p_b,
			   PosRandomVec::Cov *p_c_cov = NULL, PosRandomVec::Cov *p_b_cov = NULL) const;
    void getRelToBaseFrame(PosRandomVec &p_b) const;
    void getRelToBaseFrame(PosRandomVec const &p_c, PosRandomVec &p_b) const;

  private:

    FramePtr F_b_;
    

  };
 
}





#endif
