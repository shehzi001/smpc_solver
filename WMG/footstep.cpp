/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:41:58 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include <cmath> // sin, cos

#include "footstep.h"


/****************************************
 * FUNCTIONS 
 ****************************************/

/**
 * @brief Defines a footstep at a given position with a given orientation.
 *
 * @param[in] angle_ absolute rotation angle.
 * @param[in] posture_ absolute position and orientation of the foot.
 * @param[in] ZMPref_ absolute reference ZMP position for the foot.
 * @param[in] time_period_ amount of time to spend in the step (ms.).
 * @param[in] type_ type of the step.
 * @param[in] d_ ZMP constraints as defined in RectangularConstraint_ZMP#RectangularConstraint_ZMP.
 */
footstep::footstep(
        const double angle_, 
        const Transform<double, 3>& posture_,
        const Vector3d& ZMPref_,
        const unsigned int time_period_, 
        const fs_type type_, 
        const double *d_) : 
    RectangularConstraint_ZMP(d_),
    ZMPref(ZMPref_)
{
    posture = posture_;
    type = type_;
    angle = angle_; 
    ca = cos(angle); 
    sa = sin(angle);
    rotate_translate(ca, sa, x(), y());
    time_left = time_period = time_period_;
}


/**
 * @return x coordinate
 */
double footstep::x()
{
    return (posture.translation()[0]);
}


/**
 * @return y coordinate
 */
double footstep::y()
{
    return (posture.translation()[1]);
}



/**
 * @brief Correct position of the footstep.
 *
 * @param[in] new_posture new posture of the step.
 */
void footstep::changePosture (const double * new_posture)
{
    posture.matrix() = Matrix4d::Map(new_posture);
    Matrix3d rotation = posture.matrix().corner(TopLeft,3,3);
    angle = rotation.eulerAngles(0,1,2)[2];
    ca = cos(angle);
    sa = sin(angle);
    rotate_translate(ca, sa, x(), y());
}
