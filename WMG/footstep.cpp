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
 * @param[in] Position absolute position of the foot.
 * @param[in] ZMPref_ absolute reference ZMP position for the foot.
 * @param[in] repeat_times_ number of times this step appears in the preview window.
 * @param[in] type_ type of the step.
 * @param[in] d_ ZMP constraints as defined in RectangularConstraint_ZMP#RectangularConstraint_ZMP.
 */
FootStep::FootStep(
        const double angle_, 
        const Point2D Position,
        const Point2D ZMPref_,
        const int repeat_times_, 
        const fs_type type_, 
        const double *d_) : 
    Point2D(Position), 
    RectangularConstraint_ZMP(d_),
    ZMPref(ZMPref_)
{
    type = type_;
    angle = angle_; 
    ca = cos(angle); 
    sa = sin(angle);
    rotate_translate(ca, sa, *this);
    repeat_times = repeat_times_;
    repeat_counter = repeat_times;
}
