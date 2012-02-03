/**
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:40:40 MSD
 */


#ifndef FOOTSTEP_H
#define FOOTSTEP_H

/****************************************
 * INCLUDES 
 ****************************************/

#include "point2d.h"
#include "rect_constraint.h"
#include "WMG.h"



/****************************************
 * TYPEDEFS 
 ****************************************/

/// @addtogroup gWMG_INTERNALS
/// @{

/** \brief Defines a footstep. */
class FootStep : public Point2D, public RectangularConstraint_ZMP
{
    public:
        FootStep(
                const double, 
                const Point2D&,
                const Point2D&,
                const int, 
                const fs_type, 
                const double *);

        void correct(const double, const double);

        /// Angle (relative to the world frame) of a footstep [rad.].
        double angle;

        /// cos(angle).
        double ca;

        /// sin(angle).
        double sa;

        /// number of times this step appears in the preview window.
        int repeat_times;

        /// a counter of repeats,on initialization = repeat_times.
        int repeat_counter;

        /// type of the step.
        fs_type type;

        /// Reference ZMP point
        Point2D ZMPref;
};
///@}
#endif /*FOOTSTEP_H*/
