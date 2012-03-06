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

#include "rect_constraint.h"
#include "WMG.h"



/****************************************
 * TYPEDEFS 
 ****************************************/


/// @addtogroup gWMG_INTERNALS
/// @{

/** \brief Defines a footstep. */
class footstep : public RectangularConstraint_ZMP
{
    public:
        footstep (
                const double, 
                const Transform<double, 3>&,
                const Vector3d&,
                const unsigned int, 
                const fs_type, 
                const double *);
        footstep (const footstep&);
        ~footstep();

        void changePosture(const double *, const bool);
        double x();
        double y();



        /// Angle (relative to the world frame) of a footstep [rad.].
        double angle;

        /// cos(angle).
        double ca;

        /// sin(angle).
        double sa;

        /// the period of time spent in this support
        unsigned int time_period;

        /// the amount of time left in this support (=time_period on initialization)
        unsigned int time_left;

        /// type of the step.
        fs_type type;

        /// Reference ZMP point
        Vector3d ZMPref;


        Transform<double, 3>* posture;
};
///@}
#endif /*FOOTSTEP_H*/
