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


/****************************************
 * DEFINES
 ****************************************/



/****************************************
 * TYPEDEFS 
 ****************************************/

/// @addtogroup gWMG_INTERNALS
/// @{

/** \brief Defines a footstep. */
class FootStep
{
    public:
        FootStep();
        FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL);
        FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL, double *_d);
        void set(double _angle, Point2D _p, int _nSS, int _n, int _RL);
        void print();



        /** \brief Position (in the world frame) of a footstep [meter]. */
        Point2D p;

        /** \brief Angle (relative to the world frame) of a footstep [rad.]. */
        double angle;

        /** \brief cos(angle). */
        double ca;

        /** \brief sin(angle). */
        double sa;

        /** \brief Constraints defining the Polyogn of Support (PoS). */
        RectangularConstraint_ZMP ctr;

        /** \brief Offset from the "point of interest of a constraint" [meter] defined in the local frame. 
                
            \note This is used in order to define ZMP_ref
         */
        Point2D Offset;

        /** \brief If RL = 1 the left foot is in support. If RL = -1 the right foot is in support. */
        int RL;

        /** \brief Number of (preview window) iterations in Single Support. */
        int nSS;

        /** \brief Total number of (preview window) iterations, i.e., nSS + nDS. */
        int n;
};


/****************************************
 * PROTOTYPES 
 ****************************************/

///@}
#endif /*FOOTSTEP_H*/
