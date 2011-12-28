/**
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:35:26 MSD
 */


#ifndef RECT_CONSTRAINT_H
#define RECT_CONSTRAINT_H

/****************************************
 * INCLUDES 
 ****************************************/

#include <vector>

#include "point2d.h"



/****************************************
 * TYPEDEFS 
 ****************************************/
/// @addtogroup gWMG_INTERNALS
/// @{

/// \brief Defines rectangular constraints (of the form D*z <= d) for the ZMP.
class RectangularConstraint_ZMP
{
    public:
        RectangularConstraint_ZMP(const double *);
        void rotate_translate(const double, const double, const Point2D&);
        void Constraints2Vert();



        /** 
         * \brief Matrix of the constraints D*z <= d (where z is a 2D point). 
         *
         * \note D is a [4 x 2] matrix stored column-wise (Fortran style). 
         * D is always initialized as (later it could be rotated using D*R'):
         * \verbatim
            D[0] =  1.0; D[4] =  0.0;
            D[1] =  0.0; D[5] =  1.0;
            D[2] = -1.0; D[6] =  0.0;
            D[3] =  0.0; D[7] = -1.0;

            D = [ 1  0;
                  0  1;
                 -1  0;
                  0 -1] * [cos(a) -sin(a)
                           sin(a)  cos(a)].
           \endverbatim
         */
        double D[4*2];

        /** Size of the support polygon for a single support.
            Vector d of the constraints D*z <= d (where z is a 2D point).
         */
        double d[4];

        /// Size of the support polygon for a single support (no rotation / translation).
        double d_orig[4];

        /// Absolute coordinates of vertices.
        std::vector<Point2D> vert;
};

///@}
#endif /*RECT_CONSTRAINT_H*/

