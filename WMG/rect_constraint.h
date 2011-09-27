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
 * DEFINES
 ****************************************/



/****************************************
 * TYPEDEFS 
 ****************************************/
/// @addtogroup gWMG_INTERNALS
/// @{

/** \brief Defines rectangular constraints (of the form D*z <= d) for
    the ZMP.

    \note D is a [4 x 2] matrix stored column-wise (Fortran style). 
    D is always initialized as (later it could be rotated using D*R'):
    \verbatim
    D[0] =  1.0; D[4] =  0.0;
    D[1] =  0.0; D[5] =  1.0;
    D[2] = -1.0; D[6] =  0.0;
    D[3] =  0.0; D[7] = -1.0;
    \endverbatim
*/
class RectangularConstraint_ZMP
{
    public:
        RectangularConstraint_ZMP();
        void set_size(double * _d);
        void print();
        void rotate_translate(double ca, double sa, Point2D p);
        void Constraints2Vert();



        /** \brief Matrix of the constraints D*z <= d (where z is a 2D point). This variable is formed
            internaly.
            
            \note D is a [4 x 2] matrix stored column-wise (Fortran style). D takes the form 
            \verbatim
            D = [ 1  0;
                  0  1;
                 -1  0;
                  0 -1] * [cos(a) -sin(a)
                           sin(a)  cos(a)].
            \endverbatim
         */
        double D[4*2];

        /** \brief Size of the support polygon for a single support (USER INPUT)
                
            \verbatim
             ------------------------------------------------
             |                    |                         |
             |                    |d0(2)                    |
             |                    |                         |
             |      d0(3)         |          d0(1)          |
             |------------------- p ------------------------|
             |                    |                         |
             |                    |                         |
             |                    |d0(4)                    |
             |                    |                         |
             ------------------------------------------------
            \endverbatim
         */
        double d0[4];

        /** \brief Vector of the constraints D*z <= d (where z is a 2D point). This variable is formed internaly. 
                
            \note d is a [4 x 1] vector. d takes the form #d = #d0 - #D*p. p is a 2D reference point
            (expressed with respect to the world frame) for the rectangular constraint (see the sketch for
            #d0).
                
         */
        double d[4];
         
        std::vector<Point2D> vert;
};


/****************************************
 * PROTOTYPES 
 ****************************************/
///@}
#endif /*RECT_CONSTRAINT_H*/

