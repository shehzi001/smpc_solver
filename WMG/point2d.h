/**
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 16:02:38 MSD
 */


#ifndef POINT2D_H
#define POINT2D_H

/****************************************
 * INCLUDES 
 ****************************************/



/****************************************
 * TYPEDEFS 
 ****************************************/
/// @addtogroup gWMG_INTERNALS
/// @{

/** \brief Defines a 2D point. */
class Point2D
{
    public:
        Point2D();
        Point2D(double _x, double _y);
        void print();
        void set(double _x, double _y);



        /** \brief x position [meter]*/
        double x;

        /** \brief y position [meter] */
        double y;
};
///@}
#endif /*POINT2D_H*/
