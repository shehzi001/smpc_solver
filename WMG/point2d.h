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
        /**
         * @brief A constructor.
         *
         * @param[in] x_ x coordinate.
         * @param[in] y_ y coordinate.
         */
        Point2D(const double x_, const double y_)
        {
            x = x_;
            y = y_;
        };

    
        /**
         * @brief New point is located at 'offset' from 'base' in the direction
         *  specified using cos and sin of the rotation angle.
         *
         * @param[in] base base point
         * @param[in] offset offset from the base point
         * @param[in] sa sin of the rotation angle
         * @param[in] ca cos of the rotation angle
         */
        Point2D(const Point2D& base, 
                const Point2D& offset, 
                const double sa, 
                const double ca)
        {
            x = base.x + offset.x*ca - offset.y*sa;
            y = base.y + offset.x*sa + offset.y*ca;
        };


        /** \brief x position [meter]*/
        double x;

        /** \brief y position [meter] */
        double y;
};
///@}
#endif /*POINT2D_H*/
