/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 16:06:46 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include <cstdio> // printf

#include "point2d.h"

/****************************************
 * FUNCTIONS 
 ****************************************/

Point2D::Point2D() {}

Point2D::Point2D(double _x, double _y)  
{
    set(_x, _y);
}

void Point2D::print()
{
    printf(" p = [% f, % f] \n\n", x, y);
}

void Point2D::set(double _x, double _y)
{
    x = _x;
    y = _y;
}
