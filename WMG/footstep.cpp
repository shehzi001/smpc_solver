/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:41:58 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include <cmath> // sin, cos
#include <cstdio> // printf

#include "footstep.h"


/****************************************
 * FUNCTIONS 
 ****************************************/

/** \brief Default constructor. */
FootStep::FootStep()
{
    p.set(0.0, 0.0);
    angle = 0.0; ca = 0.0; sa = 0.0;
    Offset.set(0.0, 0.0);
    nSS = 0;
    n = 0;
    RL = 0;
}



/** \brief Defines a footstep at a given position with a given orientation. 
        
    \param[in] _p 2D position (of the "point of interest")
    \param[in] _angle orientation
    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
 */
FootStep::FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL)
{
    set(_angle, _p, _nSS, _n, _RL);
    Offset.set(0.0, 0.0);
}



/** \brief Defines a footstep at a given position with a given orientation. 
        
    \param[in] _p 2D position (of the "point of interest")
    \param[in] _angle orientation
    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
    \param[in] _d Vector of the PoS constraints (assumed to be [4 x 1])
 */
FootStep::FootStep(double _angle, Point2D _p, int _nSS, int _n, int _RL, double *_d)
{
    ctr.set_size(_d);
    set(_angle, _p, _nSS, _n, _RL);
    Offset.set(0.0, 0.0);
}



/** \brief Sets the position and orientation of a footstep.
        
    \param[in] _p 2D Position (of the point of interest)
    \param[in] _angle Orientation
    \param[in] _nSS Number of (preview window) iterations in Single Support.
    \param[in] _n Total number of (preview window) iterations, i.e., nSS + nDS.
    \param[in] _RL Right or left foot (1 for left, -1 for right).
 */
void FootStep::set(double _angle, Point2D _p, int _nSS, int _n, int _RL)  
{
    p.set(_p.x, _p.y);
    angle = _angle; ca = cos(angle); sa = sin(angle);
    ctr.rotate_translate(ca, sa, p);
    nSS = _nSS;
    n = _n;
    RL = _RL;
}



/** \brief Prints information about a footstep. */
void FootStep::print()
{
    printf("----------------------------------\n");
    printf(" nSS = %d, n = %d, RL = %d  \n", nSS, n, RL);
    printf(" angle = % f, \n", angle);
    printf(" Offset = [% f, %f] \n", Offset.x, Offset.y);
    p.print();
    ctr.print();
    printf("----------------------------------\n");
}
