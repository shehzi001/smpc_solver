/** 
 * @file
 * @author Alexander Sherikov
 * @date 27.09.2011 18:37:55 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "rect_constraint.h"

/****************************************
 * FUNCTIONS 
 ****************************************/

/** \brief Default constructor.
    \brief NAO constraint with safety margin.
 */
RectangularConstraint_ZMP::RectangularConstraint_ZMP()
{
    D[0] =  1.0; D[4] =  0.0; d0[0] = 0.09;  d[0] = 0.09; 
    D[1] =  0.0; D[5] =  1.0; d0[1] = 0.025; d[1] = 0.025;
    D[2] = -1.0; D[6] =  0.0; d0[2] = 0.03;  d[2] = 0.03;
    D[3] =  0.0; D[7] = -1.0; d0[3] = 0.025; d[3] = 0.025;

    Constraints2Vert();
}



/** \brief Set the foot size.
    \param[in] _d Vector of the PoS constraints (assumed to be [4 x 1]).
 */
void RectangularConstraint_ZMP::set_size(double * _d)
{
    d0[0] = _d[0]; 
    d0[1] = _d[1]; 
    d0[2] = _d[2]; 
    d0[3] = _d[3];
    
    Constraints2Vert();
}



/** \brief Prints the constraints. */
void RectangularConstraint_ZMP::print()
{
    printf("           D                 d \n");
    printf(" |% f % f|  % f\n", D[0], D[4], d[0]);
    printf(" |% f % f|  % f\n", D[1], D[5], d[1]);
    printf(" |% f % f|  % f\n", D[2], D[6], d[2]);
    printf(" |% f % f|  % f\n\n", D[3], D[7], d[3]);
    printf("          vert \n");
    for (int i=0; i<4; i++)
    {
        printf(" % f % f \n", vert[i].x, vert[i].y);
    }
}



/** \brief translates from [0;0] to "p" and rotates from 0 to "angle". The 0 initial angle implies
    that D = [eye(2); -eye(2)].

    \param[in] ca cos(angle)
    \param[in] sa sin(angle)
    \param[in] p 2D point 

    \note This is used when the constraints are initialized (only then the orientation can be set).
 */
void RectangularConstraint_ZMP::rotate_translate(double ca, double sa, Point2D p)
{
    // D = D*R'
    D[0] =  ca; D[4] =  sa;
    D[1] = -sa; D[5] =  ca;
    D[2] = -ca; D[6] = -sa;
    D[3] =  sa; D[7] = -ca;
    
    // d = d0 - D*p
    d[0] = d0[0] + D[0]*p.x + D[4]*p.y;
    d[1] = d0[1] + D[1]*p.x + D[5]*p.y;
    d[2] = d0[2] + D[2]*p.x + D[6]*p.y;
    d[3] = d0[3] + D[3]*p.x + D[7]*p.y;

    Constraints2Vert();
}



/** \brief Computes the vertices of a polygon from the constraints (D*x <= d)
        
    \note It is assumed that D is a [4 x 2] matrix stored column-wise.
    
    The numbering of the vertices is (see the definition of #d0)
    \verbatim
    3---------2
    |         |
    |         |
    4---------1
    \endverbatim
    
    \note I am too lazy to check whether det = 0 :).
*/
void RectangularConstraint_ZMP::Constraints2Vert()
{
    double det;
    vert.clear();
    
    // Indexes of ctr.D
    // 0 4
    // 1 5 
    // 2 6
    // 3 7
    
    // |0 4|   0     | 7 -4|
    // |3 7|   3     |-3  0|
    det = D[0]*D[7] - D[3]*D[4];
    vert.push_back(Point2D( D[7]/det*d[0] - D[4]/det*d[3],-D[3]/det*d[0] + D[0]/det*d[3])); 
    
    // |0 4|   0     | 5 -4|
    // |1 5|   1     |-1  0|
    det = D[0]*D[5] - D[4]*D[1]; 
    vert.push_back(Point2D( D[5]/det*d[0] - D[4]/det*d[1],-D[1]/det*d[0] + D[0]/det*d[1])); 
    
    // |1 5|   1     | 6 -5|
    // |2 6|   2     |-2  1|
    det = D[1]*D[6] - D[5]*D[2]; 
    vert.push_back(Point2D( D[6]/det*d[1] - D[5]/det*d[2],-D[2]/det*d[1] + D[1]/det*d[2])); 
    
    // |2 6|   2     | 7 -6|
    // |3 7|   3     |-3  2|
    det = D[2]*D[7] - D[3]*D[6]; 
    vert.push_back(Point2D( D[7]/det*d[2] - D[6]/det*d[3],-D[3]/det*d[2] + D[2]/det*d[3])); 
}

