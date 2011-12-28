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

/**
 * @brief A constructor.
 *
 * @param[in] d_ a [4x1] array of constraints.
 *
 * The constraints are 4 positive numbers, their meaning is illustrated
 * by the following picture:
    \verbatim
     ------------------------------------------------
     |                    |                         |
     |                    |d(2)                     |
     |                    |                         |
     |      d(3)          |          d(1)           |
     |------------------- p ------------------------|
     |                    |                         |
     |                    |                         |
     |                    |d(4)                     |
     |                    |                         |
     ------------------------------------------------
    \endverbatim
    At creation p is assumed to be [0;0] and rotation angle = 0.
 */
RectangularConstraint_ZMP::RectangularConstraint_ZMP(const double *d_)
{
    D[0] =  1.0; D[4] =  0.0;
    D[1] =  0.0; D[5] =  1.0;
    D[2] = -1.0; D[6] =  0.0;
    D[3] =  0.0; D[7] = -1.0;

    d[0] = d_[0]; 
    d[1] = d_[1]; 
    d[2] = d_[2]; 
    d[3] = d_[3];

    d_orig[0] = d[0];
    d_orig[1] = d[1];
    d_orig[2] = d[2];
    d_orig[3] = d[3];

    Constraints2Vert(); // determine coordinates of vertices
}



/** 
 * \brief translates from [0;0] to "p" and rotates from 0 to "angle". The 0 
 * initial angle implies that D = [eye(2); -eye(2)].

    \param[in] ca cos(angle)
    \param[in] sa sin(angle)
    \param[in] p a 2D reference point 
                
    \note This is used when the constraints are initialized (only then the orientation can be set).
 */
void RectangularConstraint_ZMP::rotate_translate(const double ca, const double sa, const Point2D& p)
{
    // D = D*R'
    D[0] =  ca; D[4] =  sa;
    D[1] = -sa; D[5] =  ca;
    D[2] = -ca; D[6] = -sa;
    D[3] =  sa; D[7] = -ca;

    // d = d - D*p
    d[0] += D[0]*p.x + D[4]*p.y;
    d[1] += D[1]*p.x + D[5]*p.y;
    d[2] += D[2]*p.x + D[6]*p.y;
    d[3] += D[3]*p.x + D[7]*p.y;

    Constraints2Vert(); // determine coordinates of vertices
}



/** \brief Computes the vertices of a polygon from the constraints (D*x <= d)
        
    The numbering of the vertices is 
    \verbatim
    3---------2
    |         |
    |         |
    4---------1
    \endverbatim
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

    //   D     d       inv (D)
    // |0 4|   0    1/det * | 7 -4|
    // |3 7|   3            |-3  0|
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

