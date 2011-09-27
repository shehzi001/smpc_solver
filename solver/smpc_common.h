/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:37:28 MSD
 */


#ifndef SMPC_COMMON_H
#define SMPC_COMMON_H

/****************************************
 * INCLUDES 
 ****************************************/

#include <cstddef>
#include "common_const.h"

/****************************************
 * DEFINES
 ****************************************/

/// @addtogroup gINTERNALS
/// @{

/// The size of one dimension of 3x3 matrix.
#define MATRIX_DIM_SIZE 3
/// The number of elements in 3x3 matrix.
#define MATRIX_SIZE 9

/// Allow variable chol_solve_param#T and chol_solve_param#h
#define QPAS_VARIABLE_T_h
/// Allow removal of constraints from active set
#define QPAS_DOWNDATE


/****************************************
 * TYPEDEFS 
 ****************************************/

/**
 * @brief A set of parameters used by #chol_solve class.
 */
struct chol_solve_param 
{
// static matrices and vectors
    /** State related penalty.*/
    double i2Q[3];

    /** Control related penalty. */
    double i2P;

    /** inv(H) * g */
    double *iHg;


    double *angle_cos;
    double *angle_sin;


// parameters used in generation of A and B matrices
    /** Preview sampling time  */
    double *T;
    /** h = @ref ph "hCoM/gravity". */
    double *h;

#ifdef QPAS_VARIABLE_T_h
    double *dh;
#endif
};

/****************************************
 * PROTOTYPES 
 ****************************************/
///@}
#endif /*SMPC_COMMON_H*/

