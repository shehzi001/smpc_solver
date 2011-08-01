/**
 * @file
 * @brief 
 *
 * @author Alexander Sherikov
 * @date 19.07.2011 22:37:28 MSD
 *
 * @todo add description
 *
 * @todo add doxygen config
 * @todo add tests
 * @todo add config.h
 */


#ifndef SMPC_COMMON_H
#define SMPC_COMMON_H

/****************************************
 * INCLUDES 
 ****************************************/

#include <cstddef>


/****************************************
 * DEFINES
 ****************************************/
#define NUM_STATE_VAR 6
#define NUM_CONTROL_VAR 2
#define NUM_VAR 8

#define QPAS_VARIABLE_T_h

/****************************************
 * TYPEDEFS 
 ****************************************/

/**
 * @brief 
 * @todo This is ugly. What to do?
 */
struct chol_solve_param 
{
// static matrices and vectors
    /** State related penalty.*/
    double i2Q[3];

    /** Control related penalty. */
    double iP;

    /** inv(H) * g */
    double *iHg;


    double *angle_cos;
    double *angle_sin;


// parameters used in generation of A and B matrices
#ifdef QPAS_VARIABLE_T_h
    double *T;
    double *h;
    double *dh;
#else
    /** Preview sampling time  */
    double T;

    /** h = #hCoM/#gravity. */
    double h;
#endif
};

/****************************************
 * PROTOTYPES 
 ****************************************/

#endif /*SMPC_COMMON_H*/

