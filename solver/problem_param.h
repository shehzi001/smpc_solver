/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:37:28 MSD
 */


#ifndef PROBLEM_PARAM_H
#define PROBLEM_PARAM_H

/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"

/****************************************
 * DEFINES
 ****************************************/

/****************************************
 * TYPEDEFS 
 ****************************************/

/**
 * @brief A set of problem parameters.
 */
class problem_parameters
{
    public:
        problem_parameters (const int, const double, const double, const double, const double);
        ~problem_parameters();

        void set_state_parameters (const double*, const double*, const double*);


        /** Number of iterations in a preview window. */
        int N;

    // static matrices and vectors
        ///@{
        /** State related penalty.*/
        double i2Q[3];
        ///@}

        ///@{
        /** Control related penalty. */
        double i2P;
        ///@}


        double *angle_cos;
        double *angle_sin;


    // parameters used in generation of A and B matrices
        /** Preview sampling time  */
        const double *T;
        /** h = @ref ph "hCoM/gravity". */
        const double *h;

#ifdef SMPC_VARIABLE_T_h
        double *B;
        double *A6;
#else
        double B[3];
        double A6;
#endif
};

/****************************************
 * PROTOTYPES 
 ****************************************/
///@}
#endif /*PROBLEM_PARAM_H*/

