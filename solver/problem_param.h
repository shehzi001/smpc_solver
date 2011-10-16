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

class state_parameters
{
    public:
        double cos;
        double sin;

#ifdef SMPC_VARIABLE_T_h
        // parameters used in generation of A and B matrices
        /** Preview sampling time  */
        double T;
        /** h = @ref ph "hCoM/gravity". */
        double h;

        double A6;

        double B[2];
#endif
};


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

        state_parameters *spar;

#ifndef SMPC_VARIABLE_T_h
        double T;
        double h;
        double B[2];
        double A6;
#endif
};

/****************************************
 * PROTOTYPES 
 ****************************************/
///@}
#endif /*PROBLEM_PARAM_H*/

