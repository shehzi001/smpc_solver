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
namespace AS
{
    class state_parameters
    {
        public:
            // parameters used in generation of A and B matrices
            /** Preview sampling time  */
            double T;
            /** h = @ref ph "hCoM/gravity". */
            double h;

            double A3;
            double A6;

            double B[3];
    };


    /**
     * @brief A set of problem parameters.
     */
    class problem_parameters
    {
        public:
            problem_parameters (const int, const double, const double, const double, const double);
            ~problem_parameters();

            void set_state_parameters (const double*, const double*, const double);


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

            /// Height of the CoM at initial state divided by the gravity, this initial state
            /// precede the first state in the preview window.
            double h_initial;
    };
}
///@}
#endif /*PROBLEM_PARAM_H*/

