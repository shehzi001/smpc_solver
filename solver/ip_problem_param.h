/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:37:28 MSD
 */


#ifndef IP_PROBLEM_PARAM_H
#define IP_PROBLEM_PARAM_H

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
namespace IP
{
    class state_parameters
    {
        public:
            double cos;
            double sin;

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

            void set_state_parameters (const double*, const double*, const double, const double*);


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
    };
}
///@}
#endif /*IP_PROBLEM_PARAM_H*/

