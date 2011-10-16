/** 
 * @file
 * @author Alexander Sherikov
 * @date 15.10.2011 23:38:54 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "problem_param.h"

#include <cmath> //cos,sin

/****************************************
 * FUNCTIONS 
 ****************************************/

problem_parameters::problem_parameters (
    const int N_,
    const double Alpha,
    const double Beta,
    const double Gamma,
    const double regularization)
{
    N = N_;

    i2Q[0] = 1/(2*(Beta/2));
    i2Q[1] = 1/(2*(Alpha/2));
    i2Q[2] = 1/(2*regularization);

    i2P = 1/(2 * (Gamma/2));

    spar = new state_parameters[N];
}



problem_parameters::~problem_parameters()
{
    if (spar != NULL)
        delete [] spar;
}



/** @brief Initializes quadratic problem.
    @param[in] T_ Sampling time (for the moment it is assumed to be constant) [sec.]
    @param[in] h_ Height of the Center of Mass divided by gravity
    @param[in] angle Rotation angle for each state in the preview window
 */
void problem_parameters::set_state_parameters (
    const double* T_,
    const double* h_,
    const double* angle)
{
    int i;

#ifndef SMPC_VARIABLE_T_h
    T = T_[0];
    h = h_[0];

    A6 = B[1] = T*T/2;
    B[0] = B[1]*T/3 - h*T;
#endif

    for (i = 0; i < N; i++)
    {
        spar[i].cos = cos(angle[i]);
        spar[i].sin = sin(angle[i]);

#ifdef SMPC_VARIABLE_T_h
        spar[i].T = T_[i];
        spar[i].h = h_[i];

        spar[i].B[1] = T_[i]*T_[i]/2;
        spar[i].B[0] = spar[i].B[1]*T_[i]/3 - h_[i]*T_[i];

        if (i == 0)
        {
            /// @todo We need delta_h here.
            spar[i].A6 = T_[i]*T_[i]/2;
        }
        else
        {
            spar[i].A6 = T_[i]*T_[i]/2 - (h_[i] - h_[i-1]);
        }
#endif
    }
}
