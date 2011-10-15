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

    angle_cos = new double[N];
    angle_sin = new double[N];
    T = NULL;
    h = NULL;

#ifdef SMPC_VARIABLE_T_h
    A6 = new double[N-1];
    B = new double[2*N];
#endif
}



problem_parameters::~problem_parameters()
{
    if (angle_cos != NULL)
        delete angle_cos;

    if (angle_sin != NULL)
        delete angle_sin;

#ifdef SMPC_VARIABLE_T_h
    if (B != NULL)
        delete B;
    if (A6 != NULL)
        delete A6;
#endif
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

    T = T_;
    h = h_;

    for (i = 0; i < N; i++)
    {
        angle_cos[i] = cos(angle[i]);
        angle_sin[i] = sin(angle[i]);
    }

#ifdef SMPC_VARIABLE_T_h
    for (i = 0; i < N; i++)
    {
        B[i*2+1] = T[i]*T[i]/2;
        B[i*2] = B[i*2+1]*T[i]/3 - h[i]*T[i];
    }
    for (i = 0; i < N-1; i++)
    {
        A6[i] = T[i+1]*T[i+1]/2 - (h[i+1] - h[i]);
    }
#else
    A6 = B[1] = T[0]*T[0]/2;
    B[0] = B[1]*T[0]/3 - h[0]*T[0];
#endif
}
