/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/
#include "qp_solver.h"
#include "state_handling.h"

#include <string.h> // in order to use memcpy and memmove
#include <cmath> // cos, sin

/****************************************
 * FUNCTIONS
 ****************************************/


/** @brief Constructor: initialization of the constant parameters

    @param[in] N_ Number of sampling times in a preview window
    @param[in] Alpha Velocity gain
    @param[in] Beta Position gain
    @param[in] Gamma Jerk gain
    @param[in] regularization regularization
    @param[in] tol_ tolerance
*/
qp_solver::qp_solver(
        const int N_, 
        const double Alpha, 
        const double Beta, 
        const double Gamma, 
        const double regularization, 
        const double tol_) :
    problem_parameters (N_, Alpha, Beta, Gamma, regularization) 
{
    tol = tol_;

    gain_alpha = Alpha;
    gain_beta  = Beta;
    gain_gamma = Gamma;

    alpha = 1;
    
    dX = new double[NUM_VAR*N]();
}


/** @brief Initializes quadratic problem.

    @param[in] T_ Sampling time (for the moment it is assumed to be constant) [sec.]
    @param[in] h_ Height of the Center of Mass divided by gravity
    @param[in] angle Rotation angle for each state in the preview window
*/
void qp_solver::set_state_parameters(
        const double* T_, 
        const double* h_, 
        const double* angle)
{
    T = T_;
    h = h_;
#ifdef SMPC_VARIABLE_T_h
    for (int i = 0; i < N-1; i++)
    {
        dh[i] = h[i+1] - h[i];
    }
#endif

    for (int i = 0; i < N; i++)
    {
        angle_cos[i] = cos(angle[i]);
        angle_sin[i] = sin(angle[i]);
    }
}



/**
 * @brief Generates an initial feasible point. 
 * First we perform a change of variable to @ref pX_tilde "X_tilde"
 * generate a feasible point, and then we go back to @ref pX_bar "X_bar".
 *
 * @param[in] x_coord x coordinates of points satisfying constraints
 * @param[in] y_coord y coordinates of points satisfying constraints
 * @param[in] X_tilde current state
 * @param[in,out] X_ initial guess / solution of optimization problem
 */
void qp_solver::form_init_fp (
        const double *x_coord, 
        const double *y_coord, 
        const double *X_tilde,
        double* X_)
{
    X = X_;

    double *control = &X[NUM_STATE_VAR*N];
    double *cur_state = X;
    const double *prev_state = X_tilde;

#ifndef SMPC_VARIABLE_T_h    
    //------------------------------------
    double Ti = T[0];
    double T2 = Ti*Ti/2;
    double hi = h[0];


    /* Control matrix. */
    double B[3] = {T2*Ti/3 - hi*Ti, T2, Ti};


    /* inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if T^3/6-h*T is
     * not equal to zero). The two elements one the main diagonal are equal, and only one of them 
     * is stored, which is equal to
        1/(T^3/6 - h*T)
     */
    double iCpB = 1/(B[0]);


    /* inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
        iCpB_CpA = [a b c 0 0 0;
                    0 0 0 a b c];

        a = iCpB
        b = iCpB*T
        c = iCpB*T^2/2
     * Only a,b and c are stored.
     */
    double iCpB_CpA[3] = {iCpB, iCpB*Ti, iCpB*T2};
    //------------------------------------
#endif /*SMPC_VARIABLE_T_h*/

    
    for (int i=0; i<N; i++)
    {
#ifdef SMPC_VARIABLE_T_h
        //------------------------------------
        double Ti = T[i];
        double T2 = Ti*Ti/2;
        double hi = h[i];


        /* Control matrix. */
        double B[3] = {T2*Ti/3 - hi*Ti, T2, Ti};


        /* inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if T^3/6-h*T is
         * not equal to zero). The two elements one the main diagonal are equal, and only one of them 
         * is stored, which is equal to
            1/(T^3/6 - h*T)
         */
        double iCpB = 1/(B[0]);


        /* inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
            iCpB_CpA = [a b c 0 0 0;
                        0 0 0 a b c];

            a = iCpB
            b = iCpB*T
            c = iCpB*T^2/2
         * Only a,b and c are stored.
         */
        double iCpB_CpA[3] = {iCpB, iCpB*Ti, iCpB*T2};
        //------------------------------------
#endif /*SMPC_VARIABLE_T_h*/


        control[0] = -iCpB_CpA[0]*prev_state[0] - iCpB_CpA[1]*prev_state[1] - iCpB_CpA[2]*prev_state[2] + iCpB*x_coord[i];
        control[1] = -iCpB_CpA[0]*prev_state[3] - iCpB_CpA[1]*prev_state[4] - iCpB_CpA[2]*prev_state[5] + iCpB*y_coord[i];

        cur_state[0] = prev_state[0] + Ti*prev_state[1] + T2*prev_state[2] + B[0]*control[0];
        cur_state[1] =                    prev_state[1] + Ti*prev_state[2] + B[1]*control[0];
        cur_state[2] =                                       prev_state[2] + B[2]*control[0];
        cur_state[3] = prev_state[3] + Ti*prev_state[4] + T2*prev_state[5] + B[0]*control[1];
        cur_state[4] =                    prev_state[4] + Ti*prev_state[5] + B[1]*control[1];
        cur_state[5] =                                       prev_state[5] + B[2]*control[1];


        prev_state = &X[NUM_STATE_VAR*i];
        cur_state = &X[NUM_STATE_VAR*(i+1)];
        control = &control[NUM_CONTROL_VAR];
    }


    // go back to bar states
    cur_state = X;
    for (int i=0; i<N; i++)
    {
        state_handling::tilde_to_bar (angle_sin[i], angle_cos[i], cur_state);
        cur_state = &cur_state[NUM_STATE_VAR];
    }
}
