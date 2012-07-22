/**
 * @file
 * @author Alexander Sherikov
 * @date 26.04.2012 19:36:23 MSD
 */


#ifndef QP_H
#define QP_H

/****************************************
 * TEMPLATES
 ****************************************/
#include "smpc_common.h"
#include "state_handling.h"


/****************************************
 * TEMPLATES
 ****************************************/

/**
 * @brief Generates an initial feasible point. 
 *
 * @param[in] ppar problem parameters
 * @param[in] x_coord x coordinates of points satisfying constraints
 * @param[in] y_coord y coordinates of points satisfying constraints
 * @param[in] init_state current state
 * @param[in] tilde_state if true the state is assumed to be in @ref pX_tilde "X_tilde" form
 * @param[in,out] X initial guess / solution of optimization problem
 */
template <class PP>
void form_init_fp_tilde (
        const PP &ppar,
        const double *x_coord, 
        const double *y_coord, 
        const double *init_state,
        const bool tilde_state,
        double* X)
{
    double *control = &X[SMPC_NUM_STATE_VAR*ppar.N];
    double *cur_state = X;
    double X_tilde[6] = {
        init_state[0], init_state[1], init_state[2],
        init_state[3], init_state[4], init_state[5]};
    if (!tilde_state)
    {
        state_handling::orig_to_tilde (ppar.h_initial, X_tilde);
    }
    const double *prev_state = X_tilde;

    
    for (int i = 0; i < ppar.N; ++i)
    {
        //------------------------------------
        /* inv(Cp*B). This is a [2 x 2] diagonal matrix (which is invertible if T^3/6-h*T is
         * not equal to zero). The two elements on the main diagonal are equal, and only one of them 
         * is stored, which is equal to
            1/(T^3/6 - h*T)
         */
        double iCpB = 1/(ppar.spar[i].B[0]);

        /* inv(Cp*B)*Cp*A. This is a [2 x 6] matrix with the following structure
            iCpB_CpA = [a b c 0 0 0;
                        0 0 0 a b c];

            a = iCpB
            b = iCpB*T
            c = iCpB*T^2/2
         * Only a,b and c are stored.
         */
        double iCpB_CpA[3] = {iCpB, iCpB*ppar.spar[i].A3, iCpB*ppar.spar[i].A6};
        //------------------------------------


        control[0] = -iCpB_CpA[0]*prev_state[0] - iCpB_CpA[1]*prev_state[1] - iCpB_CpA[2]*prev_state[2] + iCpB*x_coord[i];
        control[1] = -iCpB_CpA[0]*prev_state[3] - iCpB_CpA[1]*prev_state[4] - iCpB_CpA[2]*prev_state[5] + iCpB*y_coord[i];

        cur_state[0] = prev_state[0] + ppar.spar[i].A3*prev_state[1] + ppar.spar[i].A6*prev_state[2] + ppar.spar[i].B[0]*control[0];
        cur_state[1] =                                 prev_state[1] + ppar.spar[i].A3*prev_state[2] + ppar.spar[i].B[1]*control[0];
        cur_state[2] =                                                                 prev_state[2] + ppar.spar[i].B[2]*control[0];
        cur_state[3] = prev_state[3] + ppar.spar[i].A3*prev_state[4] + ppar.spar[i].A6*prev_state[5] + ppar.spar[i].B[0]*control[1];
        cur_state[4] =                                 prev_state[4] + ppar.spar[i].A3*prev_state[5] + ppar.spar[i].B[1]*control[1];
        cur_state[5] =                                                                 prev_state[5] + ppar.spar[i].B[2]*control[1];


        prev_state = &X[SMPC_NUM_STATE_VAR*i];
        cur_state = &X[SMPC_NUM_STATE_VAR*(i+1)];
        control = &control[SMPC_NUM_CONTROL_VAR];
    }
}

#endif /*QP_H*/

