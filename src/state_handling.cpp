/** 
 * @file
 * @brief 
 *
 * @author Alexander Sherikov
 * @date 14.09.2011 16:55:42 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "state_handling.h"



/****************************************
 * FUNCTIONS 
 ****************************************/

namespace state_handling
{
    /**
     * @brief Converts state from #X_tilde to #X_bar.
     *
     * @param[in,out] state the state (#X_tilde).
     * @param[in] sinA sin of the rotation angle.
     * @param[in] cosA cos of the rotation angle.
     */
    void tilde_to_bar (double *state, double sinA, double cosA)
    {
        double tmp =  cosA*state[0] + sinA*state[3];
        state[3]   = -sinA*state[0] + cosA*state[3];
        state[0]   = tmp;
    }

    /**
     * @brief Converts state from #X_tilde to #X_bar.
     *
     * @param[in,out] state the state (#X_bar).
     * @param[in] sinA sin of the rotation angle.
     * @param[in] cosA cos of the rotation angle.
     */
    void bar_to_tilde (double *state, double sinA, double cosA)
    {
        double tmp = cosA*state[0] - sinA*state[3];
        state[3]   = sinA*state[0] + cosA*state[3];
        state[0]   = tmp;
    }


    /**
     * @brief Converts state from #X_tilde to original variables.
     *
     * @param[in,out] state the state (#X_tilde).
     * @param[in] h #hCoM/#gravity.
     */
    void tilde_to_orig (double *state, double h)
    {
        state[0] = state[0] + h * state[2];
        state[0] = state[3] + h * state[5];
    }


    /**
     * @brief Returns the next state as #X_tilde.
     *
     * @param[in,out] state the state (#NUM_STATE_VAR elements).
     * @param[in] X a solution.
     * @param[in] csp parameters. 
     */
    void get_next_state_tilde (double *state, double *X, chol_solve_param csp)
    {
        for (int i = 0; i < NUM_STATE_VAR; i++)
        {
            state[i] = X[i];
        }

        state_handling::bar_to_tilde (state, csp.angle_sin[0], csp.angle_cos[0]);
    }


    /**
     * @brief Returns the next state represented by original variables.
     *
     * @param[in,out] state the state (#NUM_STATE_VAR elements).
     * @param[in] X a solution.
     * @param[in] csp parameters. 
     */
    void get_next_state (double *state, double *X, chol_solve_param csp)
    {
        get_next_state_tilde (state, X, csp);
        state_handling::tilde_to_orig (state, csp.h[0]);
    }
}
