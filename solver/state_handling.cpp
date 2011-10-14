/** 
 * @file
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
     * @brief Converts state from @ref pX_tilde "X_tilde" to @ref pX_bar "X_bar".
     *
     * @param[in] sinA sin of the rotation angle.
     * @param[in] cosA cos of the rotation angle.
     * @param[in,out] state the state (@ref pX_tilde "X_tilde").
     */
    void tilde_to_bar (const double sinA, const double cosA, double *state)
    {
        double tmp =  cosA*state[0] + sinA*state[3];
        state[3]   = -sinA*state[0] + cosA*state[3];
        state[0]   = tmp;
    }

    /**
     * @brief Converts state from @ref pX_tilde "X_tilde" to @ref pX_bar "X_bar".
     *
     * @param[in] sinA sin of the rotation angle.
     * @param[in] cosA cos of the rotation angle.
     * @param[in,out] state the state (@ref pX_bar "X_bar").
     */
    void bar_to_tilde (const double sinA, const double cosA, double *state)
    {
        double tmp = cosA*state[0] - sinA*state[3];
        state[3]   = sinA*state[0] + cosA*state[3];
        state[0]   = tmp;
    }


    /**
     * @brief Converts state from @ref pX_tilde "X_tilde" to original variables.
     *
     * @param[in] h @ref ph "hCoM/gravity".
     * @param[in,out] state the state (@ref pX_tilde "X_tilde").
     */
    void tilde_to_orig (const double h, double *state)
    {
        state[0] = state[0] + h * state[2];
        state[3] = state[3] + h * state[5];
    }


    /**
     * @brief Returns the next state as @ref pX_tilde "X_tilde".
     *
     * @param[in] sp parameters. 
     * @param[in] X a solution.
     * @param[in,out] state the state (#NUM_STATE_VAR elements).
     */
    void get_next_state_tilde (const problem_parameters* sp, const double *X, double *state)
    {
        for (int i = 0; i < NUM_STATE_VAR; i++)
        {
            state[i] = X[i];
        }

        state_handling::bar_to_tilde (sp->angle_sin[0], sp->angle_cos[0], state);
    }


    /**
     * @brief Returns the next state represented by original variables.
     *
     * @param[in] sp parameters. 
     * @param[in] X a solution.
     * @param[in,out] state the state (#NUM_STATE_VAR elements).
     */
    void get_next_state (const problem_parameters* sp, const double *X, double *state)
    {
        get_next_state_tilde (sp, X, state);
        state_handling::tilde_to_orig (sp->h[0], state);
    }
}
