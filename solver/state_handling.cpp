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
     * @brief Converts state from original variables to @ref pX_tilde "X_tilde".
     *
     * @param[in] h @ref ph "hCoM/gravity".
     * @param[in,out] state the state (position, velocity and acceleration of CoM).
     */
    void orig_to_tilde (const double h, double *state)
    {
        state[0] = state[0] - h * state[2];
        state[3] = state[3] - h * state[5];
    }


    /**
     * @brief Returns the next state as @ref pX_tilde "X_tilde".
     *
     * @param[in] sp parameters. 
     * @param[in] X a solution.
     * @param[in] ind index of the state.
     * @param[in,out] state the state (#SMPC_NUM_STATE_VAR elements).
     */
    void get_state_tilde (const problem_parameters& sp, const double *X, const int ind, double *state)
    {
        int index;
        if (ind >= sp.N)
        {
            index = sp.N - 1;
        }
        else
        {
            index = ind;
        }

        for (int i = 0; i < SMPC_NUM_STATE_VAR; i++)
        {
            state[i] = X[index*SMPC_NUM_STATE_VAR + i];
        }
    }


    /**
     * @brief Returns the next state represented by original variables.
     *
     * @param[in] sp parameters. 
     * @param[in] X a solution.
     * @param[in] ind index of the state.
     * @param[in,out] state the state (#SMPC_NUM_STATE_VAR elements).
     */
    void get_state (const problem_parameters& sp, const double *X, const int ind, double *state)
    {
        int index;
        if (ind >= sp.N)
        {
            index = sp.N - 1;
        }
        else
        {
            index = ind;
        }
        get_state_tilde (sp, X, ind, state);
        state_handling::tilde_to_orig (sp.spar[index].h, state);
    }


    /**
     * @brief Returns the controls,that must be applied to reach the next state.
     *
     * @param[in] preview_window_size size of the preview window
     * @param[in] X a solution
     * @param[in] ind index of the state.
     * @param[in,out] controls the controls (2 double values)
     */
    void get_controls (const int preview_window_size, const double *X, const int ind, double *controls)
    {
        int index;
        if (ind >= preview_window_size)
        {
            index = preview_window_size-1;
        }
        else
        {
            index = ind;
        }
        controls[0] = X[preview_window_size*SMPC_NUM_STATE_VAR + index*SMPC_NUM_CONTROL_VAR + 0];
        controls[1] = X[preview_window_size*SMPC_NUM_STATE_VAR + index*SMPC_NUM_CONTROL_VAR + 1];
    }
}
