/** 
 * @file
 * @author Alexander Sherikov
 */


#include "WMG.h"


/**
 * @brief Initialize state (#A) and control (#B) matrices for inverted 
 * pendulum model.
 *
 * @param[in] sampling_time period of time T.
 */
void WMG::initABMatrices (const double sampling_time)
{
    A = new double[9];
    B = new double[3];

    A[0] = A[4] = A[8] = 1;
    A[1] = A[2] = A[5] = 0;
    A[3] = A[7] = sampling_time;
    A[6] = sampling_time * sampling_time/2;

    B[0] = sampling_time * sampling_time * sampling_time / 6;
    B[1] = sampling_time * sampling_time/2;
    B[2] = sampling_time;
}



/**
 * @brief Calculate next state using inverted pendulum model (#A and #B matrices).
 *
 * @param[in] control 1x2 vector of controls
 * @param[in,out] state 1x6 state vector
 *
 * @attention If #A or #B are not initialized, the function does nothing.
 */
void WMG::calculateNextState (smpc::control &control, smpc::state_orig &state)
{
    if ((A == NULL) || (B == NULL))
    {
        return;
    }


    state.x()  = state.x()  * A[0]
               + state.vx() * A[3]
               + state.ax() * A[6]
               + control.jx() * B[0];

    state.vx() = state.vx() * A[4]
               + state.ax() * A[7]
               + control.jx() * B[1];

    state.ax() = state.ax() * A[8]
               + control.jx() * B[2];


    state.y()  = state.y()  * A[0]
               + state.vy() * A[3]
               + state.ay() * A[6]
               + control.jy() * B[0];

    state.vy() = state.vy() * A[4]
               + state.ay() * A[7]
               + control.jy() * B[1];

    state.ay() = state.ay() * A[8]
               + control.jy() * B[2];
}

