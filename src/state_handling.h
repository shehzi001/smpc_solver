/**
 * @file
 * @author Alexander Sherikov
 * @date 14.09.2011 17:03:42 MSD
 */


#ifndef STATE_HANDLING_H
#define STATE_HANDLING_H

/****************************************
 * INCLUDES
 ****************************************/
#include "smpc_common.h"


/****************************************
 * PROTOTYPES 
 ****************************************/
/// @addtogroup gINTERNALS
/// @{

/**
 * @brief Various operations on the vector of states.
 */
namespace state_handling
{
    void tilde_to_bar (double *, double, double);
    void bar_to_tilde (double *, double, double);
    void tilde_to_orig (double *, double);

    void get_next_state_tilde (double *, double *, chol_solve_param);
    void get_next_state (double *, double *, chol_solve_param);
}
/// @}
#endif /*STATE_HANDLING_H*/
