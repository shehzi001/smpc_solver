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
    void tilde_to_bar (const double, const double, double *);
    void bar_to_tilde (const double, const double, double *);
    void tilde_to_orig (const double, double *);

    void get_next_state_tilde (const solver_parameters*, const double *, double *);
    void get_next_state (const solver_parameters*, const double *, double *);
}
/// @}
#endif /*STATE_HANDLING_H*/
