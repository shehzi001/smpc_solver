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
    void tilde_to_orig (const double, double *);
    void orig_to_tilde (const double, double *);

    void get_controls (const int, const double *, const int ind, double *);
}
/// @}
#endif /*STATE_HANDLING_H*/
