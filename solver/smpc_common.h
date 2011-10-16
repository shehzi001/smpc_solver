/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:37:28 MSD
 */


#ifndef SMPC_COMMON_H
#define SMPC_COMMON_H

/****************************************
 * INCLUDES 
 ****************************************/

#include <cstddef>
#include "common_const.h"

/****************************************
 * DEFINES
 ****************************************/

/// @addtogroup gINTERNALS
/// @{

/// Allow variable problem_parameters#T and problem_parameters#h
#define SMPC_VARIABLE_T_h
/// Allow removal of constraints from active set
#define QPAS_DOWNDATE


/****************************************
 * TYPEDEFS 
 ****************************************/

/****************************************
 * PROTOTYPES 
 ****************************************/
///@}
#endif /*SMPC_COMMON_H*/

