/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef MATRIX_E_H
#define MATRIX_E_H
/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "problem_param.h"


/****************************************
 * DEFINES
 ****************************************/

using namespace std;

/// @addtogroup gINTERNALS
/// @{

/**
 * @brief Implements multiplication of matrix E and E' by a vector.
 */
class matrix_E
{
    public:
        /*********** Constructors / Destructors ************/
        matrix_E (){};
        ~matrix_E (){};

        void form_Ex (const problem_parameters&, const double *, double *);
        void form_ETx (const problem_parameters&, const double *, double *);
};
/// @}
#endif /*MATRIX_E_H*/
