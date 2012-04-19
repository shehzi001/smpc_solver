/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef AS_MATRIX_E_H
#define AS_MATRIX_E_H
/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "as_problem_param.h"


/****************************************
 * DEFINES
 ****************************************/

using namespace std;

/// @addtogroup gAS
/// @{
namespace AS
{
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
            void form_i2HETx (const problem_parameters&, const double *, double *);
    };
}
/// @}
#endif /*AS_MATRIX_E_H*/
