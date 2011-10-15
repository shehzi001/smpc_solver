/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef CHOL_SOLVE_H
#define CHOL_SOLVE_H
/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "L_initializer.h"


/****************************************
 * DEFINES
 ****************************************/

using namespace std;

/// @addtogroup gINTERNALS
/// @{

/**
 * @brief Solves @ref pKKT "KKT system" using 
 * @ref pCholesky "Cholesky decomposition".
 */
class chol_solve : public L_initializer
{
    public:
        /*********** Constructors / Destructors ************/
        chol_solve (int);
        ~chol_solve();

        void form_Ex (const problem_parameters*, const double *, double *);
        void form_ETx (const problem_parameters*, const double *, double *);

        void solve_forward(const problem_parameters*, double *);
        void solve_backward(const problem_parameters*, double *);

// ----------------------------------------------
// variables

        /// L for equality constraints, see '@ref pDetCholesky'
        double *ecL;

        /// Vector of Lagrange multipliers
        double *nu;

        /// - (X + problem_parameters#iHg)
        double *XiHg;
};
/// @}
#endif /*CHOL_SOLVE_H*/
