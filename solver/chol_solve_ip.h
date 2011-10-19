/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef CHOL_SOLVE_IP_H
#define CHOL_SOLVE_IP_H
/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "matrix_E.h"
#include "matrix_ecL_ip.h"
#include "problem_param.h"


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
class chol_solve_ip
{
    public:
        /*********** Constructors / Destructors ************/
        chol_solve_ip (const int);
        ~chol_solve_ip();

        void solve(const problem_parameters*, const double *, const double *, const double *, double *);

    private:
        /// matrix of equality constraints
        matrix_E E;

        /// L for equality constraints, see '@ref pDetCholesky'
        matrix_ecL_ip ecL;

        /// Lagrange multipliers
        double *w;
};
/// @}
#endif /*CHOL_SOLVE_IP_H*/
