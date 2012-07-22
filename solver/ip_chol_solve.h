/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef IP_CHOL_SOLVE_H
#define IP_CHOL_SOLVE_H
/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "ip_matrix_E.h"
#include "ip_matrix_ecL.h"
#include "ip_problem_param.h"


/****************************************
 * DEFINES
 ****************************************/

using namespace std;

/// @addtogroup gIP
/// @{

namespace IP
{
    /**
     * @brief Solves @ref pKKT "KKT system" using 
     * @ref pCholesky "Cholesky decomposition".
     */
    class chol_solve
    {
        public:
            /*********** Constructors / Destructors ************/
            chol_solve (const int);
            ~chol_solve();

            void solve(const problem_parameters&, const double *, const double *, const double *, double *);

        private:
            /// matrix of equality constraints
            matrix_E E;

            /// L for equality constraints
            matrix_ecL ecL;

            /// Lagrange multipliers
            double *w;
    };
}
/// @}
#endif /*IP_CHOL_SOLVE_H*/
