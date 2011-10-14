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
class chol_solve
{
    public:
        /*********** Constructors / Destructors ************/
        chol_solve (int);
        ~chol_solve();

        void solve(const solver_parameters*, const double *, const double *, double *);


        void form_Ex (const solver_parameters*, const double *, double *);
        void form_ETx (const solver_parameters*, const double *, double *);

        void solve_forward(double *);
        void solve_backward(double *);

        void form_sa_row(const solver_parameters*, const int, const int, double *);


// ----------------------------------------------
// variables

        /// L for equality constraints, see '@ref pDetCholesky'
        double *ecL;

        /// Vector of Lagrange multipliers
        double *nu;

        /// - (X + solver_parameters#iHg)
        double *XiHg;

        /// number of states in the preview window
        int N;

        /// An instance of #L_initializer class
        L_initializer L_init;
};
/// @}
#endif /*CHOL_SOLVE_H*/
