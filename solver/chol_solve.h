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

        void solve(const solver_parameters&, const double *, const double *, double *);

        void up_resolve(const solver_parameters&, const double *, const int, const int *, const double *, double *);

#ifdef QPAS_DOWNDATE
        double * get_lambda();
        void down_resolve(const solver_parameters&, const double *, const int, const int *, const int, const double *, double *);
#endif


    private:
        void update (const solver_parameters&, const int, const int *);
        void update_z (const solver_parameters&, const double *, const int, const int *, const double *);
#ifdef QPAS_DOWNDATE
        void downdate(const solver_parameters&, const int, const int, const double *);
        void downdate_z (const solver_parameters&, const int, const int *, const int, const double *);
#endif

        void resolve (const solver_parameters&, const double *, const int, const int *, const double *, double *);

        void form_Ex (const solver_parameters&, const double *, double *);
        void form_ETx (const solver_parameters&, const double *, double *);

        void solve_forward(double *);
        void solve_backward(double *);

        void form_sa_row(const solver_parameters&, const int, const int, double *);


// ----------------------------------------------
// variables

        /// L for equality constraints, see '@ref pDetCholesky'
        double *ecL;

        /// L for inequality constraints, see '@ref pDetCholesky'
        double **icL;   

        /// All lines of #icL are stored in one chunk of memory.
        double *icL_mem;   

        /// Vector of Lagrange multipliers
        double *nu;

        /// - (X + solver_parameters#iHg)
        double *XiHg;

        /// Vector @ref pz "z".
        double *z;

        /// number of states in the preview window
        int N;

        /// An instance of #L_initializer class
        L_initializer L_init;
};
/// @}
#endif /*CHOL_SOLVE_H*/
