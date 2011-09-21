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
 * @brief Solves @ref pPD_KKT "KKT system" using 
 * @ref pCholesky "Cholesky decomposition".
 */
class chol_solve
{
    public:
        /*********** Constructors / Destructors ************/
        chol_solve (int);
        ~chol_solve();

        void solve(chol_solve_param, double *, double *);

        void up_resolve(chol_solve_param, int, int *, double *, double *);

#ifdef QPAS_DOWNDATE
        double * get_lambda();
        void down_resolve(chol_solve_param, int, int *, int, double *, double *);
#endif


    private:
        void update (chol_solve_param, int, int *);
        void update_z (chol_solve_param, int, int *, double *);
#ifdef QPAS_DOWNDATE
        void downdate(chol_solve_param, int, int, double *);
        void downdate_z (chol_solve_param, int, int *, int, double *);
#endif

        void resolve (chol_solve_param, int, int *, double *, double *);

        void form_Ex (chol_solve_param, double *, double *);
        void form_ETx (chol_solve_param, double *, double *);

        void solve_forward(double *);
        void solve_backward(double *);

        void form_sa_row(chol_solve_param, int, int, double *);


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

        /// - (X + chol_solve_param#iHg)
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
