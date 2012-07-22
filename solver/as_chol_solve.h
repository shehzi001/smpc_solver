/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef AS_CHOL_SOLVE_H
#define AS_CHOL_SOLVE_H
/****************************************
 * INCLUDES 
 ****************************************/

#include <vector>

#include "smpc_common.h"
#include "as_matrix_E.h"
#include "as_matrix_ecL.h"
#include "as_problem_param.h"
#include "as_constraint.h"


/****************************************
 * DEFINES
 ****************************************/

using namespace std;

/// @addtogroup gAS
/// @{
namespace AS
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

            void solve(const AS::problem_parameters&, const double *, double *);

            void up_resolve(const AS::problem_parameters&, const vector<AS::constraint>&, const double *, double *);

            double * get_lambda(const AS::problem_parameters&);
            void down_resolve(const AS::problem_parameters&, const vector<AS::constraint>&, const int, const double *, double *);


        private:
            void update (const AS::problem_parameters&, const AS::constraint&, const int);
            void update_z (const AS::problem_parameters&, const AS::constraint&, const int, const double *);
            void downdate(const AS::problem_parameters&, const int, const int, const double *);

            void resolve (const AS::problem_parameters&, const vector<AS::constraint>&, const double *, double *);

            void form_sa_row(const AS::problem_parameters&, const AS::constraint&, const int, double *);


    // ----------------------------------------------
    // variables
            /// Vector of Lagrange multipliers
            double *nu;

            /// matrix of equality AS::constraints
            AS::matrix_E E;

            /// L for equality AS::constraints
            AS::matrix_ecL ecL;

            /// L for inequality AS::constraints
            double **icL;   

            /// All lines of #icL are stored in one chunk of memory.
            double *icL_mem;   

            /// Vector @ref pz "z".
            double *z;
    };
}
/// @}
#endif /*AS_CHOL_SOLVE_H*/
