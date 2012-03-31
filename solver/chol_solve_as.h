/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 */


#ifndef CHOL_SOLVE_AS_H
#define CHOL_SOLVE_AS_H
/****************************************
 * INCLUDES 
 ****************************************/

#include <vector>

#include "smpc_common.h"
#include "matrix_E.h"
#include "matrix_ecL_as.h"
#include "problem_param.h"
#include "constraint.h"


/****************************************
 * DEFINES
 ****************************************/

using namespace std;

/// @addtogroup gAS
/// @{

/**
 * @brief Solves @ref pKKT "KKT system" using 
 * @ref pCholesky "Cholesky decomposition".
 */
class chol_solve_as
{
    public:
        /*********** Constructors / Destructors ************/
        chol_solve_as (const int);
        ~chol_solve_as();

        void solve(const problem_parameters&, const double *, const double *, double *);

        void up_resolve(const problem_parameters&, const double *, const std::vector<constraint>&, const int, const int *, const double *, double *);

        double * get_lambda(const problem_parameters&);
        void down_resolve(const problem_parameters&, const double *, const std::vector<constraint>&, const int, const int *, const int, const double *, double *);


    private:
        void update (const problem_parameters&, const constraint&, const int, const int);
        void update_z (const problem_parameters&, const double *, const constraint&, const int, const double *);
        void downdate(const problem_parameters&, const int, const int, const double *);

        void resolve (const problem_parameters&, const double *, const std::vector<constraint>&, const int, const int *, const double *, double *);

        void form_sa_row(const problem_parameters&, const constraint&, const int, double *);


// ----------------------------------------------
// variables
        /// Vector of Lagrange multipliers
        double *nu;

        /// - (X + problem_parameters#iHg)
        double *XiHg;

        /// matrix of equality constraints
        matrix_E E;

        /// L for equality constraints, see '@ref pDetCholesky'
        matrix_ecL_as ecL;

        /// L for inequality constraints, see '@ref pDetCholesky'
        double **icL;   

        /// All lines of #icL are stored in one chunk of memory.
        double *icL_mem;   

        /// Vector @ref pz "z".
        double *z;
};
/// @}
#endif /*CHOL_SOLVE_AS_H*/
