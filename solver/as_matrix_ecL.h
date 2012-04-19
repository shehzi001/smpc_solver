/**
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:12 MSD
 */


#ifndef AS_MATRIX_ECL_H
#define AS_MATRIX_ECL_H


/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "as_problem_param.h"


/// @addtogroup gAS
/// @{

/****************************************
 * DEFINES
 ****************************************/

/// The number of elements in 3x3 matrix.
#define MATRIX_SIZE_3x3 9

/****************************************
 * TYPEDEFS 
 ****************************************/

using namespace std;

namespace AS
{
    /**
     * @brief Initializes lower diagonal matrix @ref pCholesky "L" and 
     * performs backward and forward substitutions using this matrix.
     */
    class matrix_ecL
    {
        public:
            /*********** Constructors / Destructors ************/
            matrix_ecL(const int);
            ~matrix_ecL();

            void form (const problem_parameters&);

            void solve_backward (const int, double *) const;
            void solve_forward (const int, double *, const int start_ind = 0) const;

            double *ecL;
            double **ecL_diag;
            double **ecL_ndiag;

        private:
            void chol_dec (double *);

            void form_iQBiPB (const double *, const double *, const double, double*);
            void form_iQAT (const double, const double, const double *);
            void form_AiQATiQBiPB (const problem_parameters&, const state_parameters&, double *);

            void form_L_non_diag(const double *, double *);
            void form_L_diag(const double *, double *);


            // intermediate results used in computation of L
            double *iQAT;       /// inv(Q) * A'
    };
}
/// @}

#endif /*AS_MATRIX_ECL_H*/
