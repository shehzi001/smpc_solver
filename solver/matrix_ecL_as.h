/**
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:12 MSD
 */


#ifndef MATRIX_ECL_AS_H
#define MATRIX_ECL_AS_H


/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "problem_param.h"


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


/**
 * @brief Initializes lower diagonal matrix @ref pCholesky "L" and 
 * performs backward and forward substitutions using this matrix.
 */
class matrix_ecL_as
{
    public:
        /*********** Constructors / Destructors ************/
        matrix_ecL_as(const int);
        ~matrix_ecL_as();

        void form (const problem_parameters&);

        void solve_backward (const int, double *);
        void solve_forward (const int, double *);

        double *ecL;

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
/// @}

#endif /*MATRIX_ECL_AS_H*/
