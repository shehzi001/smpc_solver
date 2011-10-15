/**
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:12 MSD
 */


#ifndef MATRIX_ECL_H
#define MATRIX_ECL_H


/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "problem_param.h"


/****************************************
 * DEFINES
 ****************************************/


/****************************************
 * TYPEDEFS 
 ****************************************/

using namespace std;

/// @addtogroup gINTERNALS
/// @{

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

        void form (const problem_parameters*, const int);

        void solve_backward (const problem_parameters*, double *);
        void solve_forward (const problem_parameters*, double *);

        double *ecL;

    private:
        void chol_dec (double *);

        void form_iQBiPB (const double *, const double *, const double);
        void form_iQAT (const double, const double, const double *);
        void form_AiQATiQBiPB (const double, const double);

        void form_L_non_diag(const double *, double *);
        void form_L_diag(const double *, double *);


        // intermediate results used in computation of L
        double *iQBiPB;     /// inv(Q) + B * inv(P) * B'
        double *iQAT;       /// inv(Q) * A'
        double *AiQATiQBiPB;/// A * inv(Q) * A' + inv(Q) + B * inv(P) * B'
};
/// @}

#endif /*MATRIX_ECL_H*/
