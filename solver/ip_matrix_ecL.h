/**
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:12 MSD
 */


#ifndef MATRIX_ECL_IP_H
#define MATRIX_ECL_IP_H


/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"
#include "ip_problem_param.h"

using namespace std;

/// @addtogroup gIP
/// @{

namespace IP
{

/****************************************
 * DEFINES
 ****************************************/

/// The number of elements in 6x6 matrix.
#define MATRIX_SIZE_6x6 36


/****************************************
 * TYPEDEFS 
 ****************************************/
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

            void form (const problem_parameters&, const double *);

            void solve_backward (const int, double *);
            void solve_forward (const int, double *);

            double *ecL;



        private:
            void chol_dec (double *);
            void form_M (const double, const double, const double*, const double*);
            void form_MAT (const double, const double);
            void form_AMATMBiPB(const double, const double, const double *, const double, double *);

            void form_L_non_diag(const double *, double *);
            void form_L_diag(const double *, const double, double *);
            void form_L_diag(const double *, double *);


            // intermediate results used in computation of L
            double M[MATRIX_SIZE_6x6];         /// R * inv(Q) * R'
            double MAT[MATRIX_SIZE_6x6];       /// M * A'
    };
}
/// @}

#endif /*MATRIX_ECL_IP_H*/
