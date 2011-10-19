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
#include "problem_param.h"


/****************************************
 * DEFINES
 ****************************************/

/// The number of elements in 6x6 matrix.
#define MATRIX_SIZE_6x6 36
#define MATRIX_SIDE_6x6 6


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
class matrix_ecL_ip
{
    public:
        /*********** Constructors / Destructors ************/
        matrix_ecL_ip(const int);
        ~matrix_ecL_ip();

        void form (const problem_parameters*, const int, const double *);

        void solve_backward (const problem_parameters*, double *);
        void solve_forward (const problem_parameters*, double *);

        double *ecL;



    private:
        void chol_dec (double *);
        void form_M (const double, const double, const double*, const double*);
        void form_MBiPB (const double *, const double, double *);
        void form_MAT (const double, const double);
        void form_AMATMBiPB(const double, const double, double *);

        void form_L_non_diag(const double *, double *);
        void form_L_diag(double *);
        void form_L_diag(const double *, double *);


        // intermediate results used in computation of L
        double *M;         /// R * inv(Q) * R'
        double *MAT;       /// M * A'
};
/// @}

#endif /*MATRIX_ECL_IP_H*/
