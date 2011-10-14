/**
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:37:28 MSD
 */


#ifndef SMPC_COMMON_H
#define SMPC_COMMON_H

/****************************************
 * INCLUDES 
 ****************************************/

#include <cstddef>
#include "common_const.h"

/****************************************
 * DEFINES
 ****************************************/

/// @addtogroup gINTERNALS
/// @{

/// The size of one dimension of 3x3 matrix.
#define MATRIX_DIM_SIZE 3
/// The number of elements in 3x3 matrix.
#define MATRIX_SIZE 9

/// Allow variable solver_parameters#T and solver_parameters#h
#define SMPC_VARIABLE_T_h
/// Allow removal of constraints from active set
#define QPAS_DOWNDATE


/****************************************
 * TYPEDEFS 
 ****************************************/

/**
 * @brief A set of parameters used by #chol_solve class.
 */
class solver_parameters
{
    public:
        solver_parameters (
            const int N_,
            const double Alpha,
            const double Beta,
            const double Gamma,
            const double regularization)
        {
            N = N_;

            i2Q[0] = 1/(2*(Beta/2));
            i2Q[1] = 1/(2*(Alpha/2));
            i2Q[2] = 1/(2*regularization);

            i2P = 1/(2 * (Gamma/2));

            angle_cos = new double[N];
            angle_sin = new double[N];
            T = NULL;
            h = NULL;
#ifdef SMPC_VARIABLE_T_h
            dh = new double[N-1];
#endif
        };

        ~solver_parameters()
        {
            if (angle_cos != NULL)
                delete angle_cos;

            if (angle_sin != NULL)
                delete angle_sin;

#ifdef SMPC_VARIABLE_T_h
            if (dh != NULL)
                delete dh;
#endif
        };


        /** Number of iterations in a preview window. */
        int N;

    // static matrices and vectors
        ///@{
        /** State related penalty.*/
        double i2Q[3];
        ///@}

        ///@{
        /** Control related penalty. */
        double i2P;
        ///@}


        double *angle_cos;
        double *angle_sin;


    // parameters used in generation of A and B matrices
        /** Preview sampling time  */
        const double *T;
        /** h = @ref ph "hCoM/gravity". */
        const double *h;

#ifdef SMPC_VARIABLE_T_h
        double *dh;
#endif
};

/****************************************
 * PROTOTYPES 
 ****************************************/
///@}
#endif /*SMPC_COMMON_H*/

