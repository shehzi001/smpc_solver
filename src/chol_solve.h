/**
 * @file
 * @brief 
 *
 * @author Alexander Sherikov
 * @date 19.07.2011 15:56:18 MSD
 * @todo add description
 */


#ifndef CHOL_SOLVE_H
#define CHOL_SOLVE_H
/****************************************
 * INCLUDES 
 ****************************************/

#include "smpc_common.h"


/****************************************
 * DEFINES
 ****************************************/
#define MATRIX_DIM_SIZE 3
#define MATRIX_SIZE 9

using namespace std;

/**
 * @brief 
 *
 *
 * It is implicitly supposed that we have 6 state variables and 2
 * control variables. Defines are used for convenience only.
 * The structure of matrices is also hardcoded. (refer to the paper
 * for more information).
 *
 * All matrices used in computations are of size 3x3, and are stored
 * as vectors. The elements are numerated in the following way:
 *  0   3   6
 *  1   4   7
 *  2   5   8
 *
 */
class chol_solve
{
    public:
        /*********** Constructors / Destructors ************/
        chol_solve (int);
        ~chol_solve();

        void solve(chol_solve_param, double *, double *);
        void update (chol_solve_param, int, int *);
        int downdate (chol_solve_param, int, int*, double *);
        void resolve (chol_solve_param, int, int *, double *, double *, bool after_update = true);


#ifndef QPAS_DEBUG
    private:
#endif
        void form_Ex (chol_solve_param, double *, double *);
        void form_ETx (chol_solve_param, double *, double *);

        void chol_dec (double *);

        void form_iQBiPB (double *, double *, double);
        void form_iQAT (double, double, double *);
        void form_AiQATiQBiPB (double, double);

        void form_L_non_diag(double *, double *);
        void form_L_diag(double *, double *);
        void form_L(chol_solve_param);

        void solve_forward(double *);
        void solve_backward(double *);

        void form_a_row(chol_solve_param, int, int, double *);


// ----------------------------------------------
// variables

        /** L for equality constraints
           full L                           |   L (ecL) in this class
            a                               |
            0   a                           |           a
            b   0   c                   0   |           b   c
            0   b   0   c                   |   ===>        d   e
                    d   0   e               |   ===>            f   g
                    0   d   0   e           |                       ...
                0           f   0   g       |
                            0   f   0   g   |
                                    ....    |

         Here a,c,e,g - lower diagonal 3x3 matrices; b,d,f - upper triangular
         3x3 matrices.
        ------------------------------------------------------------------
          ecL is stored as a sequence of vectors of 9 elements:
         
        0  3  6 
        1  4  7
        2  5  8

        9  12 15  18 21 24
        10 13 16  19 22 25
        11 14 17  20 23 26

                    ...         ...
        */
        double *ecL;

        /** L for inequality constraints
            (N*2)x(N*NUM_VAR) matrix */
        double **icL;   

        /** All lines of #icL are stored in one chunk of memory. */
        double *icL_mem;   

        // intermediate results used in computation of L
        double *iQBiPB;     /// inv(Q) + B * inv(P) * B'
        double *iQAT;       /// inv(Q) * A'
        double *AiQATiQBiPB;/// A * inv(Q) * A' + inv(Q) + B * inv(P) * B'

        // these vectors are used in ReSolve and Solve
        // respectively, the memory is allocated once in
        // constructor
        double *nu;
        double *ViHg;

        /** Vector #z = -[E;Aw]*inv(H)*(H*x+g) \ L. */
        double *z;

        /// number of states in the preview window
        int N;
};
#endif /*CHOL_SOLVE_H*/
