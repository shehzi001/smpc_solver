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
#include "L_initializer.h"


/****************************************
 * DEFINES
 ****************************************/

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

        // these vectors are used in ReSolve and Solve
        // respectively, the memory is allocated once in
        // constructor
        double *nu;
        double *ViHg;

        /** Vector #z = -[E;Aw]*inv(H)*(H*x+g) \ L. */
        double *z;

        /// number of states in the preview window
        int N;
        L_initializer L_init;
};
#endif /*CHOL_SOLVE_H*/
