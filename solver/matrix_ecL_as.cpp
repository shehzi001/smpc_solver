/** 
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:08 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/

#include "matrix_ecL_as.h"

#include <cmath> // sqrt


/****************************************
 * FUNCTIONS 
 ****************************************/

//==============================================
// constructors / destructors

matrix_ecL_as::matrix_ecL_as (const int N)
{
    ecL = new double[MATRIX_SIZE_3x3*N + MATRIX_SIZE_3x3*(N-1)]();

    iQAT = new double[MATRIX_SIZE_3x3];
#ifndef SMPC_VARIABLE_T_h
    AiQATiQBiPB = new double[MATRIX_SIZE_3x3];
#endif
}


matrix_ecL_as::~matrix_ecL_as()
{
    if (ecL != NULL)
        delete ecL;

#ifndef SMPC_VARIABLE_T_h
    if (AiQATiQBiPB != NULL)
        delete AiQATiQBiPB;
#endif

    if (iQAT != NULL)
        delete iQAT;
}
//==============================================



/**
 * @brief Performs Cholesky decomposition of 3x3 matrix.
 *
 * @param[in,out] mx9 a pointer to matrix, the result is 
 *                    stored in the same place.
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_as::chol_dec (double *mx9)
{
    // 1st line
    mx9[0] = sqrt (mx9[0]);

    // 2nd line
    mx9[1] /= mx9[0];
    mx9[4] = sqrt(mx9[4] - mx9[1]*mx9[1]); 

    // 3rd line
    mx9[2] /= mx9[0]; 
    mx9[5] = (mx9[5] - mx9[2]*mx9[1])/mx9[4];
    mx9[8] = sqrt(mx9[8] - mx9[5]*mx9[5] - mx9[2]*mx9[2]);

    // These elements must be 0. (but they are never used)
//    mx9[3] = mx9[6] = mx9[7] = 0;
}



/**
 * @brief Forms matrix iQBiPB = 0.5 * inv(Q) + 0.5 * B * inv(P) * B.
 *
 * @param[in] B a vector of 3 elements.
 * @param[in] i2Q a vector of 3 elements, which contains
 *              diagonal elements of 0.5 * inv(Q).
 * @param[in] i2P 0.5 * inv(P) (only one number)
 * @param[out] result the result
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_as::form_iQBiPB (const double *B, const double *i2Q, const double i2P, double* result)
{
    // diagonal elements
    result[0] = i2P * B[0]*B[0] + i2Q[0];
    result[4] = i2P * B[1]*B[1] + i2Q[1];
    result[8] = i2P * B[2]*B[2] + i2Q[2];

    // symmetric elements (no need to initialize all of them)
    result[1] = /*result[3] =*/ i2P * B[0]*B[1];
    result[2] = /*result[6] =*/ i2P * B[0]*B[2];
    result[5] = /*result[7] =*/ i2P * B[1]*B[2];
}



/**
 * @brief Forms matrix iQAT = 0.5 * inv (Q) * A'
 *
 * @param[in] A3 4th and 7th elements of A.
 * @param[in] A6 6th element of A.
 * @param[in] i2Q a vector of 3 elements, which contains
 *              diagonal elements of 0.5*inv(Q).
 */
void matrix_ecL_as::form_iQAT (const double A3, const double A6, const double *i2Q)
{
    iQAT[0] = i2Q[0];
    iQAT[1] = A3 * i2Q[1];
    iQAT[2] = A6 * i2Q[2];
    iQAT[4] = i2Q[1];
    iQAT[5] = A3 * i2Q[2];
    iQAT[8] = i2Q[2];
}


/**
 * @brief Forms matrix AiQATiQBiPB = 
 *  A * inv(Q) * A' + 0.5 * inv(Q) + 0.5 * B * inv(P) * B
 *
 * @param[in] ppar problem parameters.
 * @param[in] stp state parameters.
 * @param[out] result result
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_as::form_AiQATiQBiPB (const problem_parameters *ppar, const state_parameters stp, double *result)
{
    form_iQBiPB (stp.B, ppar->i2Q, ppar->i2P, result);

    // 1st column
    result[0] += iQAT[0] + stp.A3*iQAT[1] + stp.A6*iQAT[2];
    result[1] +=                  iQAT[1] + stp.A3*iQAT[2];
    result[2] +=                                   iQAT[2];

    // 2nd column
    // symmetric elements are not initialized
//    result[3] = iQBiPB[3] + A3*iQAT[4] + A6*iQAT[5];
    result[4] += iQAT[4] + stp.A3*iQAT[5];
    result[5] +=                  iQAT[5];

    // 3rd column
//    result[6] = iQBiPB[6] + A6*iQAT[8];
//    result[7] = iQBiPB[7] + A3*iQAT[8];
    result[8] += iQAT[8];
}



/**
 * @brief Forms a 3x3 matrix L(k+1, k), which lies below the 
 *  diagonal of L.
 *
 * @param[in] ecLp previous matrix lying on the diagonal of L
 * @param[in] ecLc the result is stored here
 */
void matrix_ecL_as::form_L_non_diag(const double *ecLp, double *ecLc)
{
    /* L(k+1,k) * L(k,k)' = - inv(Q) * A'
     *
     * xxx      xxx     xxx
     *  xx  *    xx =    xx
     *   x        x       x
     *
     * all matrices are upper triangular
     */

    // main diagonal
    ecLc[0] = -iQAT[0] / ecLp[0];
    ecLc[4] = -iQAT[4] / ecLp[4];
    ecLc[8] = -iQAT[8] / ecLp[8];

    // sub-diagonal 1
    ecLc[3] = (-iQAT[1] - ecLc[0]*ecLp[1]) / ecLp[4];
    ecLc[7] = (-iQAT[5] - ecLc[4]*ecLp[5]) / ecLp[8];

    // sub-diagonal 2
    ecLc[6] = (-iQAT[2] - ecLc[0]*ecLp[2] - ecLc[3]*ecLp[5]) / ecLp[8];
}



/**
 * @brief Forms a 3x3 matrix L(k+1, k+1), which lies on the main
 *  diagonal of L.
 *
 * @param[in] ecLp upper triangular matrix matrix lying to the left
 *                 from ecLc on the same level of L
 * @param[in] ecLc the result is stored here
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_as::form_L_diag(const double *ecLp, double *ecLc)
{
#ifndef SMPC_VARIABLE_T_h
    ecLc[0] = AiQATiQBiPB[0];
    ecLc[4] = AiQATiQBiPB[4];
    ecLc[8] = AiQATiQBiPB[8];
    // symmetric nondiagonal elements (no need to initialize all of them)
    ecLc[1] = AiQATiQBiPB[1];
    ecLc[2] = AiQATiQBiPB[2];
    ecLc[5] = AiQATiQBiPB[5];
#endif

    // L(k+1,k+1) = (- L(k+1,k) * L(k+1,k)') + (A * inv(Q) * A' + inv(Q) + B * inv(P) * B)
    // diagonal elements
    ecLc[0] -= ecLp[0]*ecLp[0] + ecLp[3]*ecLp[3] + ecLp[6]*ecLp[6];
    ecLc[4] -= ecLp[4]*ecLp[4] + ecLp[7]*ecLp[7];
    ecLc[8] -= ecLp[8]*ecLp[8];
    // symmetric nondiagonal elements (no need to initialize all of them)
    ecLc[1] -= /*ecLc[3] =*/ ecLp[3]*ecLp[4] + ecLp[6]*ecLp[7];
    ecLc[2] -= /*ecLc[6] =*/ ecLp[6]*ecLp[8];
    ecLc[5] -= /*ecLc[7] =*/ ecLp[7]*ecLp[8];

    // chol (L(k+1,k+1))
    chol_dec (ecLc);
}



/**
 * @brief Builds matrix L.
 *
 * @param[in] ppar parameters.
 * @param[in] N number of states in preview window.
 */
void matrix_ecL_as::form (const problem_parameters* ppar, const int N)
{
    int i;
    state_parameters stp;


    stp = ppar->spar[0];

    // form all matrices
#ifndef SMPC_VARIABLE_T_h
    form_iQAT (stp.A3, stp.A6, ppar->i2Q);
    form_AiQATiQBiPB (ppar, stp, AiQATiQBiPB);
#endif

    // the first matrix on diagonal
    form_iQBiPB (stp.B, ppar->i2Q, ppar->i2P, ecL);
    chol_dec (ecL);


    // offsets
    double *ecL_cur = &ecL[MATRIX_SIZE_3x3];
    double *ecL_prev = ecL;
    for (i = 1; i < N; i++)
    {
#ifdef SMPC_VARIABLE_T_h
        stp = ppar->spar[i];
        form_iQAT (stp.A3, stp.A6, ppar->i2Q);
#endif

        // form (b), (d), (f) ... 
        form_L_non_diag(ecL_prev, ecL_cur);
        // update offsets
        ecL_cur = &ecL_cur[MATRIX_SIZE_3x3];
        ecL_prev = &ecL_prev[MATRIX_SIZE_3x3];

        // form (c), (e), (g) ...
#ifdef SMPC_VARIABLE_T_h
        form_AiQATiQBiPB (ppar, stp, ecL_cur);
#endif
        form_L_diag(ecL_prev, ecL_cur);
        // update offsets
        ecL_cur = &ecL_cur[MATRIX_SIZE_3x3];
        ecL_prev = &ecL_prev[MATRIX_SIZE_3x3];
    }
}


/**
 * @brief Solve system ecL * x = b using forward substitution.
 *
 * @param[in] N number of states in the preview window
 * @param[in,out] x vector "b" as input, vector "x" as output
 *                  (N * #NUM_STATE_VAR)
 */
void matrix_ecL_as::solve_forward(const int N, double *x)
{
    int i;
    double *xc = x; // 6 current elements of x
    double *xp; // 6 elements of x computed on the previous iteration
    double *cur_ecL = &ecL[0];  // lower triangular matrix lying on the 
                                // diagonal of L
    double *prev_ecL;   // upper triangular matrix lying to the left from
                        // cur_ecL at the same level of L


    // compute the first 6 elements using forward substitution
    xc[0] /= cur_ecL[0];
    xc[3] /= cur_ecL[0];

    xc[1] -= xc[0] * cur_ecL[1];
    xc[1] /= cur_ecL[4];

    xc[4] -= xc[3] * cur_ecL[1];
    xc[4] /= cur_ecL[4];

    xc[2] -= xc[0] * cur_ecL[2] + xc[1] * cur_ecL[5];
    xc[2] /= cur_ecL[8];

    xc[5] -= xc[3] * cur_ecL[2] + xc[4] * cur_ecL[5];
    xc[5] /= cur_ecL[8];


    for (i = 1; i < N; i++)
    {
        // switch to the next level of L / next 6 elements
        xp = xc;
        xc = &xc[NUM_STATE_VAR];

        prev_ecL = &cur_ecL[MATRIX_SIZE_3x3];
        cur_ecL = &cur_ecL[2 * MATRIX_SIZE_3x3];


        // update the right part of the equation and compute elements
        xc[0] -= xp[0] * prev_ecL[0] + xp[1] * prev_ecL[3] + xp[2] * prev_ecL[6];
        xc[0] /= cur_ecL[0];

        xc[3] -= xp[3] * prev_ecL[0] + xp[4] * prev_ecL[3] + xp[5] * prev_ecL[6];
        xc[3] /= cur_ecL[0];


        xc[1] -= xp[1] * prev_ecL[4] + xp[2] * prev_ecL[7] + xc[0] * cur_ecL[1];
        xc[1] /= cur_ecL[4];

        xc[4] -= xp[4] * prev_ecL[4] + xp[5] * prev_ecL[7] + xc[3] * cur_ecL[1];
        xc[4] /= cur_ecL[4];


        xc[2] -= xp[2] * prev_ecL[8] + xc[0] * cur_ecL[2] + xc[1] * cur_ecL[5];
        xc[2] /= cur_ecL[8];

        xc[5] -= xp[5] * prev_ecL[8] + xc[3] * cur_ecL[2] + xc[4] * cur_ecL[5];
        xc[5] /= cur_ecL[8];
    }
}


/**
 * @brief Solve system ecL' * x = b using backward substitution.
 *
 * @param[in] N number of states in the preview window
 * @param[in,out] x vector "b" as input, vector "x" as output.
 */
void matrix_ecL_as::solve_backward (const int N, double *x)
{
    int i;
    double *xc = & x[(N-1)*NUM_STATE_VAR]; // current 6 elements of result
    double *xp; // 6 elements computed on the previous iteration
    
    // elements of these matrices accessed as if they were transposed
    // lower triangular matrix lying on the diagonal of L
    double *cur_ecL = &ecL[2 * (N - 1) * MATRIX_SIZE_3x3];
    // upper triangular matrix lying to the right from cur_ecL at the same level of L'
    double *prev_ecL; 


    // compute the last 6 elements using backward substitution
    xc[2] /= cur_ecL[8];
    xc[5] /= cur_ecL[8];

    xc[1] -= xc[2] * cur_ecL[5];
    xc[1] /= cur_ecL[4];
    xc[4] -= xc[5] * cur_ecL[5];
    xc[4] /= cur_ecL[4];

    xc[0] -= xc[2] * cur_ecL[2] + xc[1] * cur_ecL[1];
    xc[0] /= cur_ecL[0];
    xc[3] -= xc[5] * cur_ecL[2] + xc[4] * cur_ecL[1];
    xc[3] /= cur_ecL[0];


    for (i = N-2; i >= 0 ; i--)
    {
        xp = xc;
        xc = & x[i*NUM_STATE_VAR];

        cur_ecL = &ecL[2 * i * MATRIX_SIZE_3x3];
        prev_ecL = &cur_ecL[MATRIX_SIZE_3x3];


        // update the right part of the equation and compute elements
        xc[2] -= xp[0] * prev_ecL[6] + xp[1] * prev_ecL[7] + xp[2] * prev_ecL[8];
        xc[2] /= cur_ecL[8];

        xc[5] -= xp[3] * prev_ecL[6] + xp[4] * prev_ecL[7] + xp[5] * prev_ecL[8];
        xc[5] /= cur_ecL[8];


        xc[1] -= xp[0] * prev_ecL[3] + xp[1] * prev_ecL[4] + xc[2] * cur_ecL[5];
        xc[1] /= cur_ecL[4];

        xc[4] -= xp[3] * prev_ecL[3] + xp[4] * prev_ecL[4] + xc[5] * cur_ecL[5];
        xc[4] /= cur_ecL[4];


        xc[0] -= xp[0] * prev_ecL[0] + xc[2] * cur_ecL[2] + xc[1] * cur_ecL[1];
        xc[0] /= cur_ecL[0];

        xc[3] -= xp[3] * prev_ecL[0] + xc[5] * cur_ecL[2] + xc[4] * cur_ecL[1];
        xc[3] /= cur_ecL[0];
    }
}
