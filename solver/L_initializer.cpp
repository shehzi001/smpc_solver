/** 
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:08 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/

#include "L_initializer.h"

#include <cmath> // sqrt


/****************************************
 * FUNCTIONS 
 ****************************************/

//==============================================
// constructors / destructors

L_initializer::L_initializer ()
{
    iQBiPB = new double[MATRIX_SIZE];
    iQAT = new double[MATRIX_SIZE];
    AiQATiQBiPB = new double[MATRIX_SIZE];
}


L_initializer::~L_initializer()
{
    if (iQBiPB != NULL)
        delete iQBiPB;
    if (iQAT != NULL)
        delete iQAT;
    if (AiQATiQBiPB != NULL)
        delete AiQATiQBiPB;
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
void L_initializer::chol_dec (double *mx9)
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
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void L_initializer::form_iQBiPB (double *B, double *i2Q, double i2P)
{
    // diagonal elements
    iQBiPB[0] = i2P * B[0]*B[0] + i2Q[0];
    iQBiPB[4] = i2P * B[1]*B[1] + i2Q[1];
    iQBiPB[8] = i2P * B[2]*B[2] + i2Q[2];

    // symmetric elements (no need to initialize all of them)
    iQBiPB[1] = /*iQBiPB[3] =*/ i2P * B[0]*B[1];
    iQBiPB[2] = /*iQBiPB[6] =*/ i2P * B[0]*B[2];
    iQBiPB[5] = /*iQBiPB[7] =*/ i2P * B[1]*B[2];
}



/**
 * @brief Forms matrix iQAT = 0.5 * inv (Q) * A'
 *
 * @param[in] T 4th and 7th elements of A.
 * @param[in] A6 6th element of A.
 * @param[in] i2Q a vector of 3 elements, which contains
 *              diagonal elements of 0.5*inv(Q).
 */
void L_initializer::form_iQAT (double T, double A6, double *i2Q)
{
    iQAT[0] = i2Q[0];
    iQAT[1] = T * i2Q[1];
    iQAT[2] = A6 * i2Q[2];
    iQAT[4] = i2Q[1];
    iQAT[5] = T * i2Q[2];
    iQAT[8] = i2Q[2];
}


/**
 * @brief Forms matrix AiQATiQBiPB = 
 *  A * inv(Q) * A' + 0.5 * inv(Q) + 0.5 * B * inv(P) * B
 *
 * @param[in] T 4th and 7th elements of A.
 * @param[in] A6 6th element of A.
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void L_initializer::form_AiQATiQBiPB (double T, double A6)
{
    // 1st column
    AiQATiQBiPB[0] = iQBiPB[0] + iQAT[0] + T*iQAT[1] + A6*iQAT[2];
    AiQATiQBiPB[1] = iQBiPB[1] +             iQAT[1] +  T*iQAT[2];
    AiQATiQBiPB[2] = iQBiPB[2] +                          iQAT[2];

    // 2nd column
    // symmetric elements are not initialized
//    AiQATiQBiPB[3] = iQBiPB[3] + T*iQAT[4] + A6*iQAT[5];
    AiQATiQBiPB[4] = iQBiPB[4] +   iQAT[4] +  T*iQAT[5];
    AiQATiQBiPB[5] = iQBiPB[5] +                iQAT[5];

    // 3rd column
//    AiQATiQBiPB[6] = iQBiPB[6] + A6*iQAT[8];
//    AiQATiQBiPB[7] = iQBiPB[7] +  T*iQAT[8];
    AiQATiQBiPB[8] = iQBiPB[8] +    iQAT[8];
}



/**
 * @brief Forms a 3x3 matrix L(k+1, k), which lies below the 
 *  diagonal of L.
 *
 * @param[in] ecLp previous matrix lying on the diagonal of L
 * @param[in] ecLc the result is stored here
 */
void L_initializer::form_L_non_diag(double *ecLp, double *ecLc)
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
 * @brief Forms a 3x3 matrix L(k+1, k+1), which lies below the 
 *  diagonal of L.
 *
 * @param[in] ecLp upper triangular matrix matrix lying to the left
 *                 from ecLc on the same level of L
 * @param[in] ecLc the result is stored here
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void L_initializer::form_L_diag(double *ecLp, double *ecLc)
{
    // the first matrix L(0,0) is computed differently
    if (ecLp == NULL)
    {
        //0.5*inv(Q) + 0.5*B*inv(P)*B'
        // matrix is symmetric (no need to initialize all elements)
        ecLc[0] = iQBiPB[0];
        ecLc[1] = iQBiPB[1];
        ecLc[2] = iQBiPB[2];
        ecLc[4] = iQBiPB[4];
        ecLc[5] = iQBiPB[5];
        ecLc[8] = iQBiPB[8];
    }
    else
    {
    // L(k+1,k+1) = (- L(k+1,k) * L(k+1,k)') + (A * inv(Q) * A' + inv(Q) + B * inv(P) * B)
        // diagonal elements
        ecLc[0] = -(ecLp[0]*ecLp[0] + ecLp[3]*ecLp[3] + ecLp[6]*ecLp[6]) + AiQATiQBiPB[0];
        ecLc[4] = -(ecLp[4]*ecLp[4] + ecLp[7]*ecLp[7]) + AiQATiQBiPB[4];
        ecLc[8] = -(ecLp[8]*ecLp[8]) + AiQATiQBiPB[8];
        // symmetric nondiagonal elements (no need to initialize all of them)
        ecLc[1] = /*ecLc[3] =*/ -(ecLp[3]*ecLp[4] + ecLp[6]*ecLp[7]) + AiQATiQBiPB[1];
        ecLc[2] = /*ecLc[6] =*/ -ecLp[6]*ecLp[8] + AiQATiQBiPB[2];
        ecLc[5] = /*ecLc[7] =*/ -ecLp[7]*ecLp[8] + AiQATiQBiPB[5];
    }

    // chol (L(k+1,k+1))
    chol_dec (ecLc);
}



/**
 * @brief Builds matrix L.
 *
 * @param[in] csp parameters.
 * @param[in] N number of states in preview window.
 * @param[out] ecL the memory allocated for L.
 */
void L_initializer::form_L(chol_solve_param csp, int N, double *ecL)
{
    int i;
    int cur_offset;
    int prev_offset;

    double T = csp.T[0];
    double T2 = T*T/2;
    double B[3] = {T2*T/3 - csp.h[0]*T, T2, T};

    // form all matrices
    form_iQBiPB (B, csp.i2Q, csp.i2P);
#ifndef QPAS_VARIABLE_T_h
    double A6 = T2;
    form_iQAT (T, A6, csp.i2Q);
    form_AiQATiQBiPB (T, A6);
#endif


    // the first matrix on diagonal
    form_L_diag(NULL, ecL);

    // offsets
    cur_offset = MATRIX_SIZE;
    prev_offset = 0;
    for (i = 1; i < N; i++)
    {
#ifdef QPAS_VARIABLE_T_h
        T = csp.T[i];
        T2 = T*T/2;
        B[0] = T2*T/3 - csp.h[i]*T;
        B[1] = T2;
        B[2] = T;
        double A6 = T2 - csp.dh[i-1];

        // form all matrices
        form_iQBiPB (B, csp.i2Q, csp.i2P);
        form_iQAT (T, A6, csp.i2Q);
        form_AiQATiQBiPB (T, A6);
#endif

        // form (b), (d), (f) ... 
        form_L_non_diag(&ecL[prev_offset], &ecL[cur_offset]);
        // update offsets
        cur_offset += MATRIX_SIZE;
        prev_offset += MATRIX_SIZE;

        // form (c), (e), (g) ...
        form_L_diag(&ecL[prev_offset], &ecL[cur_offset]);
        // update offsets
        cur_offset += MATRIX_SIZE;
        prev_offset += MATRIX_SIZE;
    }
}
