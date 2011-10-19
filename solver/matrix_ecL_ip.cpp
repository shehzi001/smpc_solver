/** 
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:08 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/

#include "matrix_ecL_ip.h"

#include <cmath> // sqrt


/****************************************
 * FUNCTIONS 
 ****************************************/

//==============================================
// constructors / destructors


matrix_ecL_ip::matrix_ecL_ip (const int N)
{
    ecL = new double[MATRIX_SIZE_6x6*N + MATRIX_SIZE_6x6*(N-1)]();

    M = new double[MATRIX_SIZE_6x6];
    MAT = new double[MATRIX_SIZE_6x6];
}


matrix_ecL_ip::~matrix_ecL_ip()
{
    if (ecL != NULL)
        delete ecL;

    if (M != NULL)
        delete M;
    if (MAT != NULL)
        delete MAT;
}

//==============================================



/**
 * @brief Forms M = R*inv(hess_phi)*R'.
 *
 * @param[in] sinA sin of rotation angle.
 * @param[in] cosA cos of rotation angle.
 * @param[in] i2Q a vector of three repeating diagonal elements of inv(Q)
 * @param[in] i2hess a 2*N vector of diagonal elements of hess_phi 
 *                  (indicies of these elements are 1:3:N*NUM_STATE_VAR)
 *
 * @attention Only elements lying below the main diagonal of 4x4 matrix
 *            are initialized (other elements are not unique).
 */
void matrix_ecL_ip::form_M (
        const double sinA,
        const double cosA,
        const double *i2Q,
        const double* i2hess)
{
    /*
     * Numbers mark identical elements of M
     * 1  5
     *  2 
     *   3
     * 5  4
     *     2
     *      3
     */
    // diagonal elements
    M[0] = i2hess[0]*cosA*cosA + i2hess[1]*sinA*sinA;
    /*M[28] =*/ M[7] = i2Q[1];
    /*M[35] =*/ M[14] = i2Q[2];
    M[21] = i2hess[0]*sinA*sinA + i2hess[1]*cosA*cosA;

    // symmetric nondiagonal elements
    /*M[18] =*/ M[3] = (i2hess[0] - i2hess[1])*cosA*sinA;
}



/**
 * @brief Performs Cholesky decomposition of a matrix.
 *
 * @param[in,out] mx a pointer to a 6x6 matrix, the result is 
 *                    stored in the same place.
 *
 * @attention Only the elements below the main diagonal are used
 *              in conmputations.
 */
void matrix_ecL_ip::chol_dec (double *mx)
{
    double *curel;

    for (int i = 0; i < MATRIX_SIDE_6x6; i++) // row
    {
        // non-diagonal elements
        for (int j = 0; j < i; j++) // column
        {
            curel = &mx[j*MATRIX_SIDE_6x6 + i];
            for (int k = 0; k < j; k++) // column
            {
                *curel -= mx[k*MATRIX_SIDE_6x6 + i] * mx[k*MATRIX_SIDE_6x6 + j];
            }

            *curel /= mx[j*MATRIX_SIDE_6x6 + j];
        }

        // diagonal element
        curel = &mx[i*MATRIX_SIDE_6x6 + i];
        for (int k = 0; k < i; k++) // column
        {
            *curel -= mx[k*MATRIX_SIDE_6x6 + i] * mx[k*MATRIX_SIDE_6x6 + i];
        }
        *curel = sqrt(*curel);
    }
}



/**
 * @brief Forms matrix MBiPB = M + B * inv(2*P) * B.
 *
 * @param[in] B a vector of 3 elements.
 * @param[in] i2P 0.5 * inv(P) (only one number)
 * @param[out] result result.
 *
 * @attention Only elements lying below the main diagonal of 4x4 matrix
 *            are initialized (other elements are not unique).
 */
void matrix_ecL_ip::form_MBiPB (const double *B, const double i2P, double *result)
{
    // diagonal elements
    result[0]  = i2P * B[0]*B[0] + M[0];
    /*result[28] =*/ result[7]  = i2P * B[1]*B[1] + M[7];
    /*result[35] =*/ result[14] = i2P * B[2]*B[2] + M[14];
    result[21] = i2P * B[0]*B[0] + M[21];

    // symmetric elements (no need to initialize all of them)
    /*result[22] =*/ result[1] = i2P * B[0]*B[1];
    /*result[23] =*/ result[2] = i2P * B[0]*B[2];
    /*result[29] =*/ result[8] = i2P * B[1]*B[2];
    result[3] = M[3];
}



/**
 * @brief Forms matrix MAT = M * A'
 *
 * @param[in] A3 4th and 7th elements of A.
 * @param[in] A6 6th element of A.
 */
void matrix_ecL_ip::form_MAT (const double A3, const double A6)
{
    // 1st column
    MAT[0] = M[0];
    MAT[1] = A3 * M[7];
    MAT[2] = A6 * M[14];
    MAT[3] = M[3];

    // 2nd column
    MAT[7] = M[7];
    MAT[8] = A3 * M[14];

    // 3rd column
    MAT[14] = M[14];

    // 4th column
    MAT[21] = M[21];
    MAT[22] = MAT[1];
    MAT[23] = MAT[2];

    // 5th column
    MAT[28] = MAT[7];
    MAT[29] = MAT[8];

    // 6th column
    MAT[35] = MAT[14];
}



/**
 * @brief Forms a 6x6 matrix L(k+1, k), which lies below the 
 *  diagonal of L.
 *
 * @param[in] ecLp previous matrix lying on the diagonal of L
 * @param[in] ecLc the result is stored here
 */
void matrix_ecL_ip::form_L_non_diag(const double *ecLp, double *ecLc)
{
    /* 
     * L(k,k)   * L(k+1,k)' = -M*A'
     *
     * x         x  x       x  x
     * xx        xx x       xx
     * xxx     * xxxx   =   xxx
     * xxxx      xxxx       x  x
     * xxxxx     xxxxx         xx
     * xxxxxx    xxxxxx        xxx
     */


    // copy MAT' to ecLc
    ecLc[0]  = -MAT[0];
    ecLc[27] = ecLc[6]  = -MAT[1];
    ecLc[33] = ecLc[12] = -MAT[2];

    ecLc[3]  = ecLc[18] = -MAT[3];

    ecLc[28] = ecLc[7]  = -MAT[7];
    ecLc[34] = ecLc[13] = -MAT[8];

    ecLc[35] = ecLc[14] = -MAT[14];

    ecLc[21] = -MAT[21];

    ecLc[21] = -MAT[21];

    // reset elements
    ecLc[9]  = ecLc[15] = 0;


    for (int i = 0; i < MATRIX_SIDE_6x6; i++) // row of L(k+1,k)
    {
        for (int j = 0; j < MATRIX_SIDE_6x6; j++) // element in the row of L(k+1,k)
        {
            double *curel = &ecLc[j*MATRIX_SIDE_6x6 + i];
            for (int k = 0; k < j; k++) // elements in a row of L(k,k)
            {
                *curel -= ecLc[k*MATRIX_SIDE_6x6 + i] * ecLp[k*MATRIX_SIDE_6x6 + j];
            }

            *curel /= ecLp[j*MATRIX_SIDE_6x6 + j];
        }
    }
}


/**
 * @brief Forms a 6x6 matrix L(k+1, k+1), which lies below the 
 *  diagonal of L.
 *
 * @param[in,out] ecLc MBiPB as input / result
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_ip::form_L_diag(double *ecLc)
{
    // finish initialization of MBiPB
    ecLc[28] = ecLc[7];
    ecLc[35] = ecLc[14];
    // symmetric elements (no need to initialize all of them)
    ecLc[22] = ecLc[1];
    ecLc[23] = ecLc[2];
    ecLc[29] = ecLc[8];
    // reset elements
    ecLc[4] = ecLc[5] = ecLc[9] = ecLc[10] = ecLc[11] = 
        ecLc[15] = ecLc[16] = ecLc[17] = 0;


    // chol (L(k+1,k+1))
    chol_dec (ecLc);
}



/**
 * @brief Forms matrix AMATMBiPB = 
 *  A * M * A' + 0.5 * M + 0.5 * B * inv(P) * B
 *
 * @param[in] A3 4th and 7th elements of A (A is represented by two identical 3x3 matrices).
 * @param[in] A6 6th element of A (A is represented by two identical 3x3 matrices).
 * @param[in,out] result MBiPB as input / result
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_ip::form_AMATMBiPB(const double A3, const double A6, double *result)
{
    // 1st,4th column
    result[0] += MAT[0] + A3*MAT[1] + A6*MAT[2];
    result[1] +=             MAT[1] + A3*MAT[2];
    result[2] +=                         MAT[2];
    result[3] +=                         MAT[3];

    // 2nd column
    // symmetric elements are not initialized
    result[7] += MAT[7] +  A3*MAT[8];
    result[8] +=              MAT[8];

    // 6th column
    result[14] += MAT[14];
    result[35] = result[14];

    // 4th column
    result[21] += MAT[21] + A3*MAT[1] + A6*MAT[2];
    result[22] = result[1];
    result[23] = result[2];

    // 5th column
    result[28] = result[7];
    result[29] = result[8];

    // 6th column
    result[35] = result[14];

    // reset elements
    result[4] = result[5] = result[9] = result[10] = result[11] = 
        result[15] = result[16] = result[17] = 0;
}



/**
 * @brief Forms a 6x6 matrix L(k+1, k+1), which lies on the main
 *  diagonal of L.
 *
 * @param[in] ecLp a 6x6 matrix lying to the left from ecLc on the same 
 *                 level of L
 * @param[in,out] ecLc AMATMBiPB as input / the result is stored here
 *
 * @attention Only the elements below the main diagonal are initialized.
 */
void matrix_ecL_ip::form_L_diag (const double *ecLp, double *ecLc)
{
    /* - L(k+1,k) * L(k+1,k)' + A*M*A' + MBiPB
     * xxxxxx   x  x
     *  xxxxx   xx x
     *   xxxx   xxxx
     * xxxxxx * xxxx
     *     xx   xxxxx
     *      x   xxxxxx
     */
    // - L(k+1,k) * L(k+1,k)'
    for (int i = 0; i < MATRIX_SIDE_6x6; i++) // row
    {
        for (int j = 0; j <= i; j++) // column
        {
            double *curel = &ecLc[i + MATRIX_SIDE_6x6*j];
            for (int k = 0; k < MATRIX_SIDE_6x6; k++) // elements in a row
            {
                *curel -= ecLp[i + MATRIX_SIDE_6x6*k] * ecLp[j + MATRIX_SIDE_6x6*k];
            }
        }
    }

    // chol (L(k+1,k+1))
    chol_dec (ecLc);
}



/**
 * @brief Builds matrix L.
 *
 * @param[in] ppar      parameters.
 * @param[in] i2hess    2*N diagonal elements of inverted hessian.
 */
void matrix_ecL_ip::form (const problem_parameters* ppar, const double *i2hess)
{
    int i;
    state_parameters stp;

    stp = ppar->spar[0];

    // the first matrix on diagonal
    form_M (stp.sin, stp.cos, ppar->i2Q, i2hess);
    form_MBiPB (stp.B, ppar->i2P, ecL);
    form_L_diag (ecL);

    // offsets
    double *ecL_cur = &ecL[MATRIX_SIZE_6x6];
    double *ecL_prev = &ecL[0];
    for (i = 1; i < ppar->N; i++)
    {
        stp = ppar->spar[i];

        // form all matrices
        form_MAT (stp.A3, stp.A6);
        form_L_non_diag (ecL_prev, ecL_cur);

        // update offsets
        ecL_cur = &ecL_cur[MATRIX_SIZE_6x6];
        ecL_prev = &ecL_prev[MATRIX_SIZE_6x6];


        i2hess = &i2hess[2];
        form_M (stp.sin, stp.cos, ppar->i2Q, i2hess);
        form_MBiPB (stp.B, ppar->i2P, ecL_cur);
        form_AMATMBiPB(stp.A3, stp.A6, ecL_cur);
        form_L_diag(ecL_prev, ecL_cur);

        // update offsets
        ecL_cur = &ecL_cur[MATRIX_SIZE_6x6];
        ecL_prev = &ecL_prev[MATRIX_SIZE_6x6];
    }
}



/**
 * @brief Solve system ecL * x = b using forward substitution.
 *
 * @param[in] N number of states in the preview window
 * @param[in,out] x vector "b" as input, vector "x" as output
 *                  (N * #NUM_STATE_VAR)
 */
void matrix_ecL_ip::solve_forward(const int N, double *x)
{
    int i,j,k;
    double *xc = x; // 6 current elements of x
    double *xp; // 6 elements of x computed on the previous iteration
    double *ecL_cur = &ecL[0];  // lower triangular matrix lying on the 
                                // diagonal of L
    double *ecL_prev;   // upper triangular matrix lying to the left from
                        // ecL_cur at the same level of L


    // compute the first 6 elements using forward substitution
    for (j = 0; j < MATRIX_SIZE_6x6; j++) // row
    {
        for (k = 0; k < j; k++) // column
        {
            xc[j] -= xc[k]*ecL_cur[j+k*MATRIX_SIZE_6x6];
        }
        xc[j] /= ecL_cur[j+j*MATRIX_SIZE_6x6];
    }


    for (i = 1; i < N; i++)
    {
        // switch to the next level of L / next 6 elements
        xp = xc;
        xc = &xc[NUM_STATE_VAR];

        ecL_prev = &ecL_cur[MATRIX_SIZE_6x6];
        ecL_cur = &ecL_prev[MATRIX_SIZE_6x6];

        
        // update the right part of the equation
        for (j = 0; j < MATRIX_SIZE_6x6; j++) // column
        {
            for (k = 0; k < MATRIX_SIZE_6x6; k++) // row
            {
                xc[k] -= xp[j]*ecL_prev[j*MATRIX_SIZE_6x6 + k];
            }
        }

        // forward substitution
        for (j = 0; j < MATRIX_SIZE_6x6; j++) // row
        {
            for (k = 0; k < j; k++) // column
            {
                xc[j] -= xc[k]*ecL_cur[j+k*MATRIX_SIZE_6x6];
            }
            xc[j] /= ecL_cur[j+j*MATRIX_SIZE_6x6];
        }
    }
}


/**
 * @brief Solve system ecL' * x = b using backward substitution.
 *
 * @param[in] N number of states in the preview window
 * @param[in,out] x vector "b" as input, vector "x" as output.
 */
void matrix_ecL_ip::solve_backward (const int N, double *x)
{
    int i,j,k;
    double *xc = & x[(N-1)*NUM_STATE_VAR]; // current 6 elements of result
    double *xp; // 6 elements computed on the previous iteration
    
    // elements of these matrices accessed as if they were transposed
    // lower triangular matrix lying on the diagonal of L
    double *ecL_cur = &ecL[2 * (N - 1) * MATRIX_SIZE_6x6];
    // upper triangular matrix lying to the right from ecL_cur at the same level of L'
    double *ecL_prev; 


    // compute the last 6 elements using backward substitution
    for (j = MATRIX_SIZE_6x6-1; j >= 0; j--) // row
    {
        for (k = MATRIX_SIZE_6x6-1; k > j; k--) // column
        {
            xc[j] -= xc[k]*ecL_cur[k+j*MATRIX_SIZE_6x6];
        }
        xc[j] /= ecL_cur[j+j*MATRIX_SIZE_6x6];
    }


    for (i = N-2; i >= 0 ; i--)
    {
        xp = xc;
        xc = & x[i*NUM_STATE_VAR];

        ecL_cur = &ecL[2 * i * MATRIX_SIZE_6x6];
        ecL_prev = &ecL_cur[MATRIX_SIZE_6x6];


        // update the right part of the equation
        for (j = 0; j < MATRIX_SIZE_6x6; j++) // row
        {
            for (k = 0; k < MATRIX_SIZE_6x6; k++) // column
            {
                xc[j] -= xp[k]*ecL_prev[j*MATRIX_SIZE_6x6 + k];
            }
        }

        // backward substitution
        for (j = MATRIX_SIZE_6x6-1; j >= 0; j--) // row
        {
            for (k = MATRIX_SIZE_6x6-1; k > j; k--) // column
            {
                xc[j] -= xc[k]*ecL_cur[k+j*MATRIX_SIZE_6x6];
            }
            xc[j] /= ecL_cur[j+j*MATRIX_SIZE_6x6];
        }
    }
}
