/** 
 * @file
 * @author Alexander Sherikov
 * @date 28.08.2011 13:43:08 MSD
 */


/****************************************
 * INCLUDES 
 ****************************************/

#include "ip_matrix_ecL.h"

#include <cmath> // sqrt


/****************************************
 * FUNCTIONS 
 ****************************************/
namespace IP
{
    //==============================================
    // constructors / destructors


    matrix_ecL::matrix_ecL (const int N)
    {
        ecL = new double[MATRIX_SIZE_6x6*N + MATRIX_SIZE_6x6*(N-1)]();

        M = new double[MATRIX_SIZE_6x6];
        MAT = new double[MATRIX_SIZE_6x6];
    }


    matrix_ecL::~matrix_ecL()
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
     *                  (indicies of these elements are 1:3:N*SMPC_NUM_STATE_VAR)
     *
     * @attention Only elements lying below the main diagonal of 4x4 matrix
     *            are initialized (other elements are not unique).
     */
    void matrix_ecL::form_M (
            const double sinA,
            const double cosA,
            const double *i2Q,
            const double* i2hess)
    {
        /*      R        *       Q        *       R'      =      M
         * |c    -s    |   |a1          |   |c     s    |   |a1cc+a2ss     a1cs-a2cs    |
         * |  1        |   |   b        |   |  1        |   |          b                |
         * |    1      |   |     g      |   |    1      |   |            g              |
         * |s     c    |   |      a2    |   |-s    c    |   |a1cs-a2cs     a1ss+a2cc    |
         * |        1  |   |         b  |   |        1  |   |                        b  |
         * |          1|   |           g|   |          1|   |                          g|
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
    void matrix_ecL::chol_dec (double *mx)
    {
        mx[0] = sqrt(mx[0]);
        mx[1] /= mx[0];
        mx[2] /= mx[0];
        mx[3] /= mx[0];
        mx[4] /= mx[0];
        mx[5] /= mx[0];

        mx[7]  = sqrt(mx[7] - mx[1]*mx[1]);
        mx[8]  = (mx[8]  - mx[1]*mx[2])/mx[7];
        mx[9]  = (mx[9]  - mx[1]*mx[3])/mx[7];
        mx[10] = (mx[10] - mx[1]*mx[4])/mx[7];
        mx[11] = (mx[11] - mx[1]*mx[5])/mx[7];

        mx[14] = sqrt(mx[14] - mx[2]*mx[2] - mx[8]*mx[8]);
        mx[15] = (mx[15] - mx[2]*mx[3] - mx[8]*mx[9])/mx[14];
        mx[16] = (mx[16] - mx[2]*mx[4] - mx[8]*mx[10])/mx[14];
        mx[17] = (mx[17] - mx[2]*mx[5] - mx[8]*mx[11])/mx[14];

        mx[21] = sqrt(mx[21] - mx[3]*mx[3] - mx[9]*mx[9] - mx[15]*mx[15]);
        mx[22] = (mx[22] - mx[3]*mx[4] - mx[9]*mx[10] - mx[15]*mx[16])/mx[21];
        mx[23] = (mx[23] - mx[3]*mx[5] - mx[9]*mx[11] - mx[15]*mx[17])/mx[21];

        mx[28] = sqrt(mx[28] - mx[4]*mx[4] - mx[10]*mx[10] - mx[16]*mx[16] - mx[22]*mx[22]);
        mx[29] = (mx[29] - mx[4]*mx[5] - mx[10]*mx[11] - mx[16]*mx[17] - mx[22]*mx[23])/mx[28];

        mx[35] = sqrt(mx[35] - mx[5]*mx[5] - mx[11]*mx[11] - mx[17]*mx[17] - mx[23]*mx[23] - mx[29]*mx[29]);
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
    void matrix_ecL::form_MBiPB (const double *B, const double i2P, double *result)
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
    void matrix_ecL::form_MAT (const double A3, const double A6)
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
    void matrix_ecL::form_L_non_diag(const double *ecLp, double *ecLc)
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

        ecLc[0] = -MAT[0]/ecLp[0];
        ecLc[1] = 0;
        ecLc[2] = 0;
        ecLc[3] = -MAT[3]/ecLp[0]; // MAT[3] = MAT[18]
        ecLc[4] = 0;
        ecLc[5] = 0;

        ecLc[6]  = (-MAT[1] - ecLc[0] * ecLp[1]) / ecLp[7];
        ecLc[7]  = -MAT[7] / ecLp[7];
        ecLc[8]  = 0;
        ecLc[9]  = (0 - ecLc[3] * ecLp[1]) / ecLp[7];
        ecLc[10] = 0;
        ecLc[11] = 0;

        ecLc[12] = (-MAT[2] - ecLc[0] * ecLp[2] - ecLc[6] * ecLp[8]) / ecLp[14];
        ecLc[13] = (-MAT[8] - ecLc[7] * ecLp[8]) / ecLp[14]; 
        ecLc[14] = -MAT[14] / ecLp[14];
        ecLc[15] = (0 - ecLc[3] * ecLp[2] - ecLc[9] * ecLp[8]) / ecLp[14];
        ecLc[16] = 0;
        ecLc[17] = 0;

        ecLc[18] = (-MAT[3] - ecLc[0]*ecLp[3] - ecLc[6]*ecLp[9] - ecLc[12]*ecLp[15]) / ecLp[21];
        ecLc[19] = (0 - ecLc[7]*ecLp[9] - ecLc[13]*ecLp[15])/ecLp[21];
        ecLc[20] = (0 - ecLc[14]*ecLp[15]) / ecLp[21];
        ecLc[21] = (-MAT[21] - ecLc[3]*ecLp[3] - ecLc[9]*ecLp[9] - ecLc[15]*ecLp[15]) / ecLp[21];
        ecLc[22] = 0;
        ecLc[23] = 0;

        ecLc[24] = (0 - ecLc[0]*ecLp[4] - ecLc[6]*ecLp[10] - ecLc[12]*ecLp[16] - ecLc[18]*ecLp[22]) / ecLp[28];
        ecLc[25] = (0 - ecLc[7]*ecLp[10] - ecLc[13]*ecLp[16] - ecLc[19]*ecLp[22]) / ecLp[28];
        ecLc[26] = (0 - ecLc[14]*ecLp[16] - ecLc[20]*ecLp[22]) / ecLp[28];
        ecLc[27] = (-MAT[22] - ecLc[3]*ecLp[4] - ecLc[9]*ecLp[10] - ecLc[15]*ecLp[16] - ecLc[21]*ecLp[22]) / ecLp[28];
        ecLc[28] = -MAT[28] / ecLp[28];
        ecLc[29] = 0;

        ecLc[30] = (0 - ecLc[0]*ecLp[5] - ecLc[6]*ecLp[11] - ecLc[12]*ecLp[17] - ecLc[18]*ecLp[23] - ecLc[24]*ecLp[29]) / ecLp[35];
        ecLc[31] = (0 - ecLc[7]*ecLp[11] - ecLc[13]*ecLp[17] - ecLc[19]*ecLp[23] - ecLc[25]*ecLp[29]) / ecLp[35];
        ecLc[32] = (0 - ecLc[14]*ecLp[17] - ecLc[20]*ecLp[23] - ecLc[26]*ecLp[29]) / ecLp[35];
        ecLc[33] = (-MAT[23] - ecLc[3]*ecLp[5] - ecLc[9]*ecLp[11] - ecLc[15]*ecLp[17] - ecLc[21]*ecLp[23] - ecLc[27]*ecLp[29]) / ecLp[35];
        ecLc[34] = (-MAT[29] - ecLc[28]*ecLp[29])/ ecLp[35];
        ecLc[35] = -MAT[35] / ecLp[35];
    }


    /**
     * @brief Forms a 6x6 matrix L(k+1, k+1), which lies below the 
     *  diagonal of L.
     *
     * @param[in,out] ecLc MBiPB as input / result
     *
     * @attention Only the elements below the main diagonal are initialized.
     */
    void matrix_ecL::form_L_diag(double *ecLc)
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
    void matrix_ecL::form_AMATMBiPB(const double A3, const double A6, double *result)
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
    }



    /**
     * @brief Forms a 6x6 matrix L(k+1, k+1), which lies on the main
     *  diagonal of L.
     *
     * @param[in] p a 6x6 matrix lying to the left from ecLc on the same 
     *                 level of L
     * @param[in,out] ecLc AMATMBiPB as input / the result is stored here
     *
     * @attention Only the elements below the main diagonal are initialized.
     */
    void matrix_ecL::form_L_diag (const double *p, double *ecLc)
    {
        /* - L(k+1,k) * L(k+1,k)' + A*M*A' + MBiPB
         * xxxxxx   x  x
         *  xxxxx   xx x
         *   xxxx   xxxx
         * xxxxxx * xxxx
         *     xx   xxxxx
         *      x   xxxxxx
         */

        ecLc[0]  += - p[0] *p[0]  - p[6] *p[6]  - p[12]*p[12] - p[18]*p[18] - p[24]*p[24] - p[30]*p[30];
        ecLc[1]  +=               - p[6] *p[7]  - p[12]*p[13] - p[18]*p[19] - p[24]*p[25] - p[30]*p[31]; 
        ecLc[2]  +=                             - p[12]*p[14] - p[18]*p[20] - p[24]*p[26] - p[30]*p[32];
        ecLc[3]  += - p[0] *p[3]  - p[6] *p[9]  - p[12]*p[15] - p[18]*p[21] - p[24]*p[27] - p[30]*p[33];
        ecLc[4]   =                                                         - p[24]*p[28] - p[30]*p[34];
        ecLc[5]   =                                                                       - p[30]*p[35];
                   
        ecLc[7]  +=               - p[7] *p[7]  - p[13]*p[13] - p[19]*p[19] - p[25]*p[25] - p[31]*p[31];
        ecLc[8]  +=                             - p[13]*p[14] - p[19]*p[20] - p[25]*p[26] - p[31]*p[32];
        ecLc[9]   =               - p[7] *p[9]  - p[13]*p[15] - p[19]*p[21] - p[25]*p[27] - p[31]*p[33];
        ecLc[10]  =                                                         - p[25]*p[28] - p[31]*p[34];
        ecLc[11]  =                                                                       - p[31]*p[35];
                   
        ecLc[14] +=                             - p[14]*p[14] - p[20]*p[20] - p[26]*p[26] - p[32]*p[32];
        ecLc[15]  =                             - p[14]*p[15] - p[20]*p[21] - p[26]*p[27] - p[32]*p[33];
        ecLc[16]  =                                                         - p[26]*p[28] - p[32]*p[34];
        ecLc[17]  =                                                                       - p[32]*p[35];
                   
        ecLc[21] += - p[3] *p[3]  - p[9] *p[9]  - p[15]*p[15] - p[21]*p[21] - p[27]*p[27] - p[33]*p[33]; 
        ecLc[22] +=                                                         - p[27]*p[28] - p[33]*p[34]; 
        ecLc[23] +=                                                                       - p[33]*p[35]; 
                   
        ecLc[28] +=                                                         - p[28]*p[28] - p[34]*p[34]; 
        ecLc[29] +=                                                                       - p[34]*p[35]; 
                   
        ecLc[35] +=                                                                       - p[35]*p[35];  


        // chol (L(k+1,k+1))
        chol_dec (ecLc);
    }



    /**
     * @brief Builds matrix L.
     *
     * @param[in] ppar      parameters.
     * @param[in] i2hess    2*N diagonal elements of inverted hessian.
     */
    void matrix_ecL::form (const problem_parameters& ppar, const double *i2hess)
    {
        int i;
        state_parameters stp;

        stp = ppar.spar[0];

        // the first matrix on diagonal
        form_M (stp.sin, stp.cos, ppar.i2Q, i2hess);
        form_MBiPB (stp.B, ppar.i2P, ecL);
        form_L_diag (ecL);

        // offsets
        double *ecL_cur = &ecL[MATRIX_SIZE_6x6];
        double *ecL_prev = &ecL[0];
        for (i = 1; i < ppar.N; i++)
        {
            stp = ppar.spar[i];

            // form all matrices
            form_MAT (stp.A3, stp.A6);
            form_L_non_diag (ecL_prev, ecL_cur);

            // update offsets
            ecL_cur = &ecL_cur[MATRIX_SIZE_6x6];
            ecL_prev = &ecL_prev[MATRIX_SIZE_6x6];


            i2hess = &i2hess[2];
            form_M (stp.sin, stp.cos, ppar.i2Q, i2hess);
            form_MBiPB (stp.B, ppar.i2P, ecL_cur);
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
     *                  (N * #SMPC_NUM_STATE_VAR)
     */
    void matrix_ecL::solve_forward(const int N, double *x)
    {
        double *xc = x; // 6 current elements of x
        double *xp; // 6 elements of x computed on the previous iteration
        double *ecL_cur = &ecL[0];  // lower triangular matrix lying on the 
                                    // diagonal of L
        double *ecL_prev;   // upper triangular matrix lying to the left from
                            // ecL_cur at the same level of L


        // compute the first 6 elements using forward substitution
        xc[0] /= ecL_cur[0];
        xc[1] = (xc[1] - xc[0]*ecL_cur[1]) / ecL_cur[7];
        xc[2] = (xc[2] - xc[0]*ecL_cur[2] - xc[1]*ecL_cur[8]) / ecL_cur[14];
        xc[3] = (xc[3] - xc[0]*ecL_cur[3] - xc[1]*ecL_cur[9]  - xc[2]*ecL_cur[15]) / ecL_cur[21];
        xc[4] = (xc[4] - xc[0]*ecL_cur[4] - xc[1]*ecL_cur[10] - xc[2]*ecL_cur[16] - xc[3]*ecL_cur[22]) / ecL_cur[28];
        xc[5] = (xc[5] - xc[0]*ecL_cur[5] - xc[1]*ecL_cur[11] - xc[2]*ecL_cur[17] - xc[3]*ecL_cur[23] - xc[4]*ecL_cur[29]) / ecL_cur[35];


        for (int i = 1; i < N; i++)
        {
            // switch to the next level of L / next 6 elements
            xp = xc;
            xc = &xc[SMPC_NUM_STATE_VAR];

            ecL_prev = &ecL_cur[MATRIX_SIZE_6x6];
            ecL_cur = &ecL_prev[MATRIX_SIZE_6x6];

            
            // update the right part of the equation
            /*
             * xxxxxx
             *  xxxxx
             *   xxxx
             * xxxxxx
             *     xx
             *      x
             */
            xc[0] -= xp[0]*ecL_prev[0] + xp[1]*ecL_prev[6] + xp[2]*ecL_prev[12] + xp[3]*ecL_prev[18] + xp[4]*ecL_prev[24] + xp[5]*ecL_prev[30];
            xc[1] -=                     xp[1]*ecL_prev[7] + xp[2]*ecL_prev[13] + xp[3]*ecL_prev[19] + xp[4]*ecL_prev[25] + xp[5]*ecL_prev[31];
            xc[2] -=                                         xp[2]*ecL_prev[14] + xp[3]*ecL_prev[20] + xp[4]*ecL_prev[26] + xp[5]*ecL_prev[32];
            xc[3] -= xp[0]*ecL_prev[3] + xp[1]*ecL_prev[9] + xp[2]*ecL_prev[15] + xp[3]*ecL_prev[21] + xp[4]*ecL_prev[27] + xp[5]*ecL_prev[33];
            xc[4] -=                                                                                   xp[4]*ecL_prev[28] + xp[5]*ecL_prev[34];
            xc[5] -=                                                                                                        xp[5]*ecL_prev[35];

            // forward substitution
            xc[0] /= ecL_cur[0];
            xc[1] = (xc[1] - xc[0]*ecL_cur[1]) / ecL_cur[7];
            xc[2] = (xc[2] - xc[0]*ecL_cur[2] - xc[1]*ecL_cur[8]) / ecL_cur[14];
            xc[3] = (xc[3] - xc[0]*ecL_cur[3] - xc[1]*ecL_cur[9]  - xc[2]*ecL_cur[15]) / ecL_cur[21];
            xc[4] = (xc[4] - xc[0]*ecL_cur[4] - xc[1]*ecL_cur[10] - xc[2]*ecL_cur[16] - xc[3]*ecL_cur[22]) / ecL_cur[28];
            xc[5] = (xc[5] - xc[0]*ecL_cur[5] - xc[1]*ecL_cur[11] - xc[2]*ecL_cur[17] - xc[3]*ecL_cur[23] - xc[4]*ecL_cur[29]) / ecL_cur[35];
        }
    }


    /**
     * @brief Solve system ecL' * x = b using backward substitution.
     *
     * @param[in] N number of states in the preview window
     * @param[in,out] x vector "b" as input, vector "x" as output.
     */
    void matrix_ecL::solve_backward (const int N, double *x)
    {
        double *xc = & x[(N-1)*SMPC_NUM_STATE_VAR]; // current 6 elements of result
        double *xp; // 6 elements computed on the previous iteration
        
        // elements of these matrices accessed as if they were transposed
        // lower triangular matrix lying on the diagonal of L
        double *ecL_cur = &ecL[2 * (N - 1) * MATRIX_SIZE_6x6];
        // upper triangular matrix lying to the right from ecL_cur at the same level of L'
        double *ecL_prev; 


        // compute the last 6 elements using backward substitution
        xc[5] /= ecL_cur[35];
        xc[4] = (xc[4] - xc[5]*ecL_cur[29]) / ecL_cur[28];
        xc[3] = (xc[3] - xc[5]*ecL_cur[23] - xc[4]*ecL_cur[22]) / ecL_cur[21];
        xc[2] = (xc[2] - xc[5]*ecL_cur[17] - xc[4]*ecL_cur[16] - xc[3]*ecL_cur[15]) / ecL_cur[14];
        xc[1] = (xc[1] - xc[5]*ecL_cur[11] - xc[4]*ecL_cur[10] - xc[3]*ecL_cur[9] - xc[2]*ecL_cur[8]) / ecL_cur[7];
        xc[0] = (xc[0] - xc[5]*ecL_cur[5]  - xc[4]*ecL_cur[4]  - xc[3]*ecL_cur[3] - xc[2]*ecL_cur[2] - xc[1]*ecL_cur[1]) / ecL_cur[0];

        for (int i = N-2; i >= 0 ; i--)
        {
            xp = xc;
            xc = & x[i*SMPC_NUM_STATE_VAR];

            ecL_cur = &ecL[2 * i * MATRIX_SIZE_6x6];
            ecL_prev = &ecL_cur[MATRIX_SIZE_6x6];


            // update the right part of the equation
            /*
             * x  x
             * xx x
             * xxxx
             * xxxx
             * xxxxx
             * xxxxxx
             */
            xc[0] -= xp[0]*ecL_prev[0]                                            + xp[3]*ecL_prev[3];
            xc[1] -= xp[0]*ecL_prev[6]  + xp[1]*ecL_prev[7]                       + xp[3]*ecL_prev[9]; 
            xc[2] -= xp[0]*ecL_prev[12] + xp[1]*ecL_prev[13] + xp[2]*ecL_prev[14] + xp[3]*ecL_prev[15];
            xc[3] -= xp[0]*ecL_prev[18] + xp[1]*ecL_prev[19] + xp[2]*ecL_prev[20] + xp[3]*ecL_prev[21];
            xc[4] -= xp[0]*ecL_prev[24] + xp[1]*ecL_prev[25] + xp[2]*ecL_prev[26] + xp[3]*ecL_prev[27] + xp[4]*ecL_prev[28];
            xc[5] -= xp[0]*ecL_prev[30] + xp[1]*ecL_prev[31] + xp[2]*ecL_prev[32] + xp[3]*ecL_prev[33] + xp[4]*ecL_prev[34] + xp[5]*ecL_prev[35];

            // backward substitution
            xc[5] /= ecL_cur[35];
            xc[4] = (xc[4] - xc[5]*ecL_cur[29]) / ecL_cur[28];
            xc[3] = (xc[3] - xc[5]*ecL_cur[23] - xc[4]*ecL_cur[22]) / ecL_cur[21];
            xc[2] = (xc[2] - xc[5]*ecL_cur[17] - xc[4]*ecL_cur[16] - xc[3]*ecL_cur[15]) / ecL_cur[14];
            xc[1] = (xc[1] - xc[5]*ecL_cur[11] - xc[4]*ecL_cur[10] - xc[3]*ecL_cur[9] - xc[2]*ecL_cur[8]) / ecL_cur[7];
            xc[0] = (xc[0] - xc[5]*ecL_cur[5]  - xc[4]*ecL_cur[4]  - xc[3]*ecL_cur[3] - xc[2]*ecL_cur[2] - xc[1]*ecL_cur[1]) / ecL_cur[0];
        }
    }
}
