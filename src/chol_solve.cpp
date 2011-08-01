/** 
 * @file
 * @brief  
 *
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 * @todo add description
 *  @todo steps
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "chol_solve.h"

#include <cmath> // sqrt
#include <cstring> // memset


/****************************************
 * FUNCTIONS 
 ****************************************/

//==============================================
// constructors / destructors

chol_solve::chol_solve (int steps)
{
    N = steps;

    ecL = new double[MATRIX_SIZE*N + MATRIX_SIZE*(N-1)]();
    iQBiPB = new double[MATRIX_SIZE];
    iQAT = new double[MATRIX_SIZE];
    AiQATiQBiPB = new double[MATRIX_SIZE];

    nu = new double[NUM_VAR*N];
    ViHg = new double[NUM_VAR*N];

    z = new double[NUM_VAR*N]();

    icL = new double*[N*2];
    icL[0] = new double[NUM_VAR*N*2*N];
    for(int i = 1; i < N*2; i++)
    {
        icL[i] = &icL[0][i * NUM_VAR*N];
    }
}


chol_solve::~chol_solve()
{
    if (ecL != NULL)
        delete ecL;
    if (iQBiPB != NULL)
        delete iQBiPB;
    if (iQAT != NULL)
        delete iQAT;
    if (AiQATiQBiPB != NULL)
        delete AiQATiQBiPB;
    if (nu != NULL)
        delete nu;
    if (ViHg != NULL)
        delete ViHg;
    if (icL != NULL)
    {
        delete icL[0];
        delete icL;
    }
    if (z != NULL)
        delete z;
}



/*********** private functions ************/

/**
 * @brief Forms E*x.
 *
 * @param[in] csp parameters.
 * @param[in] x vector x (NUM_VAR * N).
 * @param[out] result vector E*x (NUM_STATE_VAR * N)
 */
void chol_solve::form_Ex (chol_solve_param csp, double *x, double *result)
{
    int i;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef QPAS_VARIABLE_T_h
    double T = csp.T;
    double T2 = T*T/2;
    double B0 = T2*T/3 - csp.h*T;
    double A6 = T2;
#endif


    for (i = 0; i < N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = csp.angle_cos[i];
        double sinA = csp.angle_sin[i];
#ifdef QPAS_VARIABLE_T_h
        double T = csp.T[i];
        double T2 = T*T/2;
        double B0 = T2*T/3 - csp.h[i]*T;
#endif

        // a pointer to 6 current elements of result
        double *res = &result[i*NUM_STATE_VAR];

        // a pointer to 6 current state variables
        double *xc = &x[i*NUM_STATE_VAR];

        // two current control variables
        double control[2] = {
            x[N*NUM_STATE_VAR + i*NUM_CONTROL_VAR + 0], 
            x[N*NUM_STATE_VAR + i*NUM_CONTROL_VAR + 1]};


        // result = -R * x + B * u
        res[0] = -(cosA * xc[0] - sinA * xc[3]) + B0 * control[0];
        res[1] = -xc[1]                         + T2 * control[0];
        res[2] = -xc[2]                         +  T * control[0];
        res[3] = -(sinA * xc[0] + cosA * xc[3]) + B0 * control[1];
        res[4] = -xc[4]                         + T2 * control[1];
        res[5] = -xc[5]                         +  T * control[1];


        if (i != 0) // no multiplication by A on the first iteration
        {
            int j = i-1;
#ifdef QPAS_VARIABLE_T_h
            double A6 = T2 - csp.dh[j];
#endif
            xc = &x[j*NUM_STATE_VAR];

            cosA = csp.angle_cos[j];
            sinA = csp.angle_sin[j];

            // result += A*R*x
            res[0] += cosA * xc[0] + T * xc[1] + A6 * xc[2] - sinA * xc[3];
            res[1] +=                    xc[1] +  T * xc[2];
            res[2] +=                                 xc[2];
            res[3] += cosA * xc[3] + T * xc[4] + A6 * xc[5] + sinA * xc[0];
            res[4] +=                    xc[4] +  T * xc[5]; 
            res[5] +=                                 xc[5];
        }
    }
}


/**
 * @brief Forms E' * x
 *
 * @param[in] csp parameters.
 * @param[in] x vector x (NUM_STATE_VAR * N).
 * @param[out] result vector E' * nu (NUM_VAR * N)
 */
void chol_solve::form_ETx (chol_solve_param csp, double *x, double *result)
{
    int i;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef QPAS_VARIABLE_T_h
    double T = csp.T;
    double T2 = T*T/2;
    double B0 = T2*T/3 - csp.h*T;
    double A3 = T;
    double A6 = T2;
#endif

    for (i = 0; i < N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = csp.angle_cos[i];
        double sinA = csp.angle_sin[i];


        // a pointer to 6 current elements of result
        double *res = &result[i*NUM_STATE_VAR];
        // a pointer to 6 current elements of nu
        double *xc = &x[i*NUM_STATE_VAR];


        // result = -R' * nu
        res[0] = -(cosA * xc[0] + sinA * xc[3]);
        res[1] = -xc[1];
        res[2] = -xc[2];
        res[3] = -(- sinA * xc[0] + cosA * xc[3]);
        res[4] = -xc[4];
        res[5] = -xc[5];


        if (i != N-1) // no multiplication by A on the last iteration
        {
#ifdef QPAS_VARIABLE_T_h
            double A3 = csp.T[i+1];
            double A6 = A3*A3/2 - csp.dh[i];
#endif

            xc = &x[i*NUM_STATE_VAR + NUM_STATE_VAR];

            // result += R' * A' * x
            res[0] += cosA * xc[0] + sinA * xc[3];
            res[1] +=   A3 * xc[0] + xc[1];
            res[2] +=   A6 * xc[0] + A3 * xc[1] + xc[2];

            res[3] += - sinA * xc[0] + cosA * xc[3];
            res[4] +=     A3 * xc[3] + xc[4];
            res[5] +=     A6 * xc[3] + A3 * xc[4] + xc[5]; 
        }


        res = &result[i*NUM_CONTROL_VAR + N*NUM_STATE_VAR];
        xc = &x[i*NUM_STATE_VAR];
#ifdef QPAS_VARIABLE_T_h
        double T = csp.T[i];
        double T2 = T*T/2;
        double B0 = T2*T/3 - csp.h[i]*T;
#endif

        // result = B' * x
        res[0] = B0 * xc[0] + T2 * xc[1] + T * xc[2];
        res[1] = B0 * xc[3] + T2 * xc[4] + T * xc[5];
    }
}



/**
 * @brief Performs Cholesky decomposition of 3x3 matrix.
 *
 * @param[in,out] mx9 a pointer to matrix, the result is 
 *                    stored in the same place.
 */
void chol_solve::chol_dec (double *mx9)
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

    // These elements must be 0.
    mx9[3] = mx9[6] = mx9[7] = 0;
}



/**
 * @brief Forms matrix iQBiPB = 0.5 * inv(Q) + 0.5 * B * inv(P) * B.
 *
 * @param[in] B a vector of 3 elements.
 * @param[in] i2Q a vector of 3 elements, which contains
 *              diagonal elements of 0.5 * inv(Q).
 * @param[in] iP 0.5 * inv(P) (only one number)
 */
void chol_solve::form_iQBiPB (double *B, double *i2Q, double iP)
{
    // diagonal elements
    iQBiPB[0] = 0.5*iP * B[0]*B[0] + i2Q[0];
    iQBiPB[4] = 0.5*iP * B[1]*B[1] + i2Q[1];
    iQBiPB[8] = 0.5*iP * B[2]*B[2] + i2Q[2];

    // symmetric elements
    iQBiPB[1] = iQBiPB[3] = 0.5*iP * B[0]*B[1];
    iQBiPB[2] = iQBiPB[6] = 0.5*iP * B[0]*B[2];

    /// @todo There is no need to initialize all non-diagonal elements
    ///       of symmetric matrices.
    iQBiPB[5] = iQBiPB[7] = 0.5*iP * B[1]*B[2];
}



/**
 * @brief Forms matrix iQAT = 0.5 * inv (Q) * A'
 *
 * @param[in] T 4th and 7th elements of A.
 * @param[in] A6 6th element of A.
 * @param[in] i2Q a vector of 3 elements, which contains
 *              diagonal elements of 0.5*inv(Q).
 */
void chol_solve::form_iQAT (double T, double A6, double *i2Q)
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
 */
void chol_solve::form_AiQATiQBiPB (double T, double A6)
{
    // 1st column
    AiQATiQBiPB[0] = iQBiPB[0] + iQAT[0] + T*iQAT[1] + A6*iQAT[2];
    AiQATiQBiPB[1] = iQBiPB[1] +             iQAT[1] +  T*iQAT[2];
    AiQATiQBiPB[2] = iQBiPB[2] +                          iQAT[2];

    // 2nd column
    AiQATiQBiPB[3] = iQBiPB[3] + T*iQAT[4] + A6*iQAT[5];
    AiQATiQBiPB[4] = iQBiPB[4] +   iQAT[4] +  T*iQAT[5];
    AiQATiQBiPB[5] = iQBiPB[5] +                iQAT[5];

    // 3rd column
    AiQATiQBiPB[6] = iQBiPB[6] + A6*iQAT[8];
    AiQATiQBiPB[7] = iQBiPB[7] +  T*iQAT[8];
    AiQATiQBiPB[8] = iQBiPB[8] +    iQAT[8];
}



/**
 * @brief Forms a 3x3 matrix L(k+1, k), which lies below the 
 *  diagonal of L.
 *
 * @param[in] ecLp previous matrix lying on the diagonal of L
 * @param[in] ecLc the result is stored here
 */
void chol_solve::form_L_non_diag(double *ecLp, double *ecLc)
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
 */
void chol_solve::form_L_diag(double *ecLp, double *ecLc)
{
    int i;


    // the first matrix L(0,0) is computed differently
    if (ecLp == NULL)
    {
        //0.5*inv(Q) + 0.5*B*inv(P)*B'
        for (i = 0; i < MATRIX_SIZE; i++)
        {
            ecLc[i] = iQBiPB[i];
        }
    }
    else
    {
    // L(k+1,k+1) = - L(k+1,k) * L(k+1,k)' + ...
        // diagonal elements
        ecLc[0] = -(ecLp[0]*ecLp[0] + ecLp[3]*ecLp[3] + ecLp[6]*ecLp[6]);
        ecLc[4] = -(ecLp[4]*ecLp[4] + ecLp[7]*ecLp[7]);
        ecLc[8] = -(ecLp[8]*ecLp[8]);
        // symmetric nondiagonal elements;
        ecLc[1] = ecLc[3] = -(ecLp[3]*ecLp[4] + ecLp[6]*ecLp[7]);
        ecLc[2] = ecLc[6] = -ecLp[6]*ecLp[8];
        ecLc[5] = ecLc[7] = -ecLp[7]*ecLp[8];

    // ... + A * inv(Q) * A' + inv(Q) + B * inv(P) * B
        for (i = 0; i < MATRIX_SIZE; i++)
        {
            ecLc[i] += AiQATiQBiPB[i];
        }
    }

    // chol (L(k+1,k+1))
    chol_dec (ecLc);
}



/**
 * @brief Builds matrix L.
 *
 * @param[in] csp parameters.
 */
void chol_solve::form_L(chol_solve_param csp)
{
    int i;
    int cur_offset;
    int prev_offset;

#ifdef QPAS_VARIABLE_T_h
    double T = csp.T[0];
    double T2 = T*T/2;
    double B[3] = {T2*T/3 - csp.h[0]*T, T2, T};
#else 
    double T = csp.T;
    double T2 = T*T/2;
    double B[3] = {T2*T/3 - csp.h*T, T2, T};
    double A6 = T2;
#endif

    // form all matrices
    form_iQBiPB (B, csp.i2Q, csp.iP);
#ifndef QPAS_VARIABLE_T_h
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
        form_iQBiPB (B, csp.i2Q, csp.iP);
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


/**
 * @brief Solve system ecL * x = b using forward substitution.
 *
 * @param[in,out] x vector "b" as input, vector "x" as output
 *                  (N * NUM_STATE_VAR)
 */
void chol_solve::solve_forward(double *x)
{
    int i;
    double *xc = &x[0]; // 6 current elements of x
    double *xp; // 6 elements of x computed on the previous step
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
        xc = &x[i * NUM_STATE_VAR];
        xp = &x[i * NUM_STATE_VAR - NUM_STATE_VAR];
        cur_ecL = &ecL[i * 2 * MATRIX_SIZE];
        prev_ecL = &ecL[i * 2 * MATRIX_SIZE - MATRIX_SIZE];


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
 * @param[in,out] x vector "b" as input, vector "x" as output.
 */
void chol_solve::solve_backward(double *x)
{
    int i;
    double *xc = & x[(N-1)*NUM_STATE_VAR]; // current 6 elements of result
    double *xp; // 6 elements computed on the previous step

    // elements of these matrices accessed as if they were transposed
    // lower triangular matrix lying on the diagonal of L
    double *cur_ecL = &ecL[2 * (N - 1) * MATRIX_SIZE]; 
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
        xc = & x[i*NUM_STATE_VAR];
        xp = & x[i*NUM_STATE_VAR + NUM_STATE_VAR];

        cur_ecL = &ecL[2 * i * MATRIX_SIZE];
        prev_ecL = &ecL[2 * i * MATRIX_SIZE + MATRIX_SIZE];


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


/**
 * @brief Forms row vector 'a'.

    Consider the following matrix (with "H" being a diagonal matrix, and "a"
    a column vector)
    
    \verbatim
    [C;a']*inv(H)*[C',a] = [C*inv(H)*C' , C*inv(H)*a;
                            a'*inv(H)*C', a'*inv(H)*a]

    Then,

    aRowVector = [a'*inv(H)*C', a'*inv(H)*a]

    Or, if C = [E;Aw] (with Aw containing the normals to the active inequality 
    constraints prior to the inclusion of "a")

    aRowVector = [a'*inv(H)*E', a'*inv(H)*Aw', a'*inv(H)*a])

    Note that a'*inv(H)*Aw' is a vector of zeros (with dimension 
    [1 x (number of active constraints)]). While a'*inv(H)*a = 1\Beta
                   -----------------------------------------------------
    aRowVector =   |      a'*inv(H)*E'      |  0 0 ... 0 | a'*inv(H)*a |
                   -----------------------------------------------------
    (dimensions)~~    6*N                      nW            1

    a'*inv(H)*E' = (E*inv(H)*a)' = (1/Beta*E*a)' (because "a" is avector 
    of zeros, with only one of 1:3:6*N equal to 1).
    Hence, a'*inv(H)*E' selects and scales (by 1/Beta) one column of E 
    (of course a column with index in 1:3:6*N).
    Such a column can have at most 4 nonzero entries.

    The row has 5 or 3 (for the last step) non-zero elements.
    \endverbatim
 *
 * @param[in] csp parameters
 * @param[in] ic_num number of constraint, for example 5 if 4 are already added 
 * @param[in] var_num number of constrained variable
 * @param[out] row 'a' row
 */
void chol_solve::form_a_row(chol_solve_param csp, int ic_num, int var_num, double *row)
{
    double aiH = csp.i2Q[0]; // a'*inv(H) = a'*inv(H)*a
    int step_num = var_num / 2; // step number
    int first_num = step_num * NUM_STATE_VAR; // first !=0 element
    double aiHcosA;
    double aiHsinA;


    // reset memory
    memset(row, 0, NUM_VAR * N *sizeof(double));

    aiHcosA = aiH * csp.angle_cos[step_num];
    aiHsinA = aiH * csp.angle_sin[step_num];

    // compute elements of 'a'
    if (var_num%2 == 0)   // constraint on z_x
    {
        // a * -R
        row[first_num] = -aiHcosA;
        row[first_num + 3] = -aiHsinA;

        // a * A'*R'
        if (step_num != N-1)
        {
            row[first_num + 6] = aiHcosA;
            row[first_num + 9] = aiHsinA;
        }
    }
    else  // constraint on z_y
    {   
        // a * -R
        row[first_num] = aiHsinA;
        row[first_num + 3] = -aiHcosA;

        // a * A'*R'
        if (step_num != N-1)
        {
            row[first_num + 6] = -aiHsinA;
            row[first_num + 9] = aiHcosA;
        }
    }

    // initialize the last element in the row
    row[ic_num + N*NUM_STATE_VAR] = aiH;
}



/*********** public functions ************/

/**
 * @brief Determines feasible descent direction.
 *
 * @param[in] csp   parameters.
 * @param[in] ix    initial guess.
 * @param[out] dx   feasible descent direction, must be allocated.
 */
void chol_solve::solve(chol_solve_param csp, double *ix, double *dx)
{
    double *s_nu = nu;
    int i;
    double i2Q[3] = {csp.i2Q[0], csp.i2Q[1], csp.i2Q[2]};


    // generate L
    form_L(csp);

    // -(V + inv(H) * g)
    //  V - initial feasible point
    for (i = 0; i < NUM_VAR * N; i++)
    {
        ViHg[i] = -(ix[i] + csp.iHg[i]);
    }

    // obtain s = E * x;
    form_Ex (csp, ViHg, s_nu);

    // obtain nu
    solve_forward(s_nu);
    // make copy of z - it is constant
    for (i= 0; i < NUM_STATE_VAR * N; i++)
    {
        z[i] = s_nu[i];
    }
    solve_backward(s_nu);

    // E' * nu
    form_ETx (csp, s_nu, dx);

    
    // dx = -iH*(grad + E'*nu)
    //
    // dx = -(V + inv(H) * g + inv(H) * E' * nu)
    //        ~~~~~~~~~~~~~~            ~~~~~~~
    // dx   -(   -ViHg       + inv(H) *   dx   ) 
    for (i = 0; i < N*NUM_VAR; i++)
    {
        if (i < N*NUM_STATE_VAR)
        {
            // dx for state variables
            dx[i] = -(-ViHg[i] + i2Q[i%3] * dx[i]);
        }
        else
        {
            // dx for control variables
            dx[i] = -(-ViHg[i] + csp.iP * dx[i]);
        }
    }
}



/**
 * @brief Adds a row corresponding to some inequality constraint to L.
 *
 * @param[in] csp parameters.
 * @param[in] nW number of added inequality constraints + 1.
 * @param[in] W indexes of added inequality constraints + one index to be added.
 */
void chol_solve::add_L_row (chol_solve_param csp, int nW, int *W)
{
    int i, j;

    int ic_num = nW-1; // index of added constraint in W
    int step_num = W[ic_num] / 2; // step number
    double cur_el; // temporary storage
    double *new_row = icL[ic_num]; // current row in icL

    int first_num = step_num * NUM_STATE_VAR; // the first !=0 element
    int last_num = ic_num + N*NUM_STATE_VAR; // the last !=0 element

    // a matrix on diagonal of ecL
    double* ecL_diag = &ecL[step_num * MATRIX_SIZE * 2];
    // a matrix below the ecL_diag
    double* ecL_ndiag = &ecL[step_num * MATRIX_SIZE * 2 + MATRIX_SIZE];


    // form row 'a' in the current row of icL
    form_a_row(csp, ic_num, W[ic_num], new_row);

    // update elements starting from the first non-zero
    // element in the row to NUM_STATE_VAR * N (size of ecL)
    // the calculation of the last elements is completed
    // in a separate loop
    for(i = first_num; i < NUM_STATE_VAR * N; i++)
    {
        cur_el = new_row[i];

        // propagate update in the row
        // each number in row 'a' causes update of only 3 elements following
        // it, they can be 1,2,6; 1,5,6; 4,5,6
        switch (i % 3)  // variables corresponding to x and y are computed
        {               // using the same matrices
            case 0:
                // determine number in the row of L
                cur_el /= ecL_diag[0];
                new_row[i + 1] -= cur_el * ecL_diag[1];
                new_row[i + 2] -= cur_el * ecL_diag[2];
                if (step_num != N-1) // this is not needed in the end of ecL
                {
                    new_row[i + 6] -= cur_el * ecL_ndiag[0];
                }
                break;

            case 1:
                // determine number in the row of L
                cur_el /= ecL_diag[4];
                new_row[i + 1] -= cur_el * ecL_diag[5];
                if (step_num != N-1) // this is not needed in the end of ecL
                {
                    new_row[i + 5] -= cur_el * ecL_ndiag[3];
                    new_row[i + 6] -= cur_el * ecL_ndiag[4];
                }
                break;

            case 2:
                // determine number in the row of L
                cur_el /= ecL_diag[8];
                if (step_num != N-1) // this is not needed in the end of ecL
                {
                    new_row[i + 4] -= cur_el * ecL_ndiag[6];
                    new_row[i + 5] -= cur_el * ecL_ndiag[7];
                    new_row[i + 6] -= cur_el * ecL_ndiag[8];
                }
                break;
        }

        // update the last (diagonal) number in the row
        new_row[last_num] -= cur_el * cur_el;

        // update elements after N*NUM_STATE_VAR using the previously added rows
        // in icL
        for (j = 0; j < ic_num; j++)
        {
            new_row[N*NUM_STATE_VAR + j] -= cur_el * icL[j][i];
        }

        // jump to the next pair of matrices in ecL.
        if ((i+1)%6 == 0)
        {
            step_num++;
            ecL_diag = &ecL[step_num * MATRIX_SIZE * 2];
            ecL_ndiag = &ecL[step_num * MATRIX_SIZE * 2 + MATRIX_SIZE];
        }
        new_row[i] = cur_el;
    }


    // update elements in the end of icL
    for(i = NUM_STATE_VAR * N; i < last_num; i++)
    {
        cur_el = new_row[i];

        // determine number in the row of L
        cur_el /= icL[i - N*NUM_STATE_VAR][i];

        // update the last (diagonal) number in the row
        new_row[last_num] -= cur_el * cur_el;

        for (j = (i - N*NUM_STATE_VAR) + 1; j < ic_num; j++)
        {
            new_row[N*NUM_STATE_VAR + j] -= cur_el * icL[j][i];
        }

        new_row[i] = cur_el;
    }

    // square root of the diagonal element
    new_row[last_num] = sqrt(new_row[last_num]);
}


/**
 * @brief Determines feasible descent direction with respect to added
 *  inequality constraints.
 *
 * @param[in] csp   parameters.
 * @param[in] nW    number of added constrains.
 * @param[in] W     indicies of added constraints.
 * @param[in] x     initial guess.
 * @param[out] dx   feasible descent direction, must be allocated.
 */
void chol_solve::resolve (chol_solve_param csp, int nW, int *W, double *x, double *dx)
{
    int i,j;

    int ic_num = nW-1;
    int zind = N*NUM_STATE_VAR + ic_num;
    double i2Q[3] = {csp.i2Q[0], csp.i2Q[1], csp.i2Q[2]};


    // sn
    double zn = -(csp.iHg[W[ic_num]*3] + x[W[ic_num]*3]);

    // zn
    for (i = 0; i < zind; i++)
    {
        zn -= z[i] * icL[ic_num][i];
        nu[i] = z[i];
    }
    nu[zind] = z[zind] = zn/icL[ic_num][zind];


    // backward substituition for icL
    for (i = zind; i >= NUM_STATE_VAR*N; i--)
    {
        double nui = nu[i];

        nui = nui / icL[i-N*NUM_STATE_VAR][i];
        for (j = i - 1; j >= 0; j--)
        {
            nu[j] -= nui * icL[i-N*NUM_STATE_VAR][j];
        }
        nu[i] = nui;
    }
    // backward substituition for ecL
    solve_backward(nu);


    // E' * nu
    form_ETx (csp, nu, dx);

    // dx = -iH*(grad + E'*nu  + A(W,:)'*lambda)
    //
    // dx = -(V + inv(H) * g + inv(H) * E' * nu)
    //        ~~~~~~~~~~~~~~            ~~~~~~~
    // dx   -(x +  iHg       + inv(H) *   dx   ) 
    for (i = 0; i < N*NUM_VAR; i++)
    {
        if (i < N*NUM_STATE_VAR)
        {
            // dx for state variables
            dx[i] = -(x[i] + csp.iHg[i] + i2Q[i%3] * dx[i]);
        }
        else
        {
            // dx for control variables
            dx[i] = -(x[i] + csp.iHg[i] + csp.iP * dx[i]);
        }
    }

    // -iH * A(W,:)' * lambda
    for (i = 0; i < nW; i++)
    {
        dx[W[i]*3] -= i2Q[0] * nu[N*NUM_STATE_VAR + i];
    }
}
