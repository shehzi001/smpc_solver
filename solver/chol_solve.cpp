/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
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

/**
 * @brief Constructor
 *
 * @param[in] preview_win_size size of the preview window.
 */
chol_solve::chol_solve (int preview_win_size)
{
    N = preview_win_size;

    ecL = new double[MATRIX_SIZE*N + MATRIX_SIZE*(N-1)]();

    nu = new double[NUM_VAR*N];
    XiHg = new double[NUM_VAR*N];

    z = new double[NUM_VAR*N];


    icL = new double*[N*2];
    icL_mem = new double[NUM_VAR*N*N*2];
    for(int i = 0; i < N*2; i++)
    {
        icL[i] = &icL_mem[i * NUM_VAR*N];
    }
}


chol_solve::~chol_solve()
{
    if (ecL != NULL)
        delete ecL;
    if (nu != NULL)
        delete nu;
    if (XiHg != NULL)
        delete XiHg;
    if (icL != NULL)
    {
        delete icL_mem;
        delete icL;
    }
    if (z != NULL)
        delete z;
}
//==============================================



/**
 * @brief Forms E*x.
 *
 * @param[in] csp parameters.
 * @param[in] x vector x (#NUM_VAR * N).
 * @param[out] result vector E*x (#NUM_STATE_VAR * N)
 */
void chol_solve::form_Ex (chol_solve_param csp, double *x, double *result)
{
    int i;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef QPAS_VARIABLE_T_h
    double T = csp.T[0];
    double T2 = T*T/2;
    double B0 = T2*T/3 - csp.h[0]*T;
    double A6 = T2;
#endif


    double *control = &x[N*NUM_STATE_VAR];
    double *res = result;

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

        // a pointer to 6 current state variables
        double *xc = &x[i*NUM_STATE_VAR];



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

        // next control variables
        control = &control[NUM_CONTROL_VAR];
        res = &res[NUM_STATE_VAR];
    }
}


/**
 * @brief Forms E' * x
 *
 * @param[in] csp parameters.
 * @param[in] x vector x (#NUM_STATE_VAR * N).
 * @param[out] result vector E' * nu (#NUM_VAR * N)
 */
void chol_solve::form_ETx (chol_solve_param csp, double *x, double *result)
{
    int i;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef QPAS_VARIABLE_T_h
    double T = csp.T[0];
    double T2 = T*T/2;
    double B0 = T2*T/3 - csp.h[0]*T;
    double A3 = T;
    double A6 = T2;
#endif


    double *res = result;
    double *control_res = &result[N*NUM_STATE_VAR];

    for (i = 0; i < N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = csp.angle_cos[i];
        double sinA = csp.angle_sin[i];


        // a pointer to 6 current elements of result
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


        xc = &x[i*NUM_STATE_VAR];
#ifdef QPAS_VARIABLE_T_h
        double T = csp.T[i];
        double T2 = T*T/2;
        double B0 = T2*T/3 - csp.h[i]*T;
#endif

        // result = B' * x
        control_res[0] = B0 * xc[0] + T2 * xc[1] + T * xc[2];
        control_res[1] = B0 * xc[3] + T2 * xc[4] + T * xc[5];


        res = &res[NUM_STATE_VAR];
        control_res = &control_res[NUM_CONTROL_VAR];
    }
}



/**
 * @brief Solve system ecL * x = b using forward substitution.
 *
 * @param[in,out] x vector "b" as input, vector "x" as output
 *                  (N * #NUM_STATE_VAR)
 */
void chol_solve::solve_forward(double *x)
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

        prev_ecL = &cur_ecL[MATRIX_SIZE];
        cur_ecL = &cur_ecL[2 * MATRIX_SIZE];


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
    double *xp; // 6 elements computed on the previous iteration

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
        xp = xc;
        xc = & x[i*NUM_STATE_VAR];

        cur_ecL = &ecL[2 * i * MATRIX_SIZE];
        prev_ecL = &cur_ecL[MATRIX_SIZE];


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
 * @brief Forms row vector 's_a' (@ref pCholUp).
 *
 * @param[in] csp parameters
 * @param[in] ic_num number of constraint, for example 5 if 4 are already added 
 * @param[in] var_num number of constrained variable
 * @param[out] row 's_a' row
 */
void chol_solve::form_sa_row(chol_solve_param csp, int ic_num, int var_num, double *row)
{
    double aiH = csp.i2Q[0]; // a'*inv(H) = a'*inv(H)*a
    int state_num = var_num / 2; // number of state in the preview window
    int first_num = state_num * NUM_STATE_VAR; // first !=0 element
    double aiHcosA;
    double aiHsinA;


    // reset memory
    memset(row, 0, (NUM_STATE_VAR * N + ic_num) * sizeof(double));

    aiHcosA = aiH * csp.angle_cos[state_num];
    aiHsinA = aiH * csp.angle_sin[state_num];

    // compute elements of 'a'
    if (var_num%2 == 0)   // constraint on z_x
    {
        // a * -R
        row[first_num] = -aiHcosA;
        row[first_num + 3] = -aiHsinA;

        // a * A'*R'
        if (state_num != N-1)
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
        if (state_num != N-1)
        {
            row[first_num + 6] = -aiHsinA;
            row[first_num + 9] = aiHcosA;
        }
    }

    // initialize the last element in the row
    row[ic_num + N*NUM_STATE_VAR] = aiH;
}



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
    L_init.form_L(csp, N, ecL);

    // -(V + inv(H) * g)
    //  V - initial feasible point
    for (i = 0; i < NUM_VAR * N; i++)
    {
        XiHg[i] = -(ix[i] + csp.iHg[i]);
    }

    // obtain s = E * x;
    form_Ex (csp, XiHg, s_nu);

    // obtain nu
    solve_forward(s_nu);
    // make copy of z - it is constant
    for (i = 0; i < NUM_STATE_VAR * N; i++)
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
    // dx   -(   -XiHg       + inv(H) *   dx   ) 
    for (i = 0; i < N*NUM_STATE_VAR; i++)
    {
        // dx for state variables
        dx[i] = -(-XiHg[i] + i2Q[i%3] * dx[i]);
    }
    for (; i < N*NUM_VAR; i++)
    {
        // dx for control variables
        dx[i] = -(-XiHg[i] + csp.i2P * dx[i]);
    }
}


/**
 * @brief A wrapper around private functions, which update Cholesky factor and 
 *  resolve the system.
 *
 * @param[in] csp   parameters.
 * @param[in] nW    number of added constrains.
 * @param[in] W     indicies of added constraints.
 * @param[in] x     initial guess.
 * @param[out] dx   feasible descent direction, must be allocated.
 */
void chol_solve::up_resolve(chol_solve_param csp, int nW, int *W, double *x, double *dx)
{
    update (csp, nW, W);
    update_z (csp, nW, W, x);
    resolve (csp, nW, W, x, dx);
}


/**
 * @brief Adds a row corresponding to some inequality constraint to L, see
 * '@ref pCholUpAlg'.
 *
 * @param[in] csp parameters.
 * @param[in] nW number of added inequality constraints + 1.
 * @param[in] W indexes of added inequality constraints + one index to be added.
 */
void chol_solve::update (chol_solve_param csp, int nW, int *W)
{
    int i, j, k;

    int ic_num = nW-1; // index of added constraint in W
    int state_num = W[ic_num] / 2; // number of state in the preview window
    double *new_row = icL[ic_num]; // current row in icL

    int last_num = ic_num + N*NUM_STATE_VAR; // the last !=0 element

    // a matrix on diagonal of ecL
    double* ecL_diag = &ecL[state_num * MATRIX_SIZE * 2];
    // a matrix below the ecL_diag
    double* ecL_ndiag = &ecL[state_num * MATRIX_SIZE * 2 + MATRIX_SIZE];


    // form row 'a' in the current row of icL
    form_sa_row(csp, ic_num, W[ic_num], new_row);

    // update elements starting from the first non-zero
    // element in the row to NUM_STATE_VAR * N (size of ecL)
    // the calculation of the last elements is completed
    // in a separate loop
    // each number in row 'a' causes update of only 3 elements following
    // it, they can be 1,2,6; 1,5,6; 4,5,6
    k = state_num*NUM_STATE_VAR;
    double* cur_pos = &new_row[k];
    for(i = state_num*2; i < N*2; i++) // variables corresponding to x and 
    {                                  // y are computed using the same matrices
        double tmp_copy_el[3];

        // ----------------------------------------------------------------
        cur_pos[0] /= ecL_diag[0];
        tmp_copy_el[0] = cur_pos[0];
        cur_pos[1] -= tmp_copy_el[0] * ecL_diag[1];
        cur_pos[1] /= ecL_diag[4];
        tmp_copy_el[1] = cur_pos[1];
        cur_pos[2] -= tmp_copy_el[0] * ecL_diag[2] + tmp_copy_el[1] * ecL_diag[5];
        

        // ----------------------------------------------------------------
        cur_pos[2] /= ecL_diag[8];
        tmp_copy_el[2] = cur_pos[2];
        if (i < 2*(N - 1))  // non-diagonal matrix of ecL cannot be
        {                   // used when the last state is processed

            // these elements can be updated here, since they are not 
            // used in computation of other elements on this iteration
            cur_pos[6] -= tmp_copy_el[0] * ecL_ndiag[0] 
                        + tmp_copy_el[1] * ecL_ndiag[3]
                        + tmp_copy_el[2] * ecL_ndiag[6];

            cur_pos[7] -= tmp_copy_el[1] * ecL_ndiag[4]
                        + tmp_copy_el[2] * ecL_ndiag[7];

            cur_pos[8] -= tmp_copy_el[2] * ecL_ndiag[8];

            if ((i+1)%2 == 0) // jump to the next pair of matrices in ecL.
            {
                ecL_diag = &ecL_diag[MATRIX_SIZE * 2];
                ecL_ndiag = &ecL_ndiag[MATRIX_SIZE * 2];
            }
        }
        // update the last (diagonal) number in the row
        new_row[last_num] -= tmp_copy_el[0] * tmp_copy_el[0] 
                           + tmp_copy_el[1] * tmp_copy_el[1] 
                           + tmp_copy_el[2] * tmp_copy_el[2];

        // update elements after N*NUM_STATE_VAR using the previously added rows
        // in icL
        for (j = 0; j < ic_num; j++)
        {
            new_row[N*NUM_STATE_VAR + j] -= tmp_copy_el[0] * icL[j][k]
                                          + tmp_copy_el[1] * icL[j][k+1]
                                          + tmp_copy_el[2] * icL[j][k+2];
        }
        cur_pos = &cur_pos[3];
        k += 3;
    }

    // update elements in the end of icL
    for(i = NUM_STATE_VAR * N; i < last_num; i++)
    {
        new_row[i] /= icL[i - N*NUM_STATE_VAR][i];
        double tmp_copy_el = new_row[i];

        // determine number in the row of L

        // update the last (diagonal) number in the row
        new_row[last_num] -= tmp_copy_el * tmp_copy_el;

        for (j = (i - N*NUM_STATE_VAR) + 1; j < ic_num; j++)
        {
            new_row[N*NUM_STATE_VAR + j] -= tmp_copy_el * icL[j][i];
        }
    }

    // square root of the diagonal element
    new_row[last_num] = sqrt(new_row[last_num]);
}



/**
 * @brief Adjust vector '@ref pz "z"' after update.
 *
 * @param[in] csp   parameters.
 * @param[in] nW    number of added constrains.
 * @param[in] W     indicies of added constraints.
 * @param[in] x     initial guess.
 */
void chol_solve::update_z (chol_solve_param csp, int nW, int *W, double *x)
{
    int ic_num = nW-1; // index of added constraint in W
    // update lagrange multipliers
    int zind = N*NUM_STATE_VAR + ic_num;
    // sn
    double zn = -(csp.iHg[W[ic_num]*3] + x[W[ic_num]*3]);

    // zn
    for (int i = 0; i < zind; i++)
    {
        zn -= z[i] * icL[ic_num][i];
        nu[i] = z[i];
    }
    nu[zind] = z[zind] = zn/icL[ic_num][zind];
    return;
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

    double i2Q[3] = {csp.i2Q[0], csp.i2Q[1], csp.i2Q[2]};


    // backward substituition for icL
    for (i = NUM_STATE_VAR*N + nW-1; i >= NUM_STATE_VAR*N; i--)
    {
        double nui = nu[i];
        double *icL_row = icL[i-N*NUM_STATE_VAR];

        nui = nui / icL_row[i];
        for (j = i - 1; j >= 0; j--)
        {
            nu[j] -= nui * icL_row[j];
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
    for (i = 0; i < N*NUM_STATE_VAR; i++)
    {
        // dx for state variables
        dx[i] = -(x[i] + csp.iHg[i] + i2Q[i%3] * dx[i]);
    }
    for (; i < N*NUM_VAR; i++)
    {
        // dx for control variables
        dx[i] = -(x[i] + csp.iHg[i] + csp.i2P * dx[i]);
    }

    // -iH * A(W,:)' * lambda
    double *lambda = &nu[N*NUM_STATE_VAR];
    for (i = 0; i < nW; i++)
    {
        dx[W[i]*3] -= i2Q[0] * lambda[i];
    }
}



#ifdef QPAS_DOWNDATE
/**
 * @brief A wrapper around private functions, which downdate Cholesky factor and 
 *  resolve the system.
 *
 * @param[in] csp   parameters.
 * @param[in] nW    number of added constrains (without removed constraint).
 * @param[in] W     indicies of added constraints (without removed constraint).
 * @param[in] ind_exclude index of excluded constraint.
 * @param[in] x     initial guess.
 * @param[out] dx   feasible descent direction, must be allocated.
 *
 * @note Downdate of vector @ref pz 'z' is described on the page '@ref pRemoveICz'.
 */
void chol_solve::down_resolve(chol_solve_param csp, int nW, int *W, int ind_exclude, double *x, double *dx)
{
    // for each element of z affected by removed constraint
    // find a base that stays the same
    double z_tmp = 0;
    for (int i = nW; i > ind_exclude; i--)
    {
        int zind = N*NUM_STATE_VAR + i;
        double zn = z[zind] * icL[i][zind];
        z[zind] = z_tmp;

        for (int j = N*NUM_STATE_VAR + ind_exclude; j < zind; j++)
        {
            zn += z[j] * icL[i][j];
        }
        z_tmp = zn;
    }
    z[N*NUM_STATE_VAR + ind_exclude] = z_tmp;


    // downdate L
    downdate (csp, nW, ind_exclude, x);


    // recompute elements of z
    for (int i = ind_exclude; i < nW; i++)
    {
        int zind = N*NUM_STATE_VAR + i;
        double zn = z[zind];

        // zn
        // start from the first !=0 element
        for (int j = N*NUM_STATE_VAR+ind_exclude; j < zind; j++)
        {
            zn -= z[j] * icL[i][j];
        }
        z[zind] = zn/icL[i][zind];
    }

    // copy z to nu
    for (int i = 0; i < N*NUM_STATE_VAR + nW; i++)
    {
        nu[i] = z[i];
    }


    resolve (csp, nW, W, x, dx);
}



/**
 * @return a pointer to the memory where current lambdas are stored.
 */
double * chol_solve::get_lambda()
{
    return(&nu[NUM_STATE_VAR*N]);
}




/**
 * @brief Delete a line from icL, see page '@ref pCholDown'.
 *
 * @param[in] csp   parameters.
 * @param[in] nW    number of added constrains.
 * @param[in] ind_exclude index of excluded constraint.
 * @param[in] x     initial guess.
 */
void chol_solve::downdate(chol_solve_param csp, int nW, int ind_exclude, double *x)
{
    // Shuffle memory pointers to avoid copying of the data.
    double * downdate_row = icL[ind_exclude];
    for (int i = ind_exclude + 1; i < nW + 1; i++)
    {
        icL[i-1] = icL[i];
    }
    icL[nW] = downdate_row;


    for (int i = ind_exclude; i < nW; i++)
    {
        int el_index = NUM_STATE_VAR*N + i;
        double *cur_el = &icL[i][el_index];
        double x1 = cur_el[0];
        double x2 = cur_el[1];
        double cosT, sinT;


        // Givens rotation matrix
        if (abs(x2) >= abs(x1))
        {
            double t = x1/x2;
            sinT = 1/sqrt(1 + t*t);
            cosT = sinT*t;
        }
        else
        {
            double t = x2/x1;
            cosT = 1/sqrt(1 + t*t);
            sinT = cosT*t;
        }


        // update elements in the current line
        cur_el[0] = cosT*x1 + sinT*x2;
        cur_el[1] = 0;

        // change sign if needed (diagonal elements of Cholesky 
        // decomposition must be positive)
        double sign = copysign(1, cur_el[0]);
        cur_el[0] = fabs(cur_el[0]);

        // update the lines below the current one.
        for (int j = i + 1; j < nW; j++)
        {
            x1 = icL[j][el_index];
            x2 = icL[j][el_index + 1];

            icL[j][el_index] = sign * (cosT*x1 + sinT*x2);
            icL[j][el_index + 1] = -sinT*x1 + cosT*x2;
        }
    }
}
#endif
