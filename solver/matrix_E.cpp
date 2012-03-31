/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "matrix_E.h"



/****************************************
 * FUNCTIONS 
 ****************************************/

/**
 * @brief Forms E*x.
 *
 * @param[in] ppar parameters.
 * @param[in] x vector x (#SMPC_NUM_VAR * N).
 * @param[out] result vector E*x (#SMPC_NUM_STATE_VAR * N)
 */
void matrix_E::form_Ex (const problem_parameters& ppar, const double *x, double *result)
{
    const double *control = &x[ppar.N*SMPC_NUM_STATE_VAR];
    // a pointer to 6 current elements of result
    double *res = result;
    // a pointer to 6 current state variables
    const double *xc = x;
    const double *xcp = NULL;


    state_parameters stp = ppar.spar[0];

    // result = -I * x + B * u
    res[0] = -xc[0] + stp.B[0] * control[0];
    res[1] = -xc[1] + stp.B[1] * control[0];
    res[2] = -xc[2] + stp.B[2] * control[0];
    res[3] = -xc[3] + stp.B[0] * control[1];
    res[4] = -xc[4] + stp.B[1] * control[1];
    res[5] = -xc[5] + stp.B[2] * control[1];


    for (int i = 1; i < ppar.N; i++)
    {
        // next control variables
        control = &control[SMPC_NUM_CONTROL_VAR];
        res = &res[SMPC_NUM_STATE_VAR];
        xcp = xc;
        xc = &xc[SMPC_NUM_STATE_VAR];

        stp = ppar.spar[i];

        // result = -I * x + B * u + A * xp
        res[0] = -xc[0] + stp.B[0] * control[0] + xcp[0] + stp.A3 * xcp[1] + stp.A6 * xcp[2];
        res[1] = -xc[1] + stp.B[1] * control[0] +                   xcp[1] + stp.A3 * xcp[2];
        res[2] = -xc[2] + stp.B[2] * control[0] +                                     xcp[2];
        res[3] = -xc[3] + stp.B[0] * control[1] + xcp[3] + stp.A3 * xcp[4] + stp.A6 * xcp[5];
        res[4] = -xc[4] + stp.B[1] * control[1] +                   xcp[4] + stp.A3 * xcp[5];
        res[5] = -xc[5] + stp.B[2] * control[1] +                                     xcp[5];
    }
}


/**
 * @brief Forms E' * x
 *
 * @param[in] ppar parameters.
 * @param[in] x vector x (#SMPC_NUM_STATE_VAR * N).
 * @param[out] result vector E' * nu (#SMPC_NUM_VAR * N)
 */
void matrix_E::form_ETx (const problem_parameters& ppar, const double *x, double *result)
{
    int i;

    // a pointer to 6 current elements of result
    double *res = result;
    double *control_res = &result[ppar.N*SMPC_NUM_STATE_VAR];
    // a pointer to 6 current elements of nu
    const double *xc = x;
    const double *xcn = &xc[SMPC_NUM_STATE_VAR];

    for (i = 0; i < ppar.N-1; i++)
    {
        double A3 = ppar.spar[i+1].A3;
        double A6 = ppar.spar[i+1].A6;


        // result = -I' * nu + A' * x
        res[0] = -xc[0] +      xcn[0];
        res[1] = -xc[1] + A3 * xcn[0] +      xcn[1];
        res[2] = -xc[2] + A6 * xcn[0] + A3 * xcn[1] + xcn[2];
        res[3] = -xc[3] +      xcn[3];                                                             
        res[4] = -xc[4] + A3 * xcn[3] +      xcn[4];
        res[5] = -xc[5] + A6 * xcn[3] + A3 * xcn[4] + xcn[5];
                       

        // result = B' * x
        state_parameters stp = ppar.spar[i];
        control_res[0] = stp.B[0] * xc[0] + stp.B[1] * xc[1] + stp.B[2] * xc[2];
        control_res[1] = stp.B[0] * xc[3] + stp.B[1] * xc[4] + stp.B[2] * xc[5];


        res = &res[SMPC_NUM_STATE_VAR];
        control_res = &control_res[SMPC_NUM_CONTROL_VAR];
        xc = xcn;
        xcn = &xcn[SMPC_NUM_STATE_VAR];
    }


    // result = -I' * nu
    res[0] = -xc[0];
    res[1] = -xc[1];
    res[2] = -xc[2];
    res[3] = -xc[3];
    res[4] = -xc[4];
    res[5] = -xc[5];


    // result = B' * x
    state_parameters stp = ppar.spar[i];
    control_res[0] = stp.B[0] * xc[0] + stp.B[1] * xc[1] + stp.B[2] * xc[2];
    control_res[1] = stp.B[0] * xc[3] + stp.B[1] * xc[4] + stp.B[2] * xc[5];
}
