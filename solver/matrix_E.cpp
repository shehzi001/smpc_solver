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
    int i;
    state_parameters stp;


    const double *control = &x[ppar.N*SMPC_NUM_STATE_VAR];
    // a pointer to 6 current elements of result
    double *res = result;

    for (i = 0; i < ppar.N; i++)
    {
        stp = ppar.spar[i];

        // a pointer to 6 current state variables
        const double *xc = &x[i*SMPC_NUM_STATE_VAR];


        // result = -R * x + B * u
        res[0] = -(stp.cos * xc[0] - stp.sin * xc[3]) + stp.B[0] * control[0];
        res[1] = -xc[1]                               + stp.B[1] * control[0];
        res[2] = -xc[2]                               + stp.B[2] * control[0];
        res[3] = -(stp.sin * xc[0] + stp.cos * xc[3]) + stp.B[0] * control[1];
        res[4] = -xc[4]                               + stp.B[1] * control[1];
        res[5] = -xc[5]                               + stp.B[2] * control[1];


        if (i != 0) // no multiplication by A on the first iteration
        {
            double cosA = ppar.spar[i-1].cos;
            double sinA = ppar.spar[i-1].sin;

            xc = &x[(i-1)*SMPC_NUM_STATE_VAR];


            // result += A*R*x
            res[0] += cosA * xc[0] + stp.A3 * xc[1] + stp.A6 * xc[2] - sinA * xc[3];
            res[1] +=                         xc[1] + stp.A3 * xc[2];
            res[2] +=                                          xc[2];
            res[3] += cosA * xc[3] + stp.A3 * xc[4] + stp.A6 * xc[5] + sinA * xc[0];
            res[4] +=                         xc[4] + stp.A3 * xc[5]; 
            res[5] +=                                          xc[5];
        }

        // next control variables
        control = &control[SMPC_NUM_CONTROL_VAR];
        res = &res[SMPC_NUM_STATE_VAR];
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
    state_parameters stp;



    double *res = result;
    double *control_res = &result[ppar.N*SMPC_NUM_STATE_VAR];

    for (i = 0; i < ppar.N; i++)
    {
        stp = ppar.spar[i];


        // a pointer to 6 current elements of result
        // a pointer to 6 current elements of nu
        const double *xc = &x[i*SMPC_NUM_STATE_VAR];


        // result = -R' * nu
        res[0] = -(stp.cos * xc[0] + stp.sin * xc[3]);
        res[1] = -xc[1];
        res[2] = -xc[2];
        res[3] = -(- stp.sin * xc[0] + stp.cos * xc[3]);
        res[4] = -xc[4];
        res[5] = -xc[5];


        if (i != ppar.N-1) // no multiplication by A on the last iteration
        {
            double A3 = ppar.spar[i+1].A3;
            double A6 = ppar.spar[i+1].A6;

            xc = &x[i*SMPC_NUM_STATE_VAR + SMPC_NUM_STATE_VAR];

            // result += R' * A' * x
            res[0] += stp.cos * xc[0] + stp.sin * xc[3];
            res[1] +=      A3 * xc[0] + xc[1];
            res[2] +=      A6 * xc[0] + A3 * xc[1] + xc[2];

            res[3] += - stp.sin * xc[0] + stp.cos * xc[3];
            res[4] +=        A3 * xc[3] + xc[4];
            res[5] +=        A6 * xc[3] + A3 * xc[4] + xc[5]; 
        }


        xc = &x[i*SMPC_NUM_STATE_VAR];

        // result = B' * x
        control_res[0] = stp.B[0] * xc[0] + stp.B[1] * xc[1] + stp.B[2] * xc[2];
        control_res[1] = stp.B[0] * xc[3] + stp.B[1] * xc[4] + stp.B[2] * xc[5];


        res = &res[SMPC_NUM_STATE_VAR];
        control_res = &control_res[SMPC_NUM_CONTROL_VAR];
    }
}
