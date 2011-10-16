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
 * @param[in] x vector x (#NUM_VAR * N).
 * @param[out] result vector E*x (#NUM_STATE_VAR * N)
 */
void matrix_E::form_Ex (const problem_parameters* ppar, const double *x, double *result)
{
    int i;
    double B[3];
    double A3,A6;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef SMPC_VARIABLE_T_h
    A3 = B[2] = ppar->T;
    B[1] = ppar->B[1];
    B[0] = ppar->B[0];
    A6 = ppar->A6;
#endif


    const double *control = &x[ppar->N*NUM_STATE_VAR];
    double *res = result;

    for (i = 0; i < ppar->N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = ppar->spar[i].cos;
        double sinA = ppar->spar[i].sin;
#ifdef SMPC_VARIABLE_T_h
        A3 = B[2] = ppar->spar[i].T;
        B[1] = ppar->spar[i].B[1];
        B[0] = ppar->spar[i].B[0];
#endif

        // a pointer to 6 current elements of result

        // a pointer to 6 current state variables
        const double *xc = &x[i*NUM_STATE_VAR];



        // result = -R * x + B * u
        res[0] = -(cosA * xc[0] - sinA * xc[3]) + B[0] * control[0];
        res[1] = -xc[1]                         + B[1] * control[0];
        res[2] = -xc[2]                         + B[2] * control[0];
        res[3] = -(sinA * xc[0] + cosA * xc[3]) + B[0] * control[1];
        res[4] = -xc[4]                         + B[1] * control[1];
        res[5] = -xc[5]                         + B[2] * control[1];


        if (i != 0) // no multiplication by A on the first iteration
        {
            int j = i-1;
#ifdef SMPC_VARIABLE_T_h
            A6 = ppar->spar[i].A6;
#endif
            xc = &x[j*NUM_STATE_VAR];

            cosA = ppar->spar[j].cos;
            sinA = ppar->spar[j].sin;

            // result += A*R*x
            res[0] += cosA * xc[0] + A3 * xc[1] + A6 * xc[2] - sinA * xc[3];
            res[1] +=                     xc[1] + A3 * xc[2];
            res[2] +=                                  xc[2];
            res[3] += cosA * xc[3] + A3 * xc[4] + A6 * xc[5] + sinA * xc[0];
            res[4] +=                     xc[4] + A3 * xc[5]; 
            res[5] +=                                  xc[5];
        }

        // next control variables
        control = &control[NUM_CONTROL_VAR];
        res = &res[NUM_STATE_VAR];
    }
}


/**
 * @brief Forms E' * x
 *
 * @param[in] ppar parameters.
 * @param[in] x vector x (#NUM_STATE_VAR * N).
 * @param[out] result vector E' * nu (#NUM_VAR * N)
 */
void matrix_E::form_ETx (const problem_parameters* ppar, const double *x, double *result)
{
    int i;
    double B[3];
    double A3,A6;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef SMPC_VARIABLE_T_h
    A3 = B[2] = ppar->T;
    B[1] = ppar->B[1];
    B[0] = ppar->B[0];
    A6 = ppar->A6;
#endif


    double *res = result;
    double *control_res = &result[ppar->N*NUM_STATE_VAR];

    for (i = 0; i < ppar->N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = ppar->spar[i].cos;
        double sinA = ppar->spar[i].sin;


        // a pointer to 6 current elements of result
        // a pointer to 6 current elements of nu
        const double *xc = &x[i*NUM_STATE_VAR];


        // result = -R' * nu
        res[0] = -(cosA * xc[0] + sinA * xc[3]);
        res[1] = -xc[1];
        res[2] = -xc[2];
        res[3] = -(- sinA * xc[0] + cosA * xc[3]);
        res[4] = -xc[4];
        res[5] = -xc[5];


        if (i != ppar->N-1) // no multiplication by A on the last iteration
        {
#ifdef SMPC_VARIABLE_T_h
            A3 = ppar->spar[i+1].T;
            A6 = ppar->spar[i+1].A6;
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
#ifdef SMPC_VARIABLE_T_h
        B[2] = ppar->spar[i].T;
        B[1] = ppar->spar[i].B[1];
        B[0] = ppar->spar[i].B[0];
#endif

        // result = B' * x
        control_res[0] = B[0] * xc[0] + B[1] * xc[1] + B[2] * xc[2];
        control_res[1] = B[0] * xc[3] + B[1] * xc[4] + B[2] * xc[5];


        res = &res[NUM_STATE_VAR];
        control_res = &control_res[NUM_CONTROL_VAR];
    }
}
