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

    // Matrices A and B are generated on the fly using these parameters.
#ifndef SMPC_VARIABLE_T_h
    double T = ppar->T[0];
    double T2 = T*T/2;
    double B0 = T2*T/3 - ppar->h[0]*T;
    double A6 = T2;
#endif


    const double *control = &x[ppar->N*NUM_STATE_VAR];
    double *res = result;

    for (i = 0; i < ppar->N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = ppar->angle_cos[i];
        double sinA = ppar->angle_sin[i];
#ifdef SMPC_VARIABLE_T_h
        double T = ppar->T[i];
        double T2 = T*T/2;
        double B0 = T2*T/3 - ppar->h[i]*T;
#endif

        // a pointer to 6 current elements of result

        // a pointer to 6 current state variables
        const double *xc = &x[i*NUM_STATE_VAR];



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
#ifdef SMPC_VARIABLE_T_h
            double A6 = T2 - ppar->dh[j];
#endif
            xc = &x[j*NUM_STATE_VAR];

            cosA = ppar->angle_cos[j];
            sinA = ppar->angle_sin[j];

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
 * @param[in] ppar parameters.
 * @param[in] x vector x (#NUM_STATE_VAR * N).
 * @param[out] result vector E' * nu (#NUM_VAR * N)
 */
void matrix_E::form_ETx (const problem_parameters* ppar, const double *x, double *result)
{
    int i;

    // Matrices A and B are generated on the fly using these parameters.
#ifndef SMPC_VARIABLE_T_h
    double T = ppar->T[0];
    double T2 = T*T/2;
    double B0 = T2*T/3 - ppar->h[0]*T;
    double A3 = T;
    double A6 = T2;
#endif


    double *res = result;
    double *control_res = &result[ppar->N*NUM_STATE_VAR];

    for (i = 0; i < ppar->N; i++)
    {
        // cos and sin of the current angle to form R
        double cosA = ppar->angle_cos[i];
        double sinA = ppar->angle_sin[i];


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
            double A3 = ppar->T[i+1];
            double A6 = A3*A3/2 - ppar->dh[i];
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
        double T = ppar->T[i];
        double T2 = T*T/2;
        double B0 = T2*T/3 - ppar->h[i]*T;
#endif

        // result = B' * x
        control_res[0] = B0 * xc[0] + T2 * xc[1] + T * xc[2];
        control_res[1] = B0 * xc[3] + T2 * xc[4] + T * xc[5];


        res = &res[NUM_STATE_VAR];
        control_res = &control_res[NUM_CONTROL_VAR];
    }
}
