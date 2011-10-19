/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "chol_solve_ip.h"


/****************************************
 * FUNCTIONS 
 ****************************************/

//==============================================
// constructors / destructors

/**
 * @brief Constructor
 *
 * @param[in] N size of the preview window.
 */
chol_solve_ip::chol_solve_ip (const int N) : ecL(N)
{
    w = new double[N*NUM_STATE_VAR];
}


chol_solve_ip::~chol_solve_ip()
{
    if (w != NULL)
        delete w;
}
//==============================================


/**
 * @brief Determines feasible descent direction.
 *
 * @param[in] ppar          parameters.
 * @param[in] i2hess_grad   negated inverted hessian * g.
 * @param[in] i2hess        diagonal elements of inverted hessian.
 * @param[in] x             initial guess.
 * @param[out] dx           feasible descent direction, must be allocated.
 */
void chol_solve_ip::solve(
        const problem_parameters* ppar, 
        const double *i2hess_grad,
        const double *i2hess,
        const double *x, 
        double *dx)
{
    double *s_w = w;
    int i,j;
    double i2Q[2] = {ppar->i2Q[1], ppar->i2Q[2]};


    // generate L
    ecL.form (ppar, i2hess);

    // obtain s = E * x;
    E.form_Ex (ppar, i2hess_grad, s_w);

    // obtain w
    ecL.solve_forward(ppar->N, s_w);
    ecL.solve_backward(ppar->N, s_w);

    // E' * w
    E.form_ETx (ppar, s_w, dx);

    
    // dx = -iH*(grad + E'*w)
    //
    // dx   -( -i2hess_grad  + inv(H) *   dx   ) 
    for (i = 0,j = 0; i < ppar->N*2; i++)
    {
        // dx for state variables
        dx[j] = i2hess_grad[j] - i2hess[i] * dx[j]; 
        j++;
        dx[j] = i2hess_grad[j] - i2Q[0] * dx[j]; 
        j++;
        dx[j] = i2hess_grad[j] - i2Q[1] * dx[j]; 
        j++;
    }
    for (i = ppar->N*NUM_STATE_VAR; i < ppar->N*NUM_VAR; i++)
    {
        // dx for control variables
        dx[i] = i2hess_grad[i] - ppar->i2P * dx[i];
    }
}
