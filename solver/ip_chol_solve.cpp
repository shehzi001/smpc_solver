/** 
 * @file
 * @author Alexander Sherikov
 * @date 19.07.2011 22:30:13 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "ip_chol_solve.h"


/****************************************
 * FUNCTIONS 
 ****************************************/
namespace IP
{
    //==============================================
    // constructors / destructors

    /**
     * @brief Constructor
     *
     * @param[in] N size of the preview window.
     */
    chol_solve::chol_solve (const int N) : ecL(N)
    {
        w = new double[N*SMPC_NUM_STATE_VAR];
    }


    chol_solve::~chol_solve()
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
    void chol_solve::solve(
            const problem_parameters& ppar, 
            const double *i2hess_grad,
            const double *i2hess,
            const double *x, 
            double *dx)
    {
        double *s_w = w;
        int i,j;


        // generate L
        ecL.form (ppar, i2hess);

        // obtain s = E * x;
        E.form_Ex (ppar, i2hess_grad, s_w);

        // obtain w
        ecL.solve_forward(ppar.N, s_w);
        ecL.solve_backward(ppar.N, s_w);

        // E' * w
        E.form_ETx (ppar, s_w, dx);

        
        // dx = -iH*(grad + E'*w)
        //
        // dx   -( -i2hess_grad  + inv(H) *   dx   ) 
        const double i2H[3] = {ppar.i2Q[1], ppar.i2Q[2], ppar.i2P};
        for (i = 0, j = 0; i < ppar.N*2; i += 2, j += SMPC_NUM_STATE_VAR)
        {
            // dx for state variables
            dx[j]   = i2hess_grad[j]   - i2hess[i] * dx[j]; 
            dx[j+1] = i2hess_grad[j+1] - i2H[0] * dx[j+1]; 
            dx[j+2] = i2hess_grad[j+2] - i2H[1] * dx[j+2]; 
            dx[j+3] = i2hess_grad[j+3] - i2hess[i+1] * dx[j+3]; 
            dx[j+4] = i2hess_grad[j+4] - i2H[0] * dx[j+4]; 
            dx[j+5] = i2hess_grad[j+5] - i2H[1] * dx[j+5]; 
        }
        for (i = ppar.N*SMPC_NUM_STATE_VAR; i < ppar.N*SMPC_NUM_VAR; i += 2)
        {
            // dx for control variables
            dx[i]   = i2hess_grad[i]   - i2H[2] * dx[i];
            dx[i+1] = i2hess_grad[i+1] - i2H[2] * dx[i+1];
        }
    }
}
