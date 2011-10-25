/** 
 * @file
 * @brief The interface class, essentially, a wrapper around #qp_as.
 *
 * @author Alexander Sherikov
 * @date 02.09.2011 00:49:30 MSD
 */



/****************************************
 * INCLUDES 
 ****************************************/

#include "qp_solver.h"
#include "qp_as.h"
#include "qp_ip.h"
#include "smpc_solver.h"
#include "state_handling.h"


/****************************************
 * FUNCTIONS 
 ****************************************/



smpc_solver::smpc_solver (
                const int N,
                const solver_type sol_type,
                const double Alpha, const double Beta, const double Gamma,
                const double regularization, const double tol)
{
    if (sol_type == SMPC_AS)
    {
        qp_sol = new qp_as (N, Alpha, Beta, Gamma, regularization, tol);
    }
    else if (sol_type == SMPC_IP)
    {
        qp_sol = new qp_ip (N, Alpha, Beta, Gamma, regularization, tol);
    }
}


smpc_solver::~smpc_solver()
{
    if (qp_sol != NULL)
    {
        delete qp_sol;
    }
}



void smpc_solver::set_ip_parameters (
        const double t,
        const double mu,
        const double bs_alpha,
        const double bs_beta,
        const int max_iter)
{
    qp_ip * qpip_solver = dynamic_cast<qp_ip*> (qp_sol);
    if (qpip_solver != NULL)
    {
        qpip_solver->set_ip_parameters (t, mu, bs_alpha, bs_beta, max_iter);
    }
}



void smpc_solver::set_parameters(
        const double* T, const double* h,
        const double* angle,
        const double* zref_x, const double* zref_y,
        const double* lb, const double* ub)
{
    if (qp_sol != NULL)
    {
        qp_sol->set_parameters(T, h, angle, zref_x, zref_y, lb, ub);
    }
}



void smpc_solver::form_init_fp (
        const double *x_coord,
        const double *y_coord,
        const double *X_tilde,
        double* X)
{
    if (qp_sol != NULL)
    {
        qp_sol->form_init_fp (x_coord, y_coord, X_tilde, X);
    }
}



int smpc_solver::solve()
{
    if (qp_sol != NULL)
    {
        return (qp_sol->solve ());
    }
    return (-1);
}


void smpc_solver::get_next_state_tilde (double *state)
{
    if (qp_sol != NULL)
    {
        state_handling::get_next_state_tilde (qp_sol, qp_sol->X, state);
    }
}



void smpc_solver::get_next_state (double *state)
{
    if (qp_sol != NULL)
    {
        state_handling::get_next_state (qp_sol, qp_sol->X, state);
    }
}
