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



smpc::solver::solver (
                const int N,
                const double Alpha, const double Beta, const double Gamma,
                const double regularization, const double tol)
{
    qp_sol = new qp_as (N, Alpha, Beta, Gamma, regularization, tol);
}


smpc::solver::solver (
                const int N,
                const int max_iter,
                const double Alpha, const double Beta, const double Gamma,
                const double regularization, 
                const double tol, const double tol_out,
                const double t,
                const double mu,
                const double bs_alpha, const double bs_beta)
{
    qp_ip * qpip_solver = new qp_ip (N, Alpha, Beta, Gamma, regularization, tol);
    qpip_solver->set_ip_parameters (t, mu, bs_alpha, bs_beta, max_iter, tol_out);
    qp_sol = qpip_solver;
}


smpc::solver::~solver()
{
    if (qp_sol != NULL)
    {
        delete qp_sol;
    }
}



void smpc::solver::set_parameters(
        const double* T, const double* h, const double h_initial,
        const double* angle,
        const double* zref_x, const double* zref_y,
        const double* lb, const double* ub)
{
    if (qp_sol != NULL)
    {
        qp_sol->set_parameters(T, h, h_initial, angle, zref_x, zref_y, lb, ub);
    }
}



void smpc::solver::form_init_fp (
        const double *x_coord,
        const double *y_coord,
        const state_orig &init_state,
        double* X)
{
    if (qp_sol != NULL)
    {
        qp_sol->form_init_fp (x_coord, y_coord, init_state.state_vector, X);
    }
}



int smpc::solver::solve()
{
    if (qp_sol != NULL)
    {
        return (qp_sol->solve ());
    }
    return (-1);
}


//************************************************************


smpc::state::state()
{
    state_vector[0] = 0.0;
    state_vector[1] = 0.0;
    state_vector[2] = 0.0;
    state_vector[3] = 0.0;
    state_vector[4] = 0.0;
    state_vector[5] = 0.0;
}


void smpc::state::set (double x_, double y_)
{
    set (x_, 0.0, 0.0, y_, 0.0, 0.0);
}


void smpc::state::set (
        double x_, double vx_, double ax_,
        double y_, double vy_, double ay_)
{
    state_vector[0] = x_;
    state_vector[1] = vx_;
    state_vector[2] = ax_;
    state_vector[3] = y_;
    state_vector[4] = vy_;
    state_vector[5] = ay_;
}


//************************************************************


void smpc::state_tilde::get_next_state (const smpc::solver &smpc_solver)
{
    get_state (smpc_solver, 0);
}


void smpc::state_tilde::get_state (const smpc::solver &smpc_solver, const int ind)
{
    if (smpc_solver.qp_sol != NULL)
    {
        state_handling::get_state_tilde (
                smpc_solver.qp_sol, 
                smpc_solver.qp_sol->X, 
                ind, 
                state_vector);
    }
}


//************************************************************


void smpc::state_orig::get_next_state (const smpc::solver &smpc_solver)
{
    get_state (smpc_solver, 0);
}


void smpc::state_orig::get_state (const smpc::solver &smpc_solver, const int ind)
{
    if (smpc_solver.qp_sol != NULL)
    {
        state_handling::get_state (
                smpc_solver.qp_sol, 
                smpc_solver.qp_sol->X, 
                ind, 
                state_vector);
    }
}


//************************************************************

smpc::control::control()
{
    control_vector[0] = 0.0;
    control_vector[1] = 0.0;
}


void smpc::control::get_first_controls (const smpc::solver &smpc_solver)
{
    get_controls (smpc_solver, 0);
}


void smpc::control::get_controls (const smpc::solver &smpc_solver, const int ind)
{
    if (smpc_solver.qp_sol != NULL)
    {
        state_handling::get_controls (
                smpc_solver.qp_sol->N, 
                smpc_solver.qp_sol->X, 
                ind, 
                control_vector);
    }
}
