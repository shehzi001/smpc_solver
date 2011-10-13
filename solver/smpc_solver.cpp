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

#include "qp_as.h"
#include "smpc_solver.h"
#include "state_handling.h"


/****************************************
 * FUNCTIONS 
 ****************************************/



smpc_solver::smpc_solver (
                const int N,
                const double Alpha, const double Beta, const double Gamma,
                const double regularization, const double tol) : 
    qpas_solver (new qp_as (N, Alpha, Beta, Gamma, regularization, tol)) {}


smpc_solver::~smpc_solver()
{
    delete qpas_solver;
}


void smpc_solver::set_parameters(
        const double* T, const double* h,
        const double* angle,
        const double* zref_x, const double* zref_y,
        const double* lb, const double* ub)
{
    qpas_solver->set_parameters(T, h, angle, zref_x, zref_y, lb, ub);
}

void smpc_solver::form_init_fp (
        const double *x_coord,
        const double *y_coord,
        const double *X_tilde,
        double* X)
{
    qpas_solver->form_init_fp (x_coord, y_coord, X_tilde, X);
}



int smpc_solver::solve()
{
    return (qpas_solver->solve ());
}


void smpc_solver::get_next_state_tilde (double *state)
{
    state_handling::get_next_state_tilde (qpas_solver->chol_param, qpas_solver->X, state);
}

void smpc_solver::get_next_state (double *state)
{
    state_handling::get_next_state (qpas_solver->chol_param, qpas_solver->X, state);
}
