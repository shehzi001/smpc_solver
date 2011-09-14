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
                int N,
                double Alpha, double Beta, double Gamma,
                double regularization, double tol) : 
    qpas_solver (new qp_as (N, Alpha, Beta, Gamma, regularization, tol)) {}


smpc_solver::~smpc_solver()
{
    delete qpas_solver;
}


void smpc_solver::init(
        double* T, double* h,
        double* angle,
        double* zref_x, double* zref_y,
        double* lb, double* ub,
        double* X)
{
    qpas_solver->init(T, h, angle, zref_x, zref_y, lb, ub, X);
}


int smpc_solver::solve()
{
    return (qpas_solver->solve ());
}


void smpc_solver::get_next_state_tilde (double *state)
{
    state_handling::get_next_state_tilde (state, qpas_solver->X, qpas_solver->chol_param);
}

void smpc_solver::get_next_state (double *state)
{
    state_handling::get_next_state (state, qpas_solver->X, qpas_solver->chol_param);
}
