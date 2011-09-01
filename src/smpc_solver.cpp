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


void smpc_solver::get_ZMP_CoM (double *ZMP_x, double *ZMP_y, double *CoM_x, double *CoM_y)
{
    qpas_solver->get_ZMP_CoM (ZMP_x, ZMP_y, CoM_x, CoM_y);
}
